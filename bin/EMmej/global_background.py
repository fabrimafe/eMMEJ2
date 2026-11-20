#!/usr/bin/env python3
"""
global_background.py - genome-wide controls for the EM (used by src/EM_12.R)

For every kind of event present in a detector output, the same event is placed at N random
positions across the genome and run through the detector. The EM uses these rows to estimate
how often each mechanism has an opportunity at a random genomic site (the genome-wide rate),
with its own lambda/beta, because the raw motif columns are kept.

An "event key" identifies what is placed at the random sites:
    DEL:<L>            deletion of L bp                    (one set of rows per length)
    INS:<SEQ>          insertion of SEQ                    (one set per inserted sequence)
    SUB:<L>:<SEQ>      L bp replaced by SEQ                (one set per length + new sequence)
Keys are taken from the alignment rows (alignment_ID != 0) of the detector output.

The random sites are run through the real detector (default: RMdetector3_2.py next to this file),
so the rules are exactly the ones used for the observed variants. Use the SAME -w / -mhl / -ex
that were used for the observed variants.

Typical use
-----------
1) no look-up table (everything computed on the fly):
   python3 global_background.py -i det.tsv -r genome.fa -o det.global.tsv -w 150 -mhl 2

2) with a look-up table (missing keys are computed and added to the table):
   python3 global_background.py -i det.tsv -r genome.fa -o det.global.tsv --table genome.table.tsv.gz

3) pre-compute a table once per genome (all deletion lengths, all short insertions):
   python3 global_background.py --precompute --max_del_len 50 --max_ins_len 3 \
           -r genome.fa --table genome.table.tsv.gz -t 8

Then:  Rscript src/EM_12.R -i det.tsv --global det.global.tsv --local_weight 0.2 ...
"""
import argparse
import itertools
import json
import os
import random
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
from pysam import FastaFile

REPO = "/home/fabio/git/eMMEJ2"
HERE = REPO + "/bin/EMmej"
sys.path.append('scripts')
sys.path.append("/home/fabio/git/eMMEJ2/scripts")
from sliding_vcf import sliding_vcf  # same trimming used by DSB_background_generator

FLANK = 1000  # the detector discards variants with POS <= 1000 or non-ACGT within +-1000 bp


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", help="detector output (tsv) of the observed variants; keys are read from it")
    p.add_argument("-r", "--ref", required=True, help="reference genome (indexed fasta)")
    p.add_argument("-o", "--out", help="output: genome-wide control rows for the keys in --input")
    p.add_argument("-n", "--n_sites", type=int, default=1000,
                   help="number of random genomic sites per key (default 1000)")
    p.add_argument("--seed", type=int, default=1, help="seed for the random sites (default 1)")
    p.add_argument("--table", help="look-up table (tsv or tsv.gz). Keys found here are not recomputed; "
                                   "new keys are computed and appended. Created if missing.")
    p.add_argument("--precompute", action="store_true",
                   help="build/extend --table for all deletion lengths and short insertions (no --input needed)")
    p.add_argument("--max_del_len", type=int, default=50, help="with --precompute: deletion lengths 1..N (default 50)")
    p.add_argument("--max_ins_len", type=int, default=3,
                   help="with --precompute: all inserted sequences up to this length (default 3 = 84 sequences)")
    p.add_argument("-t", "--threads", type=int, default=1, help="detector processes to run in parallel")
    p.add_argument("--batch", type=int, default=5000, help="virtual variants per detector call (default 5000)")
    p.add_argument("--detector", default=os.path.join(HERE, "RMdetector3_2.py"), help="detector script to use")
    # passed through to the detector: must match the run on the observed variants
    p.add_argument("-w", "--windowsize", type=int, default=150)
    p.add_argument("-mhl", "--MH_lengths", default="2")
    p.add_argument("-ex", "--extension", action="store_true", default=False)
    a = p.parse_args()
    if a.precompute and not a.table:
        p.error("--precompute needs --table")
    if not a.precompute and not (a.input and a.out):
        p.error("give -i and -o (or use --precompute)")
    return a


def row_key(indel_type, anc, der):
    """Event key of one detector row. The same rule is implemented in EM_12.R (event_key)."""
    anc, der = str(anc).upper(), str(der).upper()
    if indel_type == "DEL":
        return f"DEL:{len(anc) - 1}"
    if indel_type == "INS":
        return f"INS:{der[1:]}"
    return f"SUB:{len(anc)}:{der}"


def keys_from_detector_output(path):
    d = pd.read_csv(path, sep="\t", usecols=["ANC", "DER", "alignment_ID", "indel_type"],
                    dtype={"ANC": str, "DER": str, "indel_type": str}, low_memory=False)
    d = d[d["alignment_ID"] != 0].drop_duplicates(subset=["ANC", "DER", "indel_type"])
    return sorted({row_key(t, a, e) for t, a, e in zip(d["indel_type"], d["ANC"], d["DER"])})


def sample_sites(fa, n, seed, margin):
    """n random sites, chromosomes weighted by length, with an ACGT-only +-1 kb context."""
    rng = random.Random(seed)
    chroms = [c for c, l in zip(fa.references, fa.lengths) if l > 2 * (FLANK + margin) + 10]
    lengths = [fa.get_reference_length(c) for c in chroms]
    if not chroms:
        sys.exit("No chromosome is long enough to place random sites.")
    sites, tries = [], 0
    while len(sites) < n:
        tries += 1
        if tries > 200 * n + 10000:
            sys.exit("Could not find enough ACGT-only sites: is the genome mostly N / masked?")
        c = rng.choices(chroms, weights=lengths)[0]
        L = fa.get_reference_length(c)
        pos = rng.randint(FLANK + 2, L - FLANK - margin - 2)
        ctx = fa.fetch(c, pos - FLANK, pos + FLANK + margin).upper()
        if set(ctx) <= set("ACGT"):
            sites.append((c, pos))
    return sites


def virtual_event(fa, chrom, pos, key):
    """The event `key` placed at (chrom, pos), built as DSB_background_generator does. -> (pos, anc, der)"""
    parts = key.split(":")
    if parts[0] == "DEL":
        anc = fa.fetch(chrom, pos - 1, pos + int(parts[1])).upper()
        return pos, anc, anc[0]
    if parts[0] == "INS":
        anc = fa.fetch(chrom, pos - 1, pos).upper()
        return pos, anc, anc + parts[1]
    anc = fa.fetch(chrom, pos, pos + int(parts[1])).upper()
    anc, der, _, n_prefix, _ = sliding_vcf(anc, parts[2])
    if anc == der:
        return None  # the substitution would not change the sequence here
    return pos + n_prefix, anc, der


def run_detector(args, vcf_path, out_path):
    cmd = [sys.executable, args.detector, "-v", vcf_path, "-r", args.ref, "-o", out_path,
           "-w", str(args.windowsize), "-mhl", args.MH_lengths]
    if args.extension:
        cmd.append("-ex")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [os.path.join(REPO, "src"), os.path.join(REPO, "scripts"), env.get("PYTHONPATH", "")])
    r = subprocess.run(cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0 or not os.path.exists(out_path):
        sys.exit("Detector failed on the random sites:\n" + r.stderr[-3000:])
    return pd.read_csv(out_path, sep="\t", dtype={"CHR": str, "ANC": str, "DER": str}, low_memory=False)


def compute_keys(args, fa, keys):
    """Run the detector on N random sites for each key. Returns detector rows + column global_key."""
    max_len = max([int(k.split(":")[1]) for k in keys if not k.startswith("INS")] + [1])
    sites = sample_sites(fa, args.n_sites, args.seed, max_len)

    events = []  # (chrom, pos, anc, der, key)
    for key in keys:
        for chrom, pos in sites:
            ev = virtual_event(fa, chrom, pos, key)
            if ev is not None:
                events.append((chrom, ev[0], ev[1], ev[2], key))
    ev = pd.DataFrame(events, columns=["CHR", "POS", "ANC", "DER", "global_key"])
    uniq = ev.drop_duplicates(subset=["CHR", "POS", "ANC", "DER"]).reset_index(drop=True)

    tmp = tempfile.mkdtemp(prefix="global_bg_")
    jobs = []
    for b, start in enumerate(range(0, len(uniq), args.batch)):
        chunk = uniq.iloc[start:start + args.batch]
        vcf = os.path.join(tmp, f"batch{b}.vcf")
        out = pd.DataFrame({"CHR": chunk["CHR"], "POS": chunk["POS"],
                            "id": chunk["CHR"] + "_" + chunk["POS"].astype(str),
                            "REF": chunk["ANC"], "ALT": chunk["DER"], "ANC": 0})
        with open(vcf, "w") as fh:
            fh.write("#CHR\tPOS\tvariant_id\tREF\tALT\tANC\n")
            out.to_csv(fh, sep="\t", header=False, index=False)
        jobs.append((vcf, os.path.join(tmp, f"batch{b}.det.tsv")))

    print(f"[global_background] {len(keys)} keys x {len(sites)} sites = {len(uniq)} virtual variants "
          f"in {len(jobs)} detector run(s)", file=sys.stderr)
    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
        dets = list(pool.map(lambda j: run_detector(args, j[0], j[1]), jobs))
    det = pd.concat(dets, ignore_index=True)
    for vcf, out in jobs:
        for f in (vcf, out):
            if os.path.exists(f):
                os.remove(f)
    os.rmdir(tmp)

    det = det.drop_duplicates(subset=["CHR", "POS", "ANC", "DER"])
    res = ev.merge(det, on=["CHR", "POS", "ANC", "DER"], how="inner")
    res["alignment_ID"] = -1  # marks genome-wide control rows
    lost = len(ev) - len(res)
    if lost > 0.05 * len(ev):
        print(f"[global_background] WARNING: {lost} of {len(ev)} virtual variants were dropped by the detector",
              file=sys.stderr)
    return res


def table_meta(args):
    return {"ref": os.path.basename(args.ref), "n_sites": args.n_sites, "seed": args.seed,
            "windowsize": args.windowsize, "MH_lengths": args.MH_lengths, "extension": args.extension}


def load_table(args):
    if not args.table or not os.path.exists(args.table):
        return None
    meta_path = args.table + ".meta.json"
    if os.path.exists(meta_path):
        old, new = json.load(open(meta_path)), table_meta(args)
        diff = {k: (old.get(k), new[k]) for k in new if old.get(k) != new[k]}
        if diff:
            sys.exit(f"The table {args.table} was built with different settings (table, now): {diff}\n"
                     "Use the same settings or a different --table file.")
    return pd.read_csv(args.table, sep="\t", dtype={"CHR": str, "ANC": str, "DER": str, "global_key": str},
                       low_memory=False)


def save_table(args, table):
    table.to_csv(args.table, sep="\t", index=False)
    json.dump(table_meta(args), open(args.table + ".meta.json", "w"), indent=1)


def main():
    args = parse_args()
    fa = FastaFile(args.ref)

    if args.precompute:
        keys = [f"DEL:{L}" for L in range(1, args.max_del_len + 1)]
        for L in range(1, args.max_ins_len + 1):
            keys += ["INS:" + "".join(s) for s in itertools.product("ACGT", repeat=L)]
    else:
        keys = keys_from_detector_output(args.input)

    table = load_table(args)
    have = set(table["global_key"].unique()) if table is not None else set()
    missing = [k for k in keys if k not in have]
    print(f"[global_background] keys needed: {len(keys)}; in table: {len(keys) - len(missing)}; "
          f"to compute: {len(missing)}", file=sys.stderr)

    new = compute_keys(args, fa, missing) if missing else None
    parts = [x for x in (table, new) if x is not None]
    full = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    if args.table and new is not None:
        save_table(args, full)
        print(f"[global_background] table updated: {args.table} ({full['global_key'].nunique()} keys)",
              file=sys.stderr)
    if args.out:
        full[full["global_key"].isin(set(keys))].to_csv(args.out, sep="\t", index=False)
        print(f"[global_background] written {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()


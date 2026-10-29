import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--sim", required=True)    # file della simulazione
p.add_argument("--det", required=True)    # output della detection
p.add_argument("--out", required=True)    # prefisso di output
a = p.parse_args()

# sep=None + engine="python" rileva da solo tab/virgola
sim = pd.read_csv(a.sim, sep=None, engine="python")
det = pd.read_csv(a.det, sep=None, engine="python")

keep_ids = set(det["variant_id"].astype(str))
mask = sim["variant_id"].astype(str).isin(keep_ids)

kept = sim[mask]
removed = sim[~mask]

kept.to_csv(f"{a.out}_kept.tsv", sep="\t", index=False)
removed.to_csv(f"{a.out}_removed.tsv", sep="\t", index=False)

print(f"{a.sim}: total={len(sim)} kept={len(kept)} removed={len(removed)}")

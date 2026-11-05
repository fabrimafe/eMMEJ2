#!/usr/bin/env python3
import argparse
import os
import random

"""
prova 4: con pol slip
GGATTCACGTCGGAACGTTCTGTACGATACGATCCAGTCCTCCGACTCCTATA

"""


# --- PARAMETRI ---
N_CROMOSOMI = 1
LUNGHEZZA_UNITA = 52          # bp dell'unità ripetuta (per quelle generate)
RIPETIZIONI = 10000000             # ripetizioni di default per ogni cromosoma
SEED = 42                      # stesso seed = stesso genoma ogni volta
LARGHEZZA_RIGA = 60
OUTPUT_DEFAULT = "genoma_sintetico.fa"

# Sequenze scritte a mano (facoltative): chiave = nome cromosoma
# Quelle che non specifichi qui vengono generate casualmente.
SEQUENZE_MANUALI = {
    "chr1": "GGATTCACGTCGGAACGTTCTGTACGATACGATCCAGTCCTCCGACTCCTATA",
    # "chr2": "GGCCTTAA...",
}

# Ripetizioni personalizzate per cromosoma (facoltativo)
RIPETIZIONI_CUSTOM = {
    # "chr1": 5000,
}
# -----------------

def sequenza_casuale(lunghezza, rng):
    return "".join(rng.choice("ACGT") for _ in range(lunghezza))

def parse_args():
    parser = argparse.ArgumentParser(
        description="Genera un genoma sintetico in formato FASTA."
    )
    parser.add_argument(
        "-o", "--output",
        default=OUTPUT_DEFAULT,
        help=f"percorso del file FASTA di output (default: {OUTPUT_DEFAULT})",
    )
    return parser.parse_args()

def main():
    args = parse_args()

    # crea la cartella di destinazione se non esiste
    cartella = os.path.dirname(args.output)
    if cartella:
        os.makedirs(cartella, exist_ok=True)

    rng = random.Random(SEED)
    with open(args.output, "w") as f:
        for i in range(1, N_CROMOSOMI + 1):
            nome = f"chr{i}"
            unita = SEQUENZE_MANUALI.get(nome) or sequenza_casuale(LUNGHEZZA_UNITA, rng)
            n = RIPETIZIONI_CUSTOM.get(nome, RIPETIZIONI)
            genoma = unita.upper() * n

            f.write(f">{nome}\n")
            for j in range(0, len(genoma), LARGHEZZA_RIGA):
                f.write(genoma[j:j + LARGHEZZA_RIGA] + "\n")

            print(f"{nome}: unità {len(unita)} bp x {n} = {len(genoma)} bp")

    print(f"\nScritto: {args.output}")

if __name__ == "__main__":
    main()

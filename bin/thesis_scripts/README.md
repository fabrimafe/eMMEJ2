# "A new tool for estimating the contribution of Double Strand Break DNA repair mechanisms in shaping short structural variation across genomes" thesis scripts
## This folder contains the essential scripts that I wrote for this thesis work in Prof. Avraham Levy lab, at the department of plant and environmental science, the Weizmann institute of science.

### As this work is yet to be published, this repository is privet.
### To operate those scripts, one needs to first recreate the miniconda environment as follows:
```
conda env create --name recoveredenv --file guy_mmej_env.yml
conda activate guy_mmej_env
```
For WEXAC users:
```
module load miniconda
conda activate guy_mmej_env
```

### Other dependencies:
```
bedtools v2.26.0
vcftools v0.1.15
samtools v1.9
bcftools v1.9
```

### The general structure of this folder:

```
├── EMmej -> All EMmej related scripts.
├── genomic_data -> Scripts that we used to analyze genomic data.
│   ├── blockbootstrapping
│   ├── Genomic_data_preperation
│   └── plots_code
├── indel_simulations -> Scripts that we used for the simulations.
├── sdmmej_and_EMmej_comparison -> Scripts that we used for the comparison of EMmej and sdmmej.
└── src -> The modules that EMmej use.

```

### A more detailed structure can be found in the README of each folder



module load miniconda
module load tabix
module load bedtools
conda activate emmej
python ./scripts/EMmej/Markov_model_operator.py --vcf /home/labs/alevy/fabrizio/workspace/guy/rerun/03302023/Clark/RM_output.tsv --ref /home/labs/alevy/fabrizio/workspace/ancestralstate/msa_dro/data/Drosophila_melanogaster.fa -p markov -o /home/labs/alevy/fabrizio/workspace/guy/rerun/03302023/Clark/MM_output.tsv
./scripts/EMmej/EM_operator.py --vcf /home/labs/alevy/fabrizio/workspace/guy/rerun/03302023/Clark/MM_output.tsv --ref /home/labs/alevy/fabrizio/workspace/ancestralstate/msa_dro/data/Drosophila_melanogaster.fa -o /home/labs/alevy/fabrizio/workspace/guy/rerun/03302023/Clark/EM_output.tsv

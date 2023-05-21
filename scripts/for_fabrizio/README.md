# 16.07.2022

# This folder contains a shot version of the whole pipline for deletions only for fabrizio to look at the EM

# This folder contains two datasets, both come from the very same sims, one dataset went through the parsevcf.py script that you wrote (20220716_realigned_20220620_sims.vcf) and the other one didn't (20220716_non_realigned_20220620_sims.vcf)


# Instruction:
# 1st load my venv 
module load miniconda
conda activate guy_mmej_env

# I guess that any venv that have pandas and numpy are good as well, but if you haveing problems due to dependencies, use /home/labs/alevy/guyta/guy_master_project/guy_mmej_env.yml to to see my venv details

# To operate the script, run:
# flags:
# -v : the vcf
# -mh : the proportions of MMEJ simulations you want to check
# -w : window size of the context, keep this to 1500
# -r : reference genome (/home/labs/alevy/guyta/guy_master_project/data/indels_simulations/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta )

python3 20220715_whole_pipline_in_one_script.py -v 20220716_realigned_20220620_sims.vcf -mh 0.9 -w 1500 -r /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta

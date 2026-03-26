# This folder contains all the scripts for the block-bootstrapping
### see "EMmej analysis of genomic data"

```
├── batch_data_and_submit_EMmej_to_cluster.py -> A WEXAC specific script that devide
|              the data into blocks (using whole_genome_batching_by_blocks.py)
|               and submit a job of EMmej over each block seperatly.
├── EMmej_blockbootstrap.py -> A script that performs the block-bootstrapping itself.
├── EMmej_blockbootstrap_summarizer.py -> A script that summarizes the output of 
|               EMmej_blockbootstrap.py into a plotable tables.
├── README.md
└── whole_genome_batching_by_blocks.py -> Split data to blocks.
```
# EMmej Util scripts
## This folder contains scripts that help manipulate input data prior to EMmej operation, those steps are not part of the EMmej pipeline.

### Converting input data into the VCF-like format for EMmej analysis.
If data is in a standart VCF format, use vcf_format_converter.sh with the following arguments: <br>
<br>
available options: <br>
-h : Argumants docs  <br>
-r : indexed (using tabix) reference genome (.fa & .fai) <br>
-v : input vcf file (.vcf) <br>
-o : output path (string) <br>
-f: indication whether an AF file (vcftools --freq2 .frq format) is available (bool) <br>
<br>
In order to reformat it into the right VCF-like format.

### Split data into blocks for block-bootstraping <br>
Use whole_genome_block_bootrapping.py in order to split the data into blocks based on chromosome and positions, this step is recomended when one wants to account for the effect of blocks in the genome on the DNA DSB repair mechanisms.

arguments: <br>
	-h, --help            show this help message and exit <br>
	-v VCF, --vcf VCF     path to input vcf (string) <br>
	-o OUTPUTFILE, --outputfile OUTPUTFILE, path to output file (string)<br>
	-nb NBLOCKS, --Nblocks NBLOCKS, number of blocks to devide per chromosome (int)<br>
	-ns NSAMP, --Nsamp NSAMP, number of samples per block (int)<br>
	-Nb NBOOT, --Nboot NBOOT, number of iteration for bootstraping (int)<br>
	-anc ANCESTRAL, --ancestral ANCESTRAL, indicate whether ancestral state column is available (int)<br>

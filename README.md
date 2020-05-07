# hvAssembler

This tool searches for reads containing hypervariable genes genotype specific kmers and 
store them in different fastq file. Such files are then assembled by using the OLC assembler
cap3. The command to launch the script is the following

python3 hvAssembler.py read_1.fastq read_2.fastq condaDir outputFolderName

the condaDir is the conda enviroment directory

 In the output folder the following will be reported:
 
 kmerSpecificReads 
 A folder containing the fastq files for the reads that have been found for each gene/genotype 
 combination
 
 kmerSpecificScaffolds
 A folder containing the assembled scaffolds for the reads present in the kmerSpecificReads
 folder
 
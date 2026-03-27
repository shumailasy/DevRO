###############
# DevRO: Variant Caller: Dup Script 
#@version DevRO_v4.0.55
##
##
##########################################In this step, we scan the genome in widows of 1kbp size and search for duplicated, singleton, softclip and inverted reads using the bitwise flags in alignment files as the method described in SAMTOOLs. For duplications, we store all those loci where atleast 1 duplicated read is observed in any population, and record the information of deviant/abnormal/normal read counts for that loci in each population, median positions is calculated using the mates (abnormal reads) *details in algorithm section. This step just gives us the raw calls for loci with reads counts in each category for all the populations in question.##We use following filters at this step:#1.	mapping quality for dup read >=10.  (same criteria applied for each SV type in question for which the script is run)#2.	Each analysis is done separately for duplications, inversions, deletions and insertions. #a.	This means if we are running script for duplications scan and we come across inverted loci having no duplication reads then this loci will not be reported in output of duplications call but will be reported in inversion calls.#b.	But if we come across inverted reads where we have duplicated reads also, the loci will be stored with the information of reads counts for duplications and inversion. Vice versa in inversion script.## Results Variant Caller Duplications: 

############ ############################
##### run : perl VariantCaller_dup.pl genomeFile.txt prefix_Output PATHBAMFiles.config
###
###
##### The script takes Input file name (genome region file) and BAMs for each population and Output file prefix
## genome region file format: 3 columns Tab Seperated, #chromosome	chrSize	StartCoordinate

#### (1) you can give complete chromosome coordinates StartCoordinate as 0 and chrSize as Maxsize of chromosome // but its slow
#### (2) you can give coordinates in split windows e.g. 0-32K , next 32K-64K, ect.. 

### that how you split chr…
## window file on each chr size=32000000
########### bedtools makewindows -g hg.len -w 32000000 > hg.32000000bp.windows.bed

##split jobs
############# split -d -l 2 hg.32000000bp.windows.bed chr_
 
##### Next you need config file containing group info and the each population BAM file with path
##### 1	PATH/sample1.dedup.bam
##### 1 PATH/sample2.dedup.bam
##### 2	PATH/sample3.dedup.bam
##### 2 PATH/sample4.dedup.bam

### here column 1 tells we have 2 groups and column 2 tells the path and bam.

##########################################

######### OUTPUT format
#################################

1.bin
2.coordinates 1k window scan
3.size if raw SV
4.coordinates of mates window
5.median_mate_postion
6.coordinates_raw_SV
7.median position of singletons on forward strand
8.median position of singletons on forward strand
9.Hash seperated inverted reads counts in each population on forward strand
10.Hash seperated inverted mate counts in each population on forward strand
11.Hash seperated inverted reads counts in each population on reverse strand
12.Hash seperated inverted mate counts in each population on reverse strand 
13-16 Softclips same as 9-12
17-20 Singletons same as 9-12
21 Duplication reads counts in each population hash seperated
22 Duplication mate counts in each population hash seperated
23 Total reads counts in each population hash seperated
24 Total mate counts in each population hash seperated

###################### This gives raw calls, which are further parsed by parser.

























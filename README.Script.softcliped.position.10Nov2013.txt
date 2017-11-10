##This Script needs region file and config file as input
### region file format:
chr start end sizeChromosomeLimit LocusID ## tab seperated 5 columns.
### config file format:
group SampleName ## the sample name is followed by exact path of bam..

Running command:
perl Script.softclips_position.10Nov2013.pl region.bed MyTest config.file.txt

## EXAMPLE files provided in the directory including outfileFile

#### OUTPUT file format "AllpopSoftclip_Locus-ID.bed" ## the file will have clipped counts for poplution1 in region.. then second population .. third ... so on... in the order given in configFile.
1.chromosome	
2.SoftClippedPosition        
3.SoftClippedPosition        
4.softclipCount       
5.~/PATH/sample.bam
6.Locus-ID

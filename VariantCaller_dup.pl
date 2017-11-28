###############
# DevRO: Variant Caller: Dup Script 
##
##
## @copyright_shumaila.sayyab@gmail.com GNU @version DevRO_v4.0.55
##
##
#########################################
#In this step, we scan the genome in widows of 1kbp size and search for duplicated, singleton, softclip and inverted reads using the bitwise flags in alignment files as the method described in SAMTOOLs. For duplications, we store all those loci where atleast 1 duplicated read is observed in any population, and record the information of deviant/abnormal/normal read counts for that loci in each population, median positions is calculated using the mates (abnormal reads) *details in algorithm section. This step just gives us the raw calls for loci with reads counts in each category for all the populations in question.
#
#We use following filters at this step:
#1.	mapping quality for dup read >=10.  (same criteria applied for each SV type in question for which the script is run)
#2.	Each analysis is done separately for duplications, inversions, deletions and insertions. 
#a.	This means if we are running script for duplications scan and we come across inverted loci having no duplication reads then this loci will not be reported in output of duplications call but will be reported in inversion calls.
#b.	But if we come across inverted reads where we have duplicated reads also, the loci will be stored with the information of reads counts for duplications and inversion. Vice versa in inversion script.
#
# Results Varinat Caller Duplications is as follows explained in readme
## #### bin,new_reg_st,calc_distance_inv,new_reg_mate,median_dist,inv_reg,median_F,median_R,A_invy1,A_invy2,A_invy3,
## #### A_invy4,A_soft1,A_soft2,A_soft3,A_soft4,A_single1,A_single2,A_single3,A_single4,A_dup1,A_dup2,A_Total1,A_Total2
##

############ 
##### The script takes Input file name (genome region file) and BAMs for each population and Output file prefix
## genome region file format: 3 columns Tab Seperated, #chromosome	chrSize	StartCoordinate

#### (1) you can give complete chromosome coordinates StartCoordinate as 0 and chrSize as Maxsize of chromosome // but its slow
#### (2) you can give coordinates in split windows e.g. 0-32K , next 32K-64K, ect.. 
#e.g. chr	chr_size	startcoordinate
####

####
 
##### Next you need config file containing group info and the each population BAM file with path
##### 1	PATH/sample1.dedup.bam
##### 1 PATH/sample2.dedup.bam
##### 2	PATH/sample3.dedup.bam
##### 2 PATH/sample4.dedup.bam

### here column 1 tells we have 2 groups and column 2 tells the path and bam.


#use Math::Cephes qw(:explog);
use Statistics::Lite qw(:all);
use Statistics::Zscore;

$file1=shift or die "Usage: $0 input GenomeFile name: chr size start \n"; ## Region file (chr chr_size)
open (region,"$file1") or die "cannot open file, $!\n";
#@Region=<region>;
$outFile=shift or die "Usage: $1 needs OutputFile prefix \n";
$b1 = $outFile."_DupSignatures.Allpop.txt";
open (OUT, ">OUTPUT_DUP/$b1");
#open (OUT2, ">SingletonR_bam.seq");
$config=shift or die "Usage: $2 needs configFile, group and PATH to BAMfiles with name \n";
open (conFile,"$config") or die "cannot open file, $!\n";
@pop=();
while($configLine=<conFile>)
{
		chomp $configLine;
		@split_conline=split("\t",$configLine);
		
		my ($group,$BAMPATH)= @split_conline[0,1];
		push (@pop,"$BAMPATH");
#@pop=("BAMS_input/1777_sample.dedup.bam","BAMS_input/AngoraMale1_sample.dedup.bam");
}

$ReadLength_default=100; ## change this or update this...
$MappingQuality_default=1; ## change this or update this...
$meanInsert=400;
$sigma=130;
$InsertRange=$meanInsert+(3*$sigma);
$InsertRangeMin=$meanInsert-(3*$sigma);
########## Global Initialization of variables  ###########
my $pairedEndseq;
my $pair_proper_mapped;
my $read_mapped;
my $mate_mapped;
my $read_strand_plus;
my $mate_strand_plus;
my $first_read_inPair;
my $second_read_inPair;
my $Notprimary_alignment;
my $Failvendor_qualitycheck;
my $PCR_duplicate;

####################################################
while($reg=<region>)
{
		
		$flag_st=-1;
		$flag_end=-1;
		chomp $reg;
		@split_line=split("\t",$reg);
		my ($chr,$sp_st,$chr_size)= @split_line[0,1,2]; ## full bam file scanned in 1kbp windows
		#my ($chr,$chr_size)= @split_line[0,1];
		#push(@Array_pop_total_DoC_score,$total_DoC_score);
		$st_pos=$sp_st; $end_pos=$st_pos+1000;$new_reg_st=$chr."\:".$st_pos."\-".$end_pos;
		$bin=0;
	while($st_pos <= ($chr_size))
	{
		if($st_pos >0 and $end_pos <$chr_size)
		{
			$bin++;
			$i=0;
			@inversions=();
			@inversions_F=();
			@inversions_R=();
			@inversions_F_mate=();
			@inversions_R_mate=();
			
			@softclip_F=();
			@softclip_R=();
			@softclip_F_mate=();
			@softclip_R_mate=();
			@softy=();
			
			@singleton=();
			@singleton_F=();
			@singleton_R=();
			@singleton_F_mate=();
			@singleton_R_mate=();
			
			@duplication_Reads=();
			@duplication_Reads_mate=();
			
			@totalAlignedReads=();
			@totalAlignedReads_mate=();
			
			@MateDistances=();@MateDistances_dup=();
			@MateDistances_singF=();
			@MateDistances_singR=();
			$size_pop=$#pop;
			$bool=1;
			#@ReadSeq_singletons=();
			while($i <= $size_pop)
			{
				#Here we used $i because of the number of populations
			
					$count_softclip_start=0;
                	$count_softclip_end=0;
                	$flag_st=-1;
               		$flag_end=-1;
			        $total_AR=0;
                	$R_trans_start=0;
                	$F_trans_start=0;
                	$R_trans_end=0;
                	$F_trans_end=0;
                	$count_trans_end=0;
                	$flag_tr_st=-1;
                	$flag_tr_end=-1;
                	$invy_st_F=0;
                	$invy_st_R=0;
                	$invy_end_F=0;
                	$invy_end_R=0;
                	$plus_st=0;
                	$plus_end=0;
                	$minus_st=0;
                	$minus_end=0;
                	$singleton_F_st=0;
                	$singleton_F_end=0;
                	$F_read_pos=0;
                	$R_read_pos=0;
                	$singleton_R_st=0;
                	$singleton_R_end=0;
                	$dup_read_st=0;
                	$dup_read_end=0;
                	$count_softclip_plus=0;$count_softclip_minus=0;
                	###
					
					#@MateDistances_R=();
					####
                	######	
       			$call_st="samtools view $pop[$i] $new_reg_st"; 
       			@align_st=qx($call_st);
			
			### calc. library insert distribution for 1Million reads on chr1
       			if($bool >-1){
       			my $metrics=`samtools view -q 10 -f2 $pop[$i] chr1|cut -f9|head -1000000|awk '{if (\$1<0){\$1=-\$1}else{\$1=\$1} sum+=\$1; sumsq+=\$1*\$1} END {print sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}'`;
				my ($mean,$stdev)=split(/ /,$metrics);
				my ($mean,$stdev)=split(/\s/,$metrics);
				$stdev=~s/\n//;
				$InsertRange=int($mean+(3*$stdev));
				$InsertRangeMin=int($mean-(3*$stdev));
				$bool=-1;
       			}else{}
			### done calc. library insert distribution for few reads
			
			#processing for start region  ############################################
			foreach $line (@align_st)
			{
				@split_align_st=split("\t",$line);
                my ($readid,$flg,$seqid,$read_st,$MQuality,$cigar,$mateid,$mate_st,$insert,$readseq,$qual)= @split_align_st[0,1,2,3,4,5,6,7,8,9,10];
		             
		               initializeBitwiseflag();
		               my @bitFlags=calculateBitwiseflag($flg);
		               #extractflagInfo(@bitFlags);
		               ################# Extract Flag Information #########
					    if($bitFlags[0] == 1)
						{
							$pairedEndseq=1; # paired end
						}
						if($bitFlags[1] == 1)
						{
							$pair_proper_mapped=1; # pair properly mapped
						}
						if($bitFlags[2] == 0)
						{
							$read_mapped=1; #read seq is mapped
						}
						if($bitFlags[3] == 0)
						{
							$mate_mapped=1; #mate is mapped
						}
						if($bitFlags[4] == 0)
						{
							$read_strand_plus=1; #read strand is plus/forward
						}
						if($bitFlags[5] == 0)
						{
							$mate_strand_plus=1; #mate strand is plus/forward
						}
						if($bitFlags[6] == 1)
						{
							$first_read_inPair=1; # first read in pair
						}
						if($bitFlags[7] == 1)
						{
							$second_read_inPair=1; # second read in pair
						}
						if($bitFlags[8] == 1)
						{
							$Notprimary_alignment=1; # secondary alignment
						}
						if($bitFlags[9] == 0)
						{
							$Failvendor_qualitycheck=1; #fail vendor/platform quality check
						}
						if($bitFlags[10] == 1)
						{
							$PCR_duplicate=1; #PCR duplicate
					
						}
						@bitFlags=();
		              ####################################################
		               
		               #### this is special section added for the softclips
		               
		               ####################################################
		               ####################################################
		               @cigar_positions = ($cigar =~ m/(\d+)/g); ## parsing out the cigar for positions
					   @cigar_type=($cigar =~ m/(\D+)/g); ## parsing out the cigar for type (i.e. softclip, deletion, insertion,match)
					  	$softclip_position=$read_st;
					  	$prime5_softclip_position=0;
					  	$prime3_softclip_position=0;
					  	$flag_st=-1;
				
						for($sf=0;$sf<=$#cigar_type;$sf++)
						{
									if($cigar_type[$sf] eq 'S' and $cigar_positions[$sf] >=1)
									{
										 ##check if the softclip is 5'end or 3'end?? ###
					                     ## this can be done using the index position of cigar_type array
					                     ## if index zero contains S then softclip is at the 5'end
					                     ## else softclip is at the 3'end
					                     if($sf == 0)
					                     { ## 5'end
					                     	## the genomic position of read start is the softclip position also here
					                     	$prime5_softclip_position=$softclip_position;
					                     	$flag_st=1;
					                     }
					                     else
					                     { ## 3' end
					                     	## we will add the matched,deleted or inserted positions to read start to calc. 3'softclip
					                     	$prime3_softclip_position=$softclip_position;
					                     	$flag_st=2;
					                     }
					                     
									}
									elsif($cigar_type[$sf] eq 'M' and $cigar_positions[$sf] >=1)
					                {
					                	$softclip_position=$softclip_position+$cigar_positions[$sf];
					                }
					                elsif($cigar_type[$sf] eq 'I' and $cigar_positions[$sf] >=1)
					                {
					                	#$softclip_position=$softclip_position+$cigar_positions[$sf]-1;
					                }
					                elsif($cigar_type[$sf] eq 'D' and $cigar_positions[$sf] >=1)
					                {
					                	$softclip_position=$softclip_position-$cigar_positions[$sf];
					                }
					                else{}
										
						}
						
		               
		               #####################################################
		               #Here we are taking into account the duplications and inversion reads mq>=10 only and donot take into account translocations
		               ####################################################
				if($MQuality >= $MappingQuality_default and $mateid eq "=")
				{	
						 if($read_strand_plus == 1)
						 { 
						 	## It is forward strand
						 	if($flag_st > -1)
							{$flag_st=-1;
								$plus_st++;
							}else{$flag_st=-1;}
							if($mateid eq "=")
							{
								#print $mateid;
							}else{$F_trans_start++;}
							if($read_strand_plus == $mate_strand_plus and $mate_mapped==1 and $read_mapped==1 and $mateid eq "=" and $insert >=$InsertRange)
							{
		
								$mate_distance=$mate_st+$ReadLength_default-1; # 50 is added here assuming read length of 50 and going in the end of read			
								$invy_st_F++;								
							    push (@MateDistances,"$mate_distance");
								 
							}else{}
							if($mate_mapped==0 and $read_mapped==1)
							{
								
								$F_read_pos=$read_st+$ReadLength_default-1;
                                 $singleton_F_st++;
                                 push (@MateDistances_singF,"$F_read_pos");
                                
								#print OUT2 "F_st","\t",$line,"\n";
							}else{ }
							if((($flg ==145 and $insert >=$InsertRange) or ($flg ==81 and  $insert >=$InsertRange) or ($flg ==97 and $insert < $InsertRangeMin) or ($flg ==161 and $insert <$InsertRangeMin)) and ($mateid eq "=" and $mate_mapped==1 and$read_mapped==1)){
                                    $mate_distance_dup=$mate_st+$ReadLength_default-1;
                                    $dup_read_st++;
                                    push (@MateDistances_dup,"$mate_distance_dup");
                                                                

                            }else{} 
						 }
						 else
						 { 
						 	## It is reverse strand
						 	if($flag_st > -1)
							{$flag_st=-1;
								$minus_st++;
							}else{$flag_st=-1;}
							if($mateid eq "=")
							{
								#print $mateid;
							}else{$R_trans_start++;}
						 	if($read_strand_plus == $mate_strand_plus and $mate_mapped==1 and $read_mapped==1 and $mateid eq "=")
						 	{
						 		$mate_distance=$mate_st;  
						 		$invy_st_R++;								###the calculate median (+- 500 bp) to look up for all mates. 
							    push (@MateDistances,"$mate_distance");
									
						 	}else{}
						 	if($mate_mapped==0 and $read_mapped==1)
							{
								$singleton_R_st++;
								$R_read_pos=$read_st;
                                 push (@MateDistances_singR,"$R_read_pos");
                                
							}else{}
							if((($flg ==145  and $insert >=$InsertRange) or ($flg ==81  and $insert >=$InsertRange) or ($flg ==97 and $insert < $InsertRangeMin) or ($flg ==161 and $insert < $InsertRangeMin)) and ($mateid eq "=" and $mate_mapped==1 and $read_mapped==1))                                      
						 	{
						 		$mate_distance_dup=$mate_st;
								
									$dup_read_st++;#print OUT2 "R_st","\t",$line,"\n";
									push (@MateDistances_dup,"$mate_distance_dup");
								
						 		
						 	}else{}
						}
				}else{}
				### here we need to check strand of the query read.. (Do this for + and - strands seperately)
				 ### then in case of inversions.. softclips... single reads we push them to their respective vectors...
			
			} #print $call_st."\n";print $invy_st_F."#".$invy_st_R."\n";
			
						$total_translocations=$F_trans_start+$R_trans_start;
						
						$total_AR=$#align_st;
							
                        #push (@inversions,"$total_invy");
                        push (@inversions_F,"$invy_st_F");
                        push (@inversions_R,"$invy_st_R");
                        
                        push (@softclip_F,"$plus_st");
                        push (@softclip_R,"$minus_st");
                       # push (@softy,"$total_softy");
                        
                        push (@singleton_F,"$singleton_F_st");
                        push (@singleton_R,"$singleton_R_st");
                        
                        push (@duplication_Reads,"$dup_read_st");
                        
                       # push (@translocations,"$total_translocations");
                       	push (@totalAlignedReads,"$total_AR");
			@align_st=();
			$i++;
			}## end of population while loop

			#### for singleton reads median position
			$median_F=int(median(@MateDistances_singF));
			$median_R=int(median(@MateDistances_singR));
			
			$median_dist_invy=int(median(@MateDistances));
			
			$median_dist=int(median(@MateDistances_dup));
			$mate_500up=$median_dist-300; $mate_500down=$median_dist+300;

			if($median_dist == 0 or $mate_500up < 0){}else{
				$new_reg_mate=$chr."\:".$mate_500up."\-".$mate_500down;
				
				$i=0;
				while($i <= $size_pop)
				{
					#Here we used $i because of the number of populations
				
						#$count_softclip_start=0;
	                	$count_softclip_end=0;
	                	#$flag_st=-1;
	               		 $flag_end=-1;
				###### trans
	                	#$R_trans_start=0;
	                	#$F_trans_start=0;
	                	$R_trans_end=0;
	                	$F_trans_end=0;
	                	$count_trans_end=0;
	                	#$flag_tr_st=-1;
	                	$flag_tr_end=-1;
	                	#$invy_st_F=0;
	                	#$invy_st_R=0;
	                	$invy_end_F=0;
	                	$invy_end_R=0;
	                	#$plus_st=0;
	                	$total_AR_mate=0;
	                	$plus_end=0;
	                	#$minus_st=0;
	                	$minus_end=0;
	                	#$singleton_F_st=0;
	                	$singleton_F_end=0;
	                	#$singleton_R_st=0;
	                	$singleton_R_end=0;
	                	$dup_read_end=0;
	                	$count_softclip_plus=0;
	                	$count_softclip_minus=0;
	                	###
				$call_mate="samtools view $pop[$i] $new_reg_mate";@align_mate=qx($call_mate);#print $call_mate."\n\n";
						
						#processing for mate region############################################
                        foreach $lineEnd (@align_mate)
                        {
                                @split_align_end=split("\t",$lineEnd);
                                my ($readid1,$flg1,$seqid1,$read_st1,$MQuality1,$cigar1,$mateid1,$mate_st1,$insert1,$readseq1)= @split_align_end[0,1,2,3,4,5,6,7,8,9];
			                    initializeBitwiseflag();
		               my @bitFlags=calculateBitwiseflag($flg1);
		               #extractflagInfo(@bitFlags);
		               ################# Extract Flag Information #########
					    if($bitFlags[0] == 1)
						{
							$pairedEndseq=1; # paired end
						}
						if($bitFlags[1] == 1)
						{
							$pair_proper_mapped=1; # pair properly mapped
						}
						if($bitFlags[2] == 0)
						{
							$read_mapped=1; #read seq is mapped
						}
						if($bitFlags[3] == 0)
						{
							$mate_mapped=1; #mate is mapped
						}
						if($bitFlags[4] == 0)
						{
							$read_strand_plus=1; #read strand is plus/forward
						}
						if($bitFlags[5] == 0)
						{
							$mate_strand_plus=1; #mate strand is plus/forward
						}
						if($bitFlags[6] == 1)
						{
							$first_read_inPair=1; # first read in pair
						}
						if($bitFlags[7] == 1)
						{
							$second_read_inPair=1; # second read in pair
						}
						if($bitFlags[8] == 1)
						{
							$Notprimary_alignment=1; # secondary alignment
						}
						if($bitFlags[9] == 0)
						{
							$Failvendor_qualitycheck=1; #fail vendor/platform quality check
						}
						if($bitFlags[10] == 1)
						{
							$PCR_duplicate=1; #PCR duplicate
					
						}
						@bitFlags=();
				               ####################################################
		               
		               #### this is special section added for the softclips
		               
		               ####################################################
		               ####################################################
		               @cigar_positions = ($cigar1 =~ m/(\d+)/g); ## parsing out the cigar1 for positions
					   @cigar_type=($cigar1 =~ m/(\D+)/g); ## parsing out the cigar1 for type (i.e. softclip, deletion, insertion,match)
					  	$softclip_position=$read_st1;
					  	$prime5_softclip_position=0;
					  	$prime3_softclip_position=0;
					  	$flag_st=-1;
				
						for($sf=0;$sf<=$#cigar_type;$sf++)
						{
									if($cigar_type[$sf] eq 'S' and $cigar_positions[$sf] >=1)
									{
										 ##check if the softclip is 5'end or 3'end?? ###
					                     ## this can be done using the index position of cigar_type array
					                     ## if index zero contains S then softclip is at the 5'end
					                     ## else softclip is at the 3'end
					                     if($sf == 0)
					                     { ## 5'end
					                     	## the genomic position of read start is the softclip position also here
					                     	$prime5_softclip_position=$softclip_position;
					                     	$flag_st=1;
					                     }
					                     else
					                     { ## 3' end
					                     	## we will add the matched,deleted or inserted positions to read start to calc. 3'softclip
					                     	$prime3_softclip_position=$softclip_position;
					                     	$flag_st=2;
					                     }
					                     
									}
									elsif($cigar_type[$sf] eq 'M' and $cigar_positions[$sf] >=1)
					                {
					                	$softclip_position=$softclip_position+$cigar_positions[$sf];
					                }
					                elsif($cigar_type[$sf] eq 'I' and $cigar_positions[$sf] >=1)
					                {
					                	#$softclip_position=$softclip_position+$cigar_positions[$sf]-1;
					                }
					                elsif($cigar_type[$sf] eq 'D' and $cigar_positions[$sf] >=1)
					                {
					                	$softclip_position=$softclip_position-$cigar_positions[$sf];
					                }
					                else{}
										
						}
						
		               
		               #####################################################
		               #Here we are taking into account the duplications and inversion reads mq>=10 only and donot take into account translocations
		               ####################################################
				if($MQuality1 >= $MappingQuality_default and $mateid1 eq "=")
				{	
									#DeviantsMate($read_strand_plus);
											 if($read_strand_plus == 1)
											 { 
											 	## It is forward strand
											 	if($flag_end  > -1)
												{$flag_end =-1;
													$plus_end++;
												}else{$flag_end=-1;}
												if($mateid1 eq "=")
												{
													#print $mateid;
												}else{$F_trans_end++;}
												if($read_strand_plus == $mate_strand_plus and $mate_mapped==1 and $read_mapped==1 and $mateid1 eq "=" and $insert1 >=$InsertRange)
												{
													$invy_end_F++; 
												}else{}
												if($mate_mapped==0 and $read_mapped==1)
												{
													$singleton_F_end++;#print OUT2 "F_end","\t",$line,"\n";
												}else{}
												if((($flg1 == 145 and $insert1 >$InsertRange) or ($flg1 ==81 and $insert1 >$InsertRange) or ($flg1 ==97 and $insert1 <$InsertRangeMin) or ($flg1 ==161 and $insert1 <$InsertRangeMin)) and ($mateid1 eq "=" and $mate_mapped==1 and $read_mapped==1))
	{
											 		$dup_read_end++;#print "popid",$pop[$i],$lineEnd," \n";
											 		#print OUT2 "R_end","\t",$lineEnd,"\n";
														
											 	}else{}
											 }
											 else
											 { 
											 	## It is reverse strand
											 	if($flag_end > -1)
												{$flag_end=-1;
													$minus_end++;
												}else{$flag_end=-1;}
												if($mateid1 eq "=")
												{
													#print $mateid;
												}else{$R_trans_end++;}
											 	if($read_strand_plus == $mate_strand_plus and $mate_mapped==1 and $read_mapped==1 and $mateid1 eq "=" and $insert1 >=$InsertRange)
											 	{
											 		$invy_end_R++;
											 	}else{}
											 	if($mate_mapped==0 and $read_mapped==1)
												{
													$singleton_R_end++;#print OUT2 "R_end","\t",$line,"\n";
												}else{}
												if((($flg1 ==145 and $insert1 >$InsertRange) or ($flg1 ==81 and $insert1 >$InsertRange) or ($flg1 ==97 and $insert1 <$InsertRangeMin) or ($flg1 ==161 and $insert1 <$InsertRangeMin)) and ($mateid1 eq "=" and $mate_mapped==1 and $read_mapped==1))
	{
													#print "popid",$pop[$i],$lineEnd," \n";
											 		$dup_read_end++;#print OUT2 "F_end","\t",$lineEnd,"\n";
														
											 	}else{}
											 }
									
								}else{} ##end of mapping quality check
                        } ## end of mate region
                        
                        push (@inversions_F_mate,"$invy_end_F");
                        push (@inversions_R_mate,"$invy_end_R");
                        
                        push (@softclip_F_mate,"$plus_end");
                        push (@softclip_R_mate,"$minus_end");
                       # push (@softy,"$total_softy");
                        
                        push (@singleton_F_mate,"$singleton_F_end");
                        push (@singleton_R_mate,"$singleton_R_end");
                        
                        push (@duplication_Reads_mate,"$dup_read_end");
                        $total_AR_mate=$#align_mate;
			            #$total_translocations=$F_trans_end+$R_trans_end;
                        #push (@translocations,"$total_translocations");
                       	push (@totalAlignedReads_mate,"$total_AR_mate");
						@align_mate=();$i++;
				}## end of pop while
			}##end of median_dist check
                        
                        
			
		
		#print OUT "\n";			
		}else{}
		
		if($median_dist >0)
		{
			#print $check_total_reads."\n";
			
			$calc_distance_inv=($median_dist+300)-$st_pos;
		    $inv_reg=$chr."\:".$st_pos."\-".$mate_500down;
		     
		    print OUT $bin."\t".$new_reg_st."\t";
			print OUT $calc_distance_inv."\t".$new_reg_mate."\t".$median_dist."\t".$inv_reg."\t";
			print OUT $median_F,"\t",$median_R,"\t"; 
			
			$calc_distance_inv=0;
			$new_reg_mate=0;$median_dist=0;$inv_reg=0;

			display(@inversions_F);
			display(@inversions_F_mate);
			display(@inversions_R);
			display(@inversions_R_mate);
			
			#display(@softy);
			display(@softclip_F);
			display(@softclip_F_mate);
			display(@softclip_R);
			display(@softclip_R_mate);
			
			#display(@singleton);
			display(@singleton_F);
			display(@singleton_F_mate);
			display(@singleton_R);
			display(@singleton_R_mate);
			
			display(@duplication_Reads);
			display(@duplication_Reads_mate);
		
			#display(@translocations);
			display(@totalAlignedReads);
			display(@totalAlignedReads_mate);
			
			
			print OUT "\n";		
		}else{}## end of printing if
		
		$new_reg_mate=0;$median_dist_dup=0;$median_dist_invy=0;$median_F=0;$median_R=0;
		$st_pos=$end_pos; $end_pos=$end_pos+1000;$new_reg_st=$chr."\:".$st_pos."\-".$end_pos;
	}## end of outer while
	print $outFile." dup job is processed..... NEXT";
}

sub display{
	@dataArray=@_;
	foreach $data (@dataArray)
		{ print OUT $data."#";
		}print OUT "\t";
	@dataArray=();
}
sub sum_reads{
	@dataArray=@_;
	$total=0;
	foreach $data (@dataArray)
		{ $total=$total+$data;
		}
	@dataArray=();
	return $total;
}
sub calc_difference{
	@dataArray=@_;
	@val=();
	$counter=0;
	$total_AD=0;
	$total_FW=0;
	$sumFW=0;
	$sumAD=0;
	foreach $data (@dataArray)
	{ 
		$counter++;
		if($counter <=6)
		{
			$sumAD=$sumAD+$data;
		}
		if($counter >6 and $counter <=9)
		{
			$sumFW=$sumFW+$data;
		}else{}
		
	}@dataArray=();
	
	
	#print OUT "\t";
	push(@val,"$sumAD");
	push(@val,"$sumFW");return @val;
}
### the purpose of this function is to convert the decimal number to binary and return array of binary numbers
sub calculateBitwiseflag{
	#my ($flg) = @_; 
	$flg=$_[0];
	
	@bitwise=();#print $flg;
	$quotient= int($flg / 2);#print $quotient;
	$remainder= $flg % 2;
	
	$flg=$quotient;
	push (@bitwise,"$remainder");
	while($flg > 1)
	{
		$quotient= int($flg / 2);
		$remainder= $flg % 2;
		$flg=$quotient;
		push (@bitwise,"$remainder");
	}
	$remainder=1;
	push (@bitwise,"$remainder");
	return @bitwise;
		
}
sub Deviants{
	
	
	$read_strand=$_[0];
	### here we need to check strand of the query read.. (Do this for + (1) and - (0) strands seperately)
				 ### then in case of inversions.. softclips... single reads we push them to their respective vectors...
				 if($read_strand == 1)
				 { 
				 	## It is forward strand
				 	if($flag_st > -1)
					{$flag_st=-1;
						$plus_st++;
					}else{$flag_st=-1;}
					if($mateid eq "=")
					{
						#print $mateid;
					}else{$F_trans_start++;}
					if($read_strand_plus == $mate_strand_plus)
					{
						$invy_st_F++; 
						$mate_distance=$read_st+30+$insert; # 30 is added here assuming read length of 75 and going in the middle and the add up insert 															
						if(	$mate_distance <= $chr_size-100){						
						#print $mate_distance.." \n";
						push (@MateDistances,"$mate_distance");
						$mate_distance=0;}else{} # reinitializing the position
					}else{}
					if($mate_mapped==0 and $read_mapped==1)
					{
						$singleton_F_st++;
					}else{}
				 }
				 else
				 { 
				 	## It is reverse strand
				 	if($flag_st > -1)
					{$flag_st=-1;
						$minus_st++;
					}else{$flag_st=-1;}
					if($mateid eq "=")
					{
						#print $mateid;
					}else{$R_trans_start++;}
				 	if($read_strand_plus == $mate_strand_plus)
				 	{
				 		$invy_st_R++;
				 		$mate_distance=$read_st+30+$insert; # 30 is added here assuming read length of 75 and going in the middle and the add up insert 															 
						if(	$mate_distance <= $chr_size-100)
						{									###the calculate median (+- 500 bp) to look up for all mates. 
							push (@MateDistances,"$mate_distance");
							$mate_distance=0;}else{} # reinitializing the position
				 	}else{}
				 	if($mate_mapped==0 and $read_mapped==1)
					{
						$singleton_R_st++;
					}else{}
				 }
				
}
sub DeviantsMate{
	$read_strand=$_[0];
	### here we need to check strand of the query read.. (Do this for + (1) and - (0) strands seperately)
				 ### then in case of inversions.. softclips... single reads we push them to their respective vectors...
				
				 if($read_strand_plus == 1)
				 { 
				 	## It is forward strand
				 	if($flag_end  > -1)
					{$flag_end =-1;
						$plus_end++;
					}else{$flag_end=-1;}
					if($mateid1 eq "=")
					{
						#print $mateid;
					}else{$F_trans_end++;}
					if($read_strand_plus == $mate_strand_plus)
					{
						$invy_end_F++; 
					}else{}
					if($mate_mapped==0 and $read_mapped==1)
					{
						$singleton_F_end++;
					}else{}
				 }
				 else
				 { 
				 	## It is reverse strand
				 	if($flag_end > -1)
					{$flag_end=-1;
						$minus_end++;
					}else{$flag_end=-1;}
					if($mateid1 eq "=")
					{
						#print $mateid;
					}else{$R_trans_end++;}
				 	if($read_strand_plus == $mate_strand_plus)
				 	{
				 		$invy_end_R++;
				 	}else{}
				 	if($mate_mapped==0 and $read_mapped==1)
					{
						$singleton_R_end++;
					}else{}
				 }
}
sub translocationFunc{
	if($mateid eq "=")
	{
		#print $mateid;
	}else{$count_trans_start++;}
}
sub initializeBitwiseflag{
	 ### inializing the bam bitwise flags###
		$pairedEndseq=0;
		$pair_proper_mapped=0;
		$read_mapped=0;
		$mate_mapped=0;
		$read_strand_plus=0;
		$mate_strand_plus=0;
		$first_read_inPair=0;
		$second_read_inPair=0;
		$Notprimary_alignment=0;
		$Failvendor_qualitycheck=0;
		$PCR_duplicate=0;
	  #######################################
}
sub extractflagInfo{
	@bitFlags=@_;
	if(@bitFlags[0] == 1)
	{
		$pairedEndseq=1; # paired end
	}
	if(@bitFlags[1] == 1)
	{
		$pair_proper_mapped=1; # pair properly mapped
	}
	if(@bitFlags[2] == 0)
	{
		$read_mapped=1; #read seq is mapped
	}
	if(@bitFlags[3] == 0)
	{
		$mate_mapped=1; #mate is mapped
	}
	if(@bitFlags[4] == 0)
	{
		$read_strand_plus=1; #read strand is plus/forward
	}
	if(@bitFlags[5] == 0)
	{
		$mate_strand_plus=1; #mate strand is plus/forward
	}
	if(@bitFlags[6] == 1)
	{
		$first_read_inPair=1; # first read in pair
	}
	if(@bitFlags[7] == 1)
	{
		$second_read_inPair=1; # second read in pair
	}
	if(@bitFlags[8] == 1)
	{
		$Notprimary_alignment=1; # secondary alignment
	}
	if(@bitFlags[9] == 0)
	{
		$Failvendor_qualitycheck=1; #fail vendor/platform quality check
	}
	if(@bitFlags[10] == 1)
	{
		$PCR_duplicate=1; #PCR duplicate

	}
@bitFlags=();
	
}

close OUT;
close region;



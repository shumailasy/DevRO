use Math::Cephes qw(:explog);
use Statistics::Lite qw(:all);
use Statistics::Zscore;

$file1=shift; ## Region file (chr chr_size)
open (region,"$file1") or die "cannot open file, $!\n";
#@Region=<region>;
$outFile=shift;
$b1 = $outFile."_CountsSingletonReads.Allpop.txt";
open (OUT, ">$b1");
#open (OUT2, ">SingletonR_bam.seq");

@pop=("BAMS_input/1777_sample.dedup.bam","BAMS_input/AngoraMale1_sample.dedup.bam","BAMS_input/CL2_sample.dedup.bam","BAMS_input/CL4_sample.dedup.bam");

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
		my ($chr,$chr_size,$sp_st)= @split_line[0,1,2]; ## full bam file scanned in 1kbp windows
		#my ($chr,$chr_size)= @split_line[0,1];
		#push(@Array_pop_total_DoC_score,$total_DoC_score);
		$st_pos=$sp_st; $end_pos=$st_pos+1000;$new_reg_st=$chr."\:".$st_pos."\-".$end_pos;
		$bin=0;
	while($st_pos <= ($chr_size))
	{
		if($st_pos >0 and $end_pos <$chr_size)
		{$bin++;
			
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
			while($i <= $size_pop)
			{
				#Here we used 20 because of the number of populations (6 doms, 3 FW , 1 Inbredd and 11 WI )
			
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
       			$call_st="samtools view -q 10 $pop[$i] $new_reg_st"; 
       			@align_st=qx($call_st);
			
			
			#processing for start region############################################
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
		               ####################################################
				#@split_cigar=split('',$cigar);
				#foreach $sf (@split_cigar)
				#{
				#	if($sf eq 'S'){$flag_st=1;}else{}
				#}
				if($MQuality >= 10 and $mateid eq "=")
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
							if($read_strand_plus == $mate_strand_plus and $mate_mapped==1 and $read_mapped==1 and $mateid eq "=")
							{
								 
								$mate_distance=$mate_st+50-1; # 50 is added here assuming read length of 50 and going in the end of read			
								if($mate_distance <= $chr_size-1){
									#print $flg." ".$readid."=FLAG_F\t";
									$invy_st_F++;								
								#print $mate_distance.." \n";
								push (@MateDistances,"$mate_distance");
								}else{
								} 
							}else{}
							if($mate_mapped==0 and $read_mapped==1)
							{
								
								$F_read_pos=$read_st+50-1;
								if($F_read_pos <= $chr_size-1)
                                {	#print "popid",$pop[$i],$line,"\n";
                                	#print OUT2 "F_st","\t",$line,"\n";
                                    $singleton_F_st++;
                                    push (@MateDistances_singF,"$F_read_pos");
                                }else{}
								#print OUT2 "F_st","\t",$line,"\n";
							}else{ }
							if((($flg ==145 and $insert >0) or ($flg ==97 and $insert < 0)) and ($MQuality >= 10 and $mateid eq "="))
                            {
                                                                $mate_distance_dup=$mate_st+50-1;
                                                                if($mate_distance_dup <= $chr_size-1)
                                                                {	#print "popid",$pop[$i],$line,"\n";
                                                                	#print OUT2 "F_st","\t",$line,"\n";
                                                                        $dup_read_st++;
                                                                        push (@MateDistances_dup,"$mate_distance_dup");
                                                                }else{}

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
						 		
						 		#print $read_st."=readst".$insert."=insert".$mate_st."=matest........Reverse <<<<<<------------\n";
								if(	$mate_distance <= $chr_size-1)
								{	#print $flg." ".$readid."=FLAG_R\t";
									$invy_st_R++;								###the calculate median (+- 500 bp) to look up for all mates. 
									push (@MateDistances,"$mate_distance");
									}else{#print $invy_st_F."#".$invy_st_R."\t";
									#print "  ______________________ mate distance OUTTTTTTTTTTTTTTTT OF BOUNDDDDDDDDDDDDD\n";
									} # reinitializing the position
						 	}else{}
						 	if($mate_mapped==0 and $read_mapped==1)
							{
								$singleton_R_st++;
								$R_read_pos=$read_st;
								if($R_read_pos <= $chr_size-1)
                                {	#print "popid",$pop[$i],$line,"\n";
                                	#print OUT2 "R_st","\t",$line,"\n";
                                    $singleton_R_st++;
                                    push (@MateDistances_singR,"$R_read_pos");
                                }else{}
							}else{}
							if((($flg ==145 and $insert >0) or ($flg ==97 and $insert < 0)) and ($MQuality >= 10 and $mateid eq "="))                                      
						 	{
						 		$mate_distance_dup=$mate_st;
								if($mate_distance_dup <= $chr_size-1)
                                 {#print "popid",$pop[$i],$line," \n";
									$dup_read_st++;#print OUT2 "R_st","\t",$line,"\n";
									push (@MateDistances_dup,"$mate_distance_dup");
								}else{}
						 		
						 	}else{}
						}
				}else{}
				
			}##end of @align_st
			
                        
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

			$median_dist_dup=int(median(@MateDistances_dup));
			$median_dist_invy=int(median(@MateDistances));
	
		
		#print OUT "\n";			
		}else{}
		
			#print $check_total_reads."\n";
			print OUT $bin."\t".$new_reg_st."\t";
			print OUT $median_F."\t".$median_R."\t".$median_dist_invy."\t".$median_dist_dup."\t";
			
		
			display(@inversions_F);
			#display(@inversions_F_mate);
			display(@inversions_R);
			#display(@inversions_R_mate);
			
			#display(@softy);
			display(@softclip_F);
			#display(@softclip_F_mate);
			display(@softclip_R);
			#display(@softclip_R_mate);
			
			#display(@singleton);
			display(@singleton_F);
			#display(@singleton_F_mate);
			display(@singleton_R);
			#display(@singleton_R_mate);
			
			display(@duplication_Reads);
			#display(@duplication_Reads_mate);
		
			#display(@translocations);
			display(@totalAlignedReads);
			#display(@totalAlignedReads_mate);
			
			print OUT "\n";
			
		$new_reg_mate=0;$median_dist_dup=0;$median_dist_invy=0;$median_F=0;$median_R=0;
		$st_pos=$end_pos; $end_pos=$end_pos+1000;$new_reg_st=$chr."\:".$st_pos."\-".$end_pos;
	}## end of outer while
	print $outFile." Ref-del job is processed... NEXT";
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
	
}

close OUT;
#close OUT2;
close region;


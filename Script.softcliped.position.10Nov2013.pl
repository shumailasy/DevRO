#use Math::Cephes qw(:explog);
#use Statistics::Lite qw(:all);
#use Statistics::Zscore;

$file1=shift or die "Usage: $0 input GenomeFile name: chr size start \n"; ## Region file (chr chr_size)
open (region,"$file1") or die "cannot open file, $!\n";
#@Region=<region>;
$outFile=shift or die "Usage: $1 needs OutputFile prefix \n";
$b1 = $outFile."_softclips.general.txt";
open (OUTfile_soft, ">$b1");
$config=shift or die "Usage: $2 needs configFile, group and PATH to BAMfiles with name \n";
open (conFile,"$config") or die "cannot open file, $!\n";
@pop=();
while($configLine=<conFile>)
{
		chomp $configLine;
		@split_conline=split("\t",$configLine);
		
		my ($group,$BAMPATH)= @split_conline[0,1];
		push (@pop,"$BAMPATH");
#print $BAMPATH;
#@pop=("BAMS_input/1777_sample.dedup.bam","BAMS_input/AngoraMale1_sample.dedup.bam");
}

print "bin"."\t"."new_reg_st"."\t"."softclips Tab singletons Tab overlaps Tab total","\n";
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
	
		$count_softclip_start=0;
		$count_softclip_end=0;
		$flag_st=-1;
		$flag_end=-1;
		chomp $reg;#print OUT $reg."\t";
		@split_line=split("\t",$reg);
		my ($chr,$start1,$end1,$chr_size,$info)= @split_line[0,1,2,3,4];print $info;
		$filename1 = "Allpop"."Softclip_".$info.".bed";
    	open (OUTfile, ">$filename1");
		$st_pos=$start1; $end_pos=$end1;$new_reg_st=$chr."\:".$st_pos."\-".$end_pos;
		$bin=0;
	#while($st_pos <= ($chr_size))
	#{
		if($st_pos >0 and $end_pos <$chr_size)
		{$bin++;print $bin."\t".$info."\t".$new_reg_st."\t";
			
			$i=0;
			@inversions=();
			@inversions_F=();
			@inversions_R=();
			@softclip_F=();
			@softclip_R=();
			@inversions_F_mate=();
			@inversions_R_mate=();
			@softclip_F_mate=();
			@softclip_R_mate=();
			@softy=();
			@singleton=();
			@singleton_F=();
			@singleton_R=();
			@Overlap_reads=();
			@singleton_F_mate=();
			@singleton_R_mate=();
			@duplication_Reads=();
			@duplication_Reads_mate=();
			@translocations=();
			@totalAlignedReads=();
			@totalAlignedReads_mate=();
			@MateDistances=();@MateDistances_dup=();
			@MateDistances_singF=();@MateDistances_singR=();
			@narrow_F=();@narrow_R=();
			@narrow_dist=();
		    
		    
			@Array_hashes=();
			### extracting information for population for that region::
			$size_pop=$#pop;
			while($i <= $size_pop)
			{
				#Here we used 20 because of the number of populations (6 doms, 3 FW , 1 Inbredd and 11 WI )
					%hash_softclips = ();
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
                	$Overlaps_count=0;
                	$singleton_R_end=0;
                	$dup_read_st=0;
                	$dup_read_end=0;
                	$count_softclip_plus=0;$count_softclip_minus=0;
                	for($index=$st_pos;$index<=$end_pos;$index++)
					{
						## initialize hash with the size of region
						$val=0;
						$hash_softclips{$index}=$val;
					}					
					#@MateDistances_R=();
					####
                	######	
       			$call_st="samtools view $pop[$i] $new_reg_st"; 
       			@align_st=qx($call_st);
			#$call_end="samtools view $pop[$i] $new_reg_end";@align_end=qx($call_end);print $call_end."\n\n";
			#print $i."\n";print $align_st[0]."\n";
			%readid_array=();
			#processing each line of bam file in each population
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
					                	$softclip_position=$softclip_position+$cigar_positions[$sf];
					                }
					                else{}
										
						}
						
		               
		               #####################################################
		         $dummy_val=1;      
				if($MQuality >= 0 and $mate_mapped==1 and $read_mapped==1 and $mateid eq "=" )
				{	
					#Deviants($read_strand_plus);
						 if($read_strand_plus == 1)
						 { 
						 	## It is forward strand
						 	
								if($flag_st==1)
								{ #5'end
									
									if(exists $hash_softclips{$prime5_softclip_position}) {
								        $hash_softclips{$prime5_softclip_position}=$hash_softclips{$prime5_softclip_position}+1;
								        $readid_array{$prime5_softclip_position}=$cigar;
								       # print $i."	".$flag_st." ".$prime5_softclip_position." ".$hash_softclips{$prime5_softclip_position}."\n";												
								        $flag_st=-1;
								    }else 
								    {
								        $hash_softclips{$prime5_softclip_position}=$dummy_val;$flag_st=-1;
								        $readid_array{$prime5_softclip_position}=$cigar;
								    }
									
								}
								elsif($flag_st==2)
								{ #3'end
									$prime3_softclip_position=$prime3_softclip_position-1;
									if (exists $hash_softclips{$prime3_softclip_position}) {
								        $hash_softclips{$prime3_softclip_position}=$hash_softclips{$prime3_softclip_position}+1;$flag_st=-1;
								        $readid_array{$prime3_softclip_position}=$cigar;
								        
								    } else {
								        $hash_softclips{$prime3_softclip_position} = $dummy_val;$flag_st=-1;
								        $readid_array{$prime3_softclip_position}=$cigar;
								    }
								}else{}
								$flag_st=-1;
								$plus_st++;
							
							
						 }
						 else
						 { 
						 	## It is reverse strand
						 		if($flag_st==1)
								{ #5'end
									$flag_st=-1;
									if(exists $hash_softclips{$prime5_softclip_position}) {
								        $hash_softclips{$prime5_softclip_position}=$hash_softclips{$prime5_softclip_position}+1;
								        $readid_array{$prime5_softclip_position}=$cigar;
								    }else{
								        $hash_softclips{$prime5_softclip_position} = 1;$readid_array{$prime5_softclip_position}=$cigar;
								    }
									
								}
								elsif($flag_st==2)
								{ #3'end
									$flag_st=-1;$prime3_softclip_position=$prime3_softclip_position-1;
									if(exists $hash_softclips{$prime3_softclip_position}) {
								        $hash_softclips{$prime3_softclip_position}=$hash_softclips{$prime3_softclip_position}+1;
								        $readid_array{$prime3_softclip_position}=$cigar;
								    }else{
								        $hash_softclips{$prime3_softclip_position} = 1;
								        $readid_array{$prime3_softclip_position}=$cigar;
								    }
								}else{}
								$flag_st=-1;
								$plus_st++;
							
							
						}
				}else{}
			
			} #print $call_st."\n";print $invy_st_F."#".$invy_st_R."\n";
						$pop_i=$i+1;
						
						while (my ($k,$v)=each %hash_softclips){$dumy_st=$k-1;print OUTfile $chr."\t".$k."\t".$k."\t".$v."\t".$pop[$i]."\t".$info."\n";
						if($v >0){print $chr."\t".$dumy_st."\t".$k."\t".$v."\t"."pop".$pop_i."\t";print $readid_array{$k}."\n";
						}}
						push (@Array_hashes,\%hash_softclips);
						%hash_softclips = ();
						
                        $total_AR=$#align_st;	
                        push (@softclip_F,"$plus_st");
                       	push (@totalAlignedReads,"$total_AR");
			@align_st=();
			$i++;
			}## end of population while loop

		}else{}
	

		for($index1=$st_pos;$index1<=$end_pos;$index1++)
		{
			print OUTfile_soft $chr."\t".$index1."\t";
			 foreach (@Array_hashes){print OUTfile_soft $_->{$index1}."\t";print OUTfile_soft "\n";}
			#print OUTfile_soft "\n";
		}
		for (keys %hash_softclips)
		{
			delete $hash_softclips{$_};
		}
		
		print OUT $new_reg_st."\t";
		#display(@softy);
		display(@softclip_F);
		display(@totalAlignedReads);
		print OUT "\n";		
		close OUTfile;
		#$st_pos=$end_pos; $end_pos=$end_pos+1000;$new_reg_st=$chr."\:".$st_pos."\-".$end_pos;
	#}## end of outer while
}## end of region file loop

sub display{
	@dataArray=@_;
	foreach $data (@dataArray)
		{ print OUT $data."#";
		}print OUT "\t";
	@dataArray=();
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
close region;


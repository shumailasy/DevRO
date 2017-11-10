
#			print OUT $bin."\t".$new_reg_st."\t";
#			print OUT $median_F."\t".$median_R."\t".$median_dist_invy."\t".$median_dist_dup."\t";
#			display(@inversions_F);
#			display(@inversions_R);
#			display(@softclip_F);
#			display(@softclip_R);
#			display(@singleton_F);
#			display(@singleton_R);
#			display(@duplication_Reads);
#			display(@totalAlignedReads);

#########

$file1=shift; ## Region file (chr chr_size)
open (region,"$file1") or die "cannot open file, $!\n";

$outFile=shift;
$b1 = $outFile."_"."REFDEL_OutParser_WholeSet.bed";
open (DUPPARSER, ">$b1");
#$c1 = $outFile."_"."DELDOM_OutParser_WholeSet.bed";
#open (DUPPARSER, ">$c1");

#$c1 = $outFile."_"."DEL_OutParser.bed";
#open (DELPARSER, ">$c1");
#$INSP = $outFile."_"."INSERT_OutParser.bed";
#open (INSPARSER, ">$INSP");

while($reg=<region>){
	chomp $reg;#print OUT $reg."\t";
		@split_line=split("\t",$reg);
		
		my ($bin,$new_reg_st,$median_F,$median_R,$median_dist_invy,$median_dist_dup,$A_invy1,$A_invy2,$A_soft1,$A_soft2,$A_single1,$A_single2,$A_dup1, $A_Total1)= @split_line[0,1,2,3,4,5,6,7,8,9,10,11,12,13];
		
		
		@invy1=split("#",$A_invy1);@invy2=split("#",$A_invy2);
		@soft1=split("#",$A_soft1);@soft2=split("#",$A_soft2);
		@single1=split("#",$A_single1);@single2=split("#",$A_single2);
		@dup1=split("#",$A_dup1);
		@Total1=split("#",$A_Total1);
		#@DEL1=split("#",$A_DELETION);
		#@INS1=split("#",$A_INSERTION);
		
		### Need to edit here##


#################
		
		$AD_total=$Total1[0]+$Total1[1];
		$AW_total=$Total1[2]+$Total1[3];
		
		## Calculating the M-values on the breakpoints total reads:
		if($AD_total > 0 and $AW_total> 0){$Mval=$AW_total/$AD_total;}else{$Mval=-1;}

		$sum_total=$AD_total+$AW_total;
		
		$AD_invy=$invy1[0]+$invy1[1]+$invy2[0]+$invy2[1];
		$AW_invy=$invy1[2]+$invy1[3]+$invy2[2]+$invy2[3];

		$AD_soft=$soft1[0]+$soft1[1]+$soft2[0]+$soft2[1];
		$AW_soft=$soft1[2]+$soft1[3]+$soft2[2]+$soft2[3];
		
		$AD_single=$single1[0]+$single1[1]+$single2[0]+$single2[1];
		$AW_single=$single1[2]+$single1[3]+$single2[2]+$single2[3];
	
		$AD_dup=$dup1[0]+$dup1[1];
		$AW_dup=$dup1[2]+$dup1[3];
		
@position=split("\:",$new_reg_st);
@chromo=split("chr",$position[0]);
$chrom1=$chromo[1];
@coordinates=split("\-",$position[1]);
## check this later to see if somehing is  not miseed
$st_coord=$coordinates[0];$end_coord=$coordinates[1];if($st_coord < $end_coord){}else{$st_coord=$coordinates[1];$end_coord=$coordinates[0];}

#if($median_del > 0){$size_DEL=$median_del-$st_coord;}else{$size_DEL=-1;}
#if($median_ins > 0){$size_INS=$median_ins-$st_coord;}else{$size_INS=-1;}
######

if(($AW_soft+$AD_soft)>0){$frac_soft=$AW_soft/($AW_soft+$AD_soft);}else{$frac_soft=-1;}
if(($AD_single+$AW_single)>0){$frac_single=$AW_single/($AD_single+$AW_single);}else{$frac_single=-1;}
if(($AW_invy+$AD_invy)>0){$frac_invy=$AW_invy/($AW_invy+$AD_invy);}else{$frac_invy=-1;}
if(($AW_dup+$AD_dup)>0){$frac_dup=$AW_dup/($AW_dup+$AD_dup);}else{$frac_dup=-1;}
print $bin," ",$inv_reg," ",$frac_dup,"\n";
#if(($AW_DEL1+$AD_DEL1)>0){$frac_DEL1=$AW_DEL1/($AW_DEL1+$AD_DEL1);}else{$frac_DEL1=-1;}
#if(($AW_INS+$AD_INS)>0){$frac_INS=$AW_INS/($AW_INS+$AD_INS);}else{$frac_INS=-1;}	

$reg=~ tr/"\t"/"\:"/;

## abs. fraction of deviant reads
 $AD_sum_dev=$AD_soft+$AD_single+$AD_invy+$AD_dup;
 $AW_sum_dev=$AW_soft+$AW_single+$AW_invy+$AW_dup;
 if(($AD_total)>0){$AD_devFrac=$AD_sum_dev/($AD_total);}else{$AD_devFrac=-1;}
 if(($AW_total)>0){$AW_devFrac=$AW_sum_dev/($AW_total);}else{$AW_devFrac=-1;}
 $abs_diff_dev=abs($AW_devFrac-$AD_devFrac);

##

## need to edit here $frac_invy ?? 
##################################### To DO List 
## Check why we donot get any duplication in doms? relaxing criteria and Mvalues based on total reads..

## this is  to get almost everything pssible greater than size 0 , 1) size issue 2) invy issue 3)fraction_dup 4)frac_abs_dev
	if($frac_single >=0.9)
	{
		print DUPPARSER $chrom1,"\t",$st_coord,"\t",$end_coord,"\t",$new_reg_st,"\t",$median_F,"\t",$median_R,"\t",$median_dist_invy,"\t",$median_dist_dup,"\t",$frac_invy,"\t",$frac_soft,"\t",$frac_single,"\t",$frac_dup,"\t",$AD_invy,"\t",$AW_invy,"\t",$AD_dup,"\t",$AW_dup,"\t",$AD_single,"\t",$AW_single,"\t",$AD_soft,"\t",$AW_soft,"\t",$AD_total,"\|",$AW_total,"\t",$reg,"\t","\t",$AD_devFrac,"\t",$AW_devFrac,"\t",$abs_diff_dev,"\t","pileupsWilds.RefDEL.WholeSet","\t",$Mval,"\n";
	 }else{}
	 if($frac_single <=0.1)
	{
		#print DUPPARSER $chrom1,"\t",$st_coord,"\t",$end_coord,"\t",$new_reg_st,"\t",$median_F,"\t",$median_R,"\t",$median_dist_invy,"\t",$median_dist_dup,"\t",$frac_invy,"\t",$frac_soft,"\t",$frac_single,"\t",$frac_dup,"\t",$AD_invy,"\t",$AW_invy,"\t",$AD_dup,"\t",$AW_dup,"\t",$AD_single,"\t",$AW_single,"\t",$AD_soft,"\t",$AW_soft,"\t",$AD_total,"\|",$AW_total,"\t",$reg,"\t","\t",$AD_devFrac,"\t",$AW_devFrac,"\t",$abs_diff_dev,"\t","pileupsDOM.WholeSet","\t",$Mval,"\n";
	}else{}
}

sub Min{
	$val1=$_[0];
	$val2=$_[1];
	$min=0;
	if($val1 > $val2)
	{
		$min=$val2;
	}
	else
	{
		$min=$val1;
	}#print $min;
	return $min;
}

## Running the script
##per


$file1=shift; ## Region file (chr chr_size)
open (region,"$file1") or die "cannot open file, $!\n";

$outFile=shift;
$b1 = $outFile."_"."INV_OutParser_chr16.bed";
open (INVPARSER, ">$b1");

#$c1 = $outFile."_"."DEL_OutParser.bed";
#open (DELPARSER, ">$c1");
#$INSP = $outFile."_"."INSERT_OutParser.bed";
#open (INSPARSER, ">$INSP");

$inv=0;$del=0;$isn=0;
while($reg=<region>){
	chomp $reg;#print OUT $reg."\t";
		@split_line=split("\t",$reg);
		
		my ($bin,$new_reg_st,$calc_distance_inv,$new_reg_mate,$median_dist,$inv_reg,$median_F,$median_R,$A_invy1,$A_invy2,$A_invy3,$A_invy4,$A_soft1,$A_soft2,$A_soft3,$A_soft4,$A_single1,$A_single2,$A_single3,$A_single4,$A_dup1,$A_dup2,$A_Total1,$A_Total2,$A_DELETION,$A_INSERTION,$median_del,$median_ins,$median_dist_dup)= @split_line[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28];
		
		
		@invy1=split("#",$A_invy1);@invy2=split("#",$A_invy2);@invy3=split("#",$A_invy3);@invy4=split("#",$A_invy4);
		@soft1=split("#",$A_soft1);@soft2=split("#",$A_soft2);@soft3=split("#",$A_soft3);@soft4=split("#",$A_soft4);
		@single1=split("#",$A_single1);@single2=split("#",$A_single2);@single3=split("#",$A_single3);@single4=split("#",$A_single4);
		@dup1=split("#",$A_dup1);@dup2=split("#",$A_dup2);
		@Total1=split("#",$A_Total1);@Total2=split("#",$A_Total2);
		@DEL1=split("#",$A_DELETION);
		@INS1=split("#",$A_INSERTION);
		
		### Need to edit here##

## for Inversions checking the minimum invy reads in FF or RR in reg and mate_reg 
$wf1=Min($invy1[0],$invy2[0]);$Wf1=$wf1+$wf1; ## for forward FF Dom1
$wr1=Min($invy3[0]+$invy4[0]);$Wr1=$wr1+$wr1; ## for Reverse RR Dom1

$wf2=Min($invy1[1],$invy2[1]);$Wf2=$wf2+$wf2; ## for forward FF Dom2
$wr2=Min($invy3[1]+$invy4[1]);$Wr2=$wr2+$wr2; ## for Reverse RR Dom2

$wf3=Min($invy1[2],$invy2[2]);$Wf3=$wf3+$wf3; ## for forward FF Wild1
$wr3=Min($invy3[2]+$invy4[2]);$Wr3=$wr3+$wr3; ## for Reverse RR Wild1

$wf4=Min($invy1[3],$invy2[3]);$Wf4=$wf4+$wf4; ## for forward FF Wild2
$wr4=Min($invy3[3]+$invy4[3]);$Wr4=$wr4+$wr4; ## for Reverse RR Wild2

$flg_dom_inv=0;
$flg_wild_inv=0;

if(($Wf2+$Wr2) == 0 ){$flg_dom_inv=0;}else{$flg_dom_inv=1;}
if(($Wf4+$Wr4) == 0){$flg_wild_inv=0;}else{$flg_wild_inv=2;}


## for duplications checking the minimum dups reads in RF reg and mate_reg
$wd1=Min($dup1[0],$dup2[0]);$Wd1=$wd1+$wd1; ##dom1
$wd2=Min($dup1[1],$dup2[1]);$Wd2=$wd2+$wd2; ##dom2

$wd3=Min($dup1[2],$dup2[2]);$Wd3=$wd3+$wd3; ## wild1
$wd4=Min($dup1[3],$dup2[3]);$Wd4=$wd4+$wd4; ## wild2

$flg_dom=0;
$flg_wild=0;

if($Wd1 == 0 or $Wd2==0){$flg_dom=0;}else{$flg_dom=1;}
if($Wd3 == 0 or $Wd4==0){$flg_wild=0;}else{$flg_wild=2;}

#################
		
		$AD_total=$Total1[0]+$Total1[1]+$Total2[0]+$Total2[1];
		$AW_total=$Total1[2]+$Total1[3]+$Total2[2]+$Total2[3];
		
		## Calculating the M-values on the breakpoints total reads:
		if($AD_total > 0 and $AW_total> 0){$Mval=$AW_total/$AD_total;}else{$Mval=-1;}

		$sum_total=$AD_total+$AW_total;

		$AD_invy=$Wf2+$Wr2; ## dom1 and dom2 forward and reverse reads
		$AW_invy=$Wf4+$Wr4; ## wild1 and wild2 forward and reverse reads
		
		$raw_AD_invy=$invy1[1]+$invy2[1]+$invy3[1]+$invy4[1];
		$raw_AW_invy=$invy1[3]+$invy2[3]+$invy3[3]+$invy4[3];

		$AD_soft=$soft1[0]+$soft1[1]+$soft2[0]+$soft2[1]+$soft3[0]+$soft3[1]+$soft4[0]+$soft4[1];
		$AW_soft=$soft1[2]+$soft1[3]+$soft2[2]+$soft2[3]+$soft3[2]+$soft3[3]+$soft4[2]+$soft4[3];
		
		$AD_single=$single1[0]+$single1[1]+$single2[0]+$single2[1]+$single3[0]+$single3[1]+$single4[0]+$single4[1];
		$AW_single=$single1[2]+$single1[3]+$single2[2]+$single2[3]+$single3[2]+$single3[3]+$single4[2]+$single4[3];
		
		$AD_dup=$Wd1+$Wd2; ## dom1 and dom2 dup reads
		$AW_dup=$Wd3+$Wd4; ## wild1 and wild2 dup reads
	
		$raw_AD_dup=$dup1[0]+$dup1[1]+$dup2[0]+$dup2[1];
		$raw_AW_dup=$dup1[2]+$dup1[3]+$dup2[2]+$dup2[3];
		
		$AD_DEL1=$DEL1[0]+$DEL1[1];
		$AW_DEL1=$DEL1[2]+$DEL1[3];
		
		$AD_INS=$INS1[0]+$INS1[1];
		$AW_INS=$INS1[2]+$INS1[3];
		
@position=split("\:",$inv_reg);
@chromo=split("chr",$position[0]);
$chrom1=$chromo[1];
@coordinates=split("\-",$position[1]);
## check this later to see if somehing is  not miseed
$st_coord=$coordinates[0];$end_coord=$coordinates[1];if($st_coord < $end_coord){}else{$st_coord=$coordinates[1];$end_coord=$coordinates[0];}

if($median_del > 0){$size_DEL=$median_del-$st_coord;}else{$size_DEL=-1;}
if($median_ins > 0){$size_INS=$median_ins-$st_coord;}else{$size_INS=-1;}
######

if(($AW_soft+$AD_soft)>0){$frac_soft=$AW_soft/($AW_soft+$AD_soft);}else{$frac_soft=-1;}
if(($AD_single+$AW_single)>0){$frac_single=$AW_single/($AD_single+$AW_single);}else{$frac_single=-1;}

if(($AW_invy+$AD_invy)>0){$frac_invy=$AW_invy/($AW_invy+$AD_invy);}else{$frac_invy=-1;}
if(($raw_AW_invy+$raw_AD_invy)>0){$raw_frac_invy=$raw_AW_invy/($raw_AW_invy+$raw_AD_invy);}else{$raw_frac_invy=-1;}

if(($AW_dup+$AD_dup)>0){$frac_dup=$AW_dup/($AW_dup+$AD_dup);}else{$frac_dup=-1;}
if(($raw_AW_dup+$raw_AD_dup)>0){$raw_frac_dup=$raw_AW_dup/($raw_AW_dup+$raw_AD_dup);}else{$raw_frac_dup=-1;}
print $bin," ",$inv_reg," INVERSIONNNNNNNNN  ",$calc_distance_inv ."	".$raw_frac_invy."#".$flg_wild_inv."_".$flg_dom_inv."\n";

if(($AW_DEL1+$AD_DEL1)>0){$frac_DEL1=$AW_DEL1/($AW_DEL1+$AD_DEL1);}else{$frac_DEL1=-1;}
if(($AW_INS+$AD_INS)>0){$frac_INS=$AW_INS/($AW_INS+$AD_INS);}else{$frac_INS=-1;}	

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
	
	 if($calc_distance_inv > 0)
	{
		print "HELLLLLLLLLLO ",$flg_wild_inv."_".$flg_dom_inv."\n";
		print INVPARSER $chrom1,"\t",$st_coord,"\t",$end_coord,"\t",$new_reg_st,"\t",$new_reg_mate,"\t",$calc_distance_inv,"\t",$frac_invy,"\t",$frac_soft,"\t",$frac_single,"\t",$frac_dup,"\t",$frac_DEL1,"\t",$frac_INS,"\t",$AD_invy,"\t",$AW_invy,"\t",$AD_dup,"\t",$AW_dup,"\t",$AD_DEL1,"\t",$AW_DEL1,"\t",$AD_INS,"\t",$AW_INS,"\t",$AD_single,"\t",$AW_single,"\t",$AD_soft,"\t",$AW_soft,"\t",$AD_total,"\|",$AW_total,"\t",$reg,"\t",$median_F,"\t",$median_R,"\t",$AD_devFrac,"\t",$AW_devFrac,"\t",$abs_diff_dev,"\t","INVERSION.WholeSet","\t",$Mval,"\t",$raw_frac_invy,"\t",$raw_frac_dup,"	LocusID",$inv,"\n";
		$inv++;
		
	}else{}

}
print "##".$inv."	processed  INVERSIONS \n";
#print $del."	processed  DELETIONS \n";
#print $isn."	processed  INSERTIONS \n";

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

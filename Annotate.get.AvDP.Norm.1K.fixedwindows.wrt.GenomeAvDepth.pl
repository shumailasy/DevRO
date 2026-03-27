
#use Math::Cephes qw(:explog);
#use Statistics::Lite qw(:all);
#use Statistics::Zscore;

$file1=shift; ## DOC file (chr chr_size)
open (region,"$file1") or die "cannot open file, $!\n";
#@Region=<region>;
$outFile=shift;
$b1 = $outFile."_"."OUT.txt";
$nor= $outFile."_"."OUT_normalized_wrt_Av.genome.txt";
open (OUT, ">$b1");
#open (OUT2, ">$nor");
@Av_DP=();
$counter=-1;

### Score for normalization w.r.t. FlemishGaint (av.depth in unknown chr)
$Domestic_FlemishGiant_1=10.82;
$Domestic_FrenchLop_1=10.72;
$WildFrench_Villemolaque_1=10.44;
$Domestic_NetherlandDwarf_1=10.37;
$WildIberian_Calzada_1=10.05;
$WildIberian_SCMora_1=10.04;
$WildFrench_LaRoque_1=10.01;
$WildIberian_Carrion_1=9.95;
$Domestic_NewZealand_1=9.88;
$WildIberian_CO_1=9.75;
$Inbred__1=9.66;
$Domestic_Champagne_1=9.6;
$WildIberian_Guadalajara_1=9.44;
$WildIberian_Castanar_1=9.4;
$WildIberian_TO_1=9.38;
$Domestic_BelgianHare_1=9.31;
$WildFrench_Caumont=9.31;
$WildIberian_M_1=9.24;
$WildIberian_Huelva=9.21;
$Domestic_Dutch_1=9.06;
$Wild_PortoSanto_1=9.02;
$WildIberian_Toledo_1=8.83;
$WildIberian_Zaragoza_1=5.3;
$Lepusamericanus_outgroup_1=4.88;

######

$Depth_for_Domestic_BelgianHare_1_sum=0;$Depth_for_Domestic_Champagne_1_sum=0;$Depth_for_Domestic_Dutch_1_sum=0;$Depth_for_Domestic_FlemishGiant_1_sum=0;$Depth_for_Domestic_FrenchLop_1_sum=0;$Depth_for_Domestic_NetherlandDwarf_1_sum=0;$Depth_for_Domestic_NewZealand_1_sum=0;$Depth_for_Inbred__1_sum=0;$Depth_for_Lepusamericanus_outgroup_1_sum=0;$Depth_for_WildFrench_Caumont_sum=0;$Depth_for_WildFrench_LaRoque_1_sum=0;$Depth_for_WildFrench_Villemolaque_1_sum=0;$Depth_for_WildIberian_CO_1_sum=0;$Depth_for_WildIberian_Calzada_1_sum=0;$Depth_for_WildIberian_Carrion_1_sum=0;$Depth_for_WildIberian_Castanar_1_sum=0;$Depth_for_WildIberian_Guadalajara_1_sum=0;$Depth_for_WildIberian_Huelva_sum=0;$Depth_for_WildIberian_M_1_sum=0;$Depth_for_WildIberian_SCMora_1_sum=0;$Depth_for_WildIberian_TO_1_sum=0;$Depth_for_WildIberian_Toledo_1_sum=0;$Depth_for_WildIberian_Zaragoza_1_sum=0;$Depth_for_Wild_PortoSanto_1_sum=0;
$temp=0;
while($line=<region>)
{
	print "$line";
	chomp $line;
	@split_line=split("\t",$line);
	my ($Locus,$Total_Depth,$Average_Depth_sample,$Depth_for_Domestic_BelgianHare_1,$Depth_for_Domestic_Champagne_1,$Depth_for_Domestic_Dutch_1,$Depth_for_Domestic_FlemishGiant_1,$Depth_for_Domestic_FrenchLop_1,$Depth_for_Domestic_NetherlandDwarf_1,$Depth_for_Domestic_NewZealand_1,$Depth_for_Inbred__1,$Depth_for_Lepusamericanus_outgroup_1,$Depth_for_WildFrench_Caumont,$Depth_for_WildFrench_LaRoque_1,$Depth_for_WildFrench_Villemolaque_1,$Depth_for_WildIberian_CO_1,$Depth_for_WildIberian_Calzada_1,$Depth_for_WildIberian_Carrion_1,$Depth_for_WildIberian_Castanar_1,$Depth_for_WildIberian_Guadalajara_1,$Depth_for_WildIberian_Huelva,$Depth_for_WildIberian_M_1,$Depth_for_WildIberian_SCMora_1,$Depth_for_WildIberian_TO_1,$Depth_for_WildIberian_Toledo_1,$Depth_for_WildIberian_Zaragoza_1,$Depth_for_Wild_PortoSanto_1)= @split_line[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26];
	$counter++;
	
	@split_locus=split("\:",$Locus);
	my ($chromo,$coordinate)= @split_locus[0,1];
	
	
$Depth_for_Domestic_BelgianHare_1_sum+=$Depth_for_Domestic_BelgianHare_1;
$Depth_for_Domestic_Champagne_1_sum+=$Depth_for_Domestic_Champagne_1;
$Depth_for_Domestic_Dutch_1_sum+=$Depth_for_Domestic_Dutch_1;
$Depth_for_Domestic_FlemishGiant_1_sum+=$Depth_for_Domestic_FlemishGiant_1;
$Depth_for_Domestic_FrenchLop_1_sum+=$Depth_for_Domestic_FrenchLop_1;
$Depth_for_Domestic_NetherlandDwarf_1_sum+=$Depth_for_Domestic_NetherlandDwarf_1;
$Depth_for_Domestic_NewZealand_1_sum+=$Depth_for_Domestic_NewZealand_1;
$Depth_for_Inbred__1_sum+=$Depth_for_Inbred__1;
$Depth_for_Lepusamericanus_outgroup_1_sum+=$Depth_for_Lepusamericanus_outgroup_1;
$Depth_for_WildFrench_Caumont_sum+=$Depth_for_WildFrench_Caumont;
$Depth_for_WildFrench_LaRoque_1_sum+=$Depth_for_WildFrench_LaRoque_1;
$Depth_for_WildFrench_Villemolaque_1_sum+=$Depth_for_WildFrench_Villemolaque_1;
$Depth_for_WildIberian_CO_1_sum+=$Depth_for_WildIberian_CO_1;
$Depth_for_WildIberian_Calzada_1_sum+=$Depth_for_WildIberian_Calzada_1;
$Depth_for_WildIberian_Carrion_1_sum+=$Depth_for_WildIberian_Carrion_1;
$Depth_for_WildIberian_Castanar_1_sum+=$Depth_for_WildIberian_Castanar_1;
$Depth_for_WildIberian_Guadalajara_1_sum+=$Depth_for_WildIberian_Guadalajara_1;
$Depth_for_WildIberian_Huelva_sum+=$Depth_for_WildIberian_Huelva;
$Depth_for_WildIberian_M_1_sum+=$Depth_for_WildIberian_M_1;
$Depth_for_WildIberian_SCMora_1_sum+=$Depth_for_WildIberian_SCMora_1;
$Depth_for_WildIberian_TO_1_sum+=$Depth_for_WildIberian_TO_1;
$Depth_for_WildIberian_Toledo_1_sum+=$Depth_for_WildIberian_Toledo_1;
$Depth_for_WildIberian_Zaragoza_1_sum+=$Depth_for_WildIberian_Zaragoza_1;
$Depth_for_Wild_PortoSanto_1_sum+=$Depth_for_Wild_PortoSanto_1;
	
	
	#if($counter==1000)
	if(($coordinate%1000) == 0)
	{
		#take average
$av1=$Depth_for_Domestic_BelgianHare_1_sum/1000;
$av2=$Depth_for_Domestic_Champagne_1_sum/1000;
$av3=$Depth_for_Domestic_Dutch_1_sum/1000;
$av4=$Depth_for_Domestic_FlemishGiant_1_sum/1000;
$av5=$Depth_for_Domestic_FrenchLop_1_sum/1000;
$av6=$Depth_for_Domestic_NetherlandDwarf_1_sum/1000;
$av7=$Depth_for_Domestic_NewZealand_1_sum/1000;
$av8=$Depth_for_Inbred__1_sum/1000;
$av9=$Depth_for_Lepusamericanus_outgroup_1_sum/1000;
$av10=$Depth_for_WildFrench_Caumont_sum/1000;
$av11=$Depth_for_WildFrench_LaRoque_1_sum/1000;
$av12=$Depth_for_WildFrench_Villemolaque_1_sum/1000;
$av13=$Depth_for_WildIberian_CO_1_sum/1000;
$av14=$Depth_for_WildIberian_Calzada_1_sum/1000;
$av15=$Depth_for_WildIberian_Carrion_1_sum/1000;
$av16=$Depth_for_WildIberian_Castanar_1_sum/1000;
$av17=$Depth_for_WildIberian_Guadalajara_1_sum/1000;
$av18=$Depth_for_WildIberian_Huelva_sum/1000;
$av19=$Depth_for_WildIberian_M_1_sum/1000;
$av20=$Depth_for_WildIberian_SCMora_1_sum/1000;
$av21=$Depth_for_WildIberian_TO_1_sum/1000;
$av22=$Depth_for_WildIberian_Toledo_1_sum/1000;
$av23=$Depth_for_WildIberian_Zaragoza_1_sum/1000;
$av24=$Depth_for_Wild_PortoSanto_1_sum/1000;
		
		print OUT $chromo,"\t",$temp,"\t",$coordinate,"\t",$Locus,":",$counter,"\t",$av1,"\t",$av2,"\t",$av3,"\t",$av4,"\t",$av5,"\t",$av6,"\t",$av7,"\t",$av8,"\t",$av9,"\t",$av10,"\t",$av11,"\t",$av12,"\t",$av13,"\t",$av14,"\t",$av15,"\t",$av16,"\t",$av17,"\t",$av18,"\t",$av19,"\t",$av20,"\t",$av21,"\t",$av22,"\t",$av23,"\t",$av24,"\n";
		
		
		### normalization
	
$normlz1=$av1/$Domestic_BelgianHare_1 if $Domestic_BelgianHare_1;
$normlz2=$av2/$Domestic_Champagne_1 if $Domestic_Champagne_1;
$normlz3=$av3/$Domestic_Dutch_1 if $Domestic_Dutch_1;
$normlz4=$av4/$Domestic_FlemishGiant_1 if $Domestic_FlemishGiant_1;
$normlz5=$av5/$Domestic_FrenchLop_1 if $Domestic_FrenchLop_1;
$normlz6=$av6/$Domestic_NetherlandDwarf_1 if $Domestic_NetherlandDwarf_1;
$normlz7=$av7/$Domestic_NewZealand_1 if $Domestic_NewZealand_1;
$normlz8=$av8/$Inbred__1 if $Inbred__1;
$normlz9=$av9/$Lepusamericanus_outgroup_1 if $Lepusamericanus_outgroup_1;
$normlz10=$av10/$WildFrench_Caumont if $WildFrench_Caumont;
$normlz11=$av11/$WildFrench_LaRoque_1 if $WildFrench_LaRoque_1;
$normlz12=$av12/$WildFrench_Villemolaque_1 if $WildFrench_Villemolaque_1;
$normlz13=$av13/$WildIberian_CO_1 if $WildIberian_CO_1;
$normlz14=$av14/$WildIberian_Calzada_1 if $WildIberian_Calzada_1;
$normlz15=$av15/$WildIberian_Carrion_1 if $WildIberian_Carrion_1;
$normlz16=$av16/$WildIberian_Castanar_1 if $WildIberian_Castanar_1;
$normlz17=$av17/$WildIberian_Guadalajara_1 if $WildIberian_Guadalajara_1;
$normlz18=$av18/$WildIberian_Huelva if $WildIberian_Huelva;
$normlz19=$av19/$WildIberian_M_1 if $WildIberian_M_1;
$normlz20=$av20/$WildIberian_SCMora_1 if $WildIberian_SCMora_1;
$normlz21=$av21/$WildIberian_TO_1 if $WildIberian_TO_1;
$normlz22=$av22/$WildIberian_Toledo_1 if $WildIberian_Toledo_1;
$normlz23=$av23/$WildIberian_Zaragoza_1 if $WildIberian_Zaragoza_1;
$normlz24=$av24/$Wild_PortoSanto_1 if $Wild_PortoSanto_1;


		#print OUT2 $chromo,"\t",$temp,"\t",$coordinate,"\t",$Locus,":",$counter,"\t",$normlz1,"\t",$normlz2,"\t",$normlz3,"\t",$normlz4,"\t",$normlz5,"\t",$normlz6,"\t",$normlz7,"\t",$normlz8,"\t",$normlz9,"\t",$normlz10,"\t",$normlz11,"\t",$normlz12,"\t",$normlz13,"\t",$normlz14,"\t",$normlz15,"\t",$normlz16,"\t",$normlz17,"\t",$normlz18,"\t",$normlz19,"\t",$normlz20,"\t",$normlz21,"\t",$normlz22,"\t",$normlz23,"\t",$normlz24,"\n";
		$Depth_for_Domestic_BelgianHare_1_sum=0;$Depth_for_Domestic_Champagne_1_sum=0;$Depth_for_Domestic_Dutch_1_sum=0;$Depth_for_Domestic_FlemishGiant_1_sum=0;$Depth_for_Domestic_FrenchLop_1_sum=0;$Depth_for_Domestic_NetherlandDwarf_1_sum=0;$Depth_for_Domestic_NewZealand_1_sum=0;$Depth_for_Inbred__1_sum=0;$Depth_for_Lepusamericanus_outgroup_1_sum=0;$Depth_for_WildFrench_Caumont_sum=0;$Depth_for_WildFrench_LaRoque_1_sum=0;$Depth_for_WildFrench_Villemolaque_1_sum=0;$Depth_for_WildIberian_CO_1_sum=0;$Depth_for_WildIberian_Calzada_1_sum=0;$Depth_for_WildIberian_Carrion_1_sum=0;$Depth_for_WildIberian_Castanar_1_sum=0;$Depth_for_WildIberian_Guadalajara_1_sum=0;$Depth_for_WildIberian_Huelva_sum=0;$Depth_for_WildIberian_M_1_sum=0;$Depth_for_WildIberian_SCMora_1_sum=0;$Depth_for_WildIberian_TO_1_sum=0;$Depth_for_WildIberian_Toledo_1_sum=0;$Depth_for_WildIberian_Zaragoza_1_sum=0;$Depth_for_Wild_PortoSanto_1_sum=0;
	
		$counter=0;
		$temp=$coordinate;
	}
	else
	{
		
		
	}
}

close OUT;
close OUT2;




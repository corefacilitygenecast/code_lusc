#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Statistics::Basic qw(:all);
use File::Basename qw(basename dirname);
use Cwd;
use Cwd 'abs_path';
my $BEGIN_TIME=time();
my $version="1.1.0";
##############################################################
my ($minus_vaf,$max_vaf,$minus_depth,$bam,$exon_bed,$black_file,$type,$database,$prefix,$outdir,$snv,$driver,$config_file);
GetOptions(
	"help|?" =>\&usage,
	"config-file=s" => \$config_file,
	"type:s"=>\$type,
	"minux_af:s"=>\$minus_vaf,
	"max_af:s"=>\$max_vaf,
	"d:s"=>\$minus_depth,
	"bam:s" =>\$bam,
	"bed:s" =>\$exon_bed,
	"rank:s"=>\$database,
	"p:s"=>\$prefix,
	"od:s"=>\$outdir,
	"snv:s"=>\$snv,
	"black:s"=>\$black_file,
	"driver:s"=>\$driver
);

use Data::Dumper;
die &usage if (!$bam || !$database || !$snv ||!$outdir || !defined $config_file);
my ($low,$high,$extract_high ) = (0,0,0);
my %Driver;
my %pos;

open D, $driver or die $!;
while (<D>) {
	chomp;s/\r//g;
#	$Driver{$_}=1;
	my @s = split/\t/ ;
	$Driver{$s[0]}{$s[5]} = 1;
}close D;
################################################
$type ||= "tumor";
if($type eq 'tumor'){
	$minus_depth||= 100;
	$minus_vaf||=5;
	$max_vaf ||=100;
	($low,$high,$extract_high ) = (3.804,8.4278,13.7693);
}elsif($type eq "plasma"){
	$minus_depth ||= 500;
	$minus_vaf ||=0.8;
	$max_vaf ||=100;
	($low,$high,$extract_high ) = (2.4728,6.3721,12.99);
}
my $config_hash = &load_config_file($config_file);
my $depth_file = "$outdir/$prefix\_depth.txt";
my $cmd = ${$config_hash}{"SAMTOOLS"}." depth -d 1000000 -b $exon_bed $bam >$depth_file";
if(-e $depth_file){
	if(-z $depth_file){
		system($cmd);
	}
}else{
	system($cmd);
}
my %Hit_List = ();
if($black_file){
	open I,$black_file;
	while (<I>) {
		chomp;next if(/^chromosome/);
		next if(/^\#/);
		my @s = split/\t/;
		$Hit_List{$s[0]}{$s[1]} = 1 ;
	}
}
read_table($snv);

##################read SNV_File#################################
sub read_table {
	my $iput = shift;
	my ($key, $tab);
	my @exac= qw(ExAC_ALL ExAC_AFR ExAC_AMR ExAC_EAS ExAC_FIN ExAC_NFE ExAC_OTH ExAC_SAS gnomAD_exome_ALL gnomAD_exome_AFR gnomAD_exome_AMR gnomAD_exome_ASJ gnomAD_exome_EAS gnomAD_exome_FIN gnomAD_exome_NFE gnomAD_exome_OTH gnomAD_exome_SAS gnomAD_genome_ALL gnomAD_genome_AFR gnomAD_genome_AMR gnomAD_genome_ASJ gnomAD_genome_EAS gnomAD_genome_FIN gnomAD_genome_NFE gnomAD_genome_OTH 1000g2015aug_all 1000g2015aug_chb 1000g2015aug_chs 1000g2015aug_afr 1000g2015aug_eas 1000g2015aug_eur 1000g2015aug_sas 1000g2015aug_amr);
#	my @exac = qw(gnomAD_exome_ALL gnomAD_exome_AFR gnomAD_exome_AMR gnomAD_exome_ASJ gnomAD_exome_EAS gnomAD_exome_FIN gnomAD_exome_NFE gnomAD_exome_OTH gnomAD_exome_SAS);
	open IPUT, "<$iput" or die;
	my ($mean,$percentage20,$base_count) = parse_depth($depth_file);
	my $txt = <IPUT>;
	open O,">$outdir/$prefix.anno.hg19_multianno_filter.txt" or die $!;
	print O "Chr\tStart\tEnd\tRef\tAlt\tGENE\tAA\tratio\tdepth\n";
	my $final = basename($iput);
	$final =~ s/final.txt/final.addTMB.txt/;
	open F,">$outdir/$final" or die $!;
	chomp($txt);
	print F "$txt\tTMB_flag\n";
	$key = read_title ($txt);
	my $tmb_num = 0;
	while (<IPUT>) {
		/^#/ and next;
		/\w/ or next;
		chomp;
		my $tmp = read_value($key, $_);
		if($$tmp{'freq'}<$minus_vaf|| $$tmp{'freq'}>$max_vaf||$$tmp{'depth'}<$minus_depth){print F "$_\tNo\n";next ;}
		if($$tmp{'FILTER'} ne 'PASS'||$$tmp{'forward_reads'}==0 ||$$tmp{'reverse_reads'}==0){print F "$_\tNo\n";next ;}
		if($$tmp{'Func_refGene'} ne 'exonic'||exists $Driver{$$tmp{'Gene_refGene'}}{$$tmp{'mutation_p'}} || $$tmp{'black_list_flag'} ne ''){print F "$_\tNo\n";next ;}
		if($$tmp{'ExonicFunc_refGene'}!~/^nonsynonymous_SNV/ && $$tmp{'ExonicFunc_refGene'}!~/^frameshift/ && $$tmp{'ExonicFunc_refGene'}!~/^stopgain/){print F "$_\tNo\n";next ;}
		if(exists $Hit_List{$$tmp{'chrom'}}{$$tmp{'pos_raw'}}){print F "$_\tNo\n";next ;}####Filter Hit List SNV
		my $flag = 0;
		for my $index (@exac) {
			if($$tmp{$index} ne '.' && $$tmp{$index}>0.002&& $$tmp{$index} !~/,/){
				$flag++;
			}elsif($$tmp{$index} =~ /,/){
				my @record = (split/,/,$$tmp{$index});
				for my $record(@record) {
					$flag++ if($record ne '.' && $record >0.002);
				}
			} #######Remove larger than 0.002
		}
		if($$tmp{'ExonicFunc_refGene'}=~/^frameshift/ && $$tmp{'freq'}<=3){print F "$_\tNo\n";next ;}
		if($$tmp{'ExonicFunc_refGene'}=~/^frameshift/){
			my $Rlen = length($$tmp{'ref'});
			my $Alen = length($$tmp{'alt'});
			if($Rlen == 2||$Alen==2){
				if ($$tmp{'freq'}<=10){print F "$_\tNo\n";next ;}
			}
		}###########inframeshift INDEL minus freq 4% && filter single Base INDEL
		if($flag >0){print F "$_\tNo\n";next ;}
		if($$tmp{'snp138'} ne '.' && $$tmp{'cosmic80'} !~/ID/) {print F "$_\tNo\n";next ;}###
		#next if($$tmp{'snp138'} ne '.' && $$tmp{'cosmic80'} eq '.');
		$tmb_num++;
		print F "$_\tYes\n";
		print O "$$tmp{'chrom'}\t$$tmp{'pos_raw'}\t$$tmp{'pos_raw'}\t$$tmp{'ref'}\t$$tmp{'alt'}\t$$tmp{'Gene_refGene'}\t$$tmp{'mutation_p'}\t$$tmp{'freq'}\t$$tmp{'depth'}\n";
	}
	close IPUT;close O;
	open T,">$outdir/$prefix.mub.txt" or die $!;
	print T "Sample_ID\ttissue\tAbsolute_Mutation_Count\tNum_Exonic_Bases_Coverage\tNormalized_MuB\taverage_depth\tpercentage20\tTMB_level\trank\n";
	my $Normalized_MuB = sprintf("%.4f",$tmb_num/$base_count*1000000);
	my $level =tmb_level($Normalized_MuB);
	my $rank = rank($database,$Normalized_MuB);
	print T "$prefix\t$type\t$tmb_num\t$base_count\t$Normalized_MuB\t$mean\t$percentage20\t$level\t$rank\n";
}
#####################################################
sub tmb_level{
	#$low,$high,$extract_high 
	my $normal = shift @_;
	if($normal<$low){
		return "Low";
	}elsif($normal<$high){
		return "Moderate";
	}elsif($normal<$extract_high){
		return "High";
	}else{
		return "Extra High"
	}
}
#######################################
sub rank {
	my ($file,$tmb) =  @_;
	open I,$file or die $!;
	my ($rank,$total) = (0,0);
	while (<I>) {
		chomp;next if(/^Normalized/);
		next if(/^$/||$_!~/$type/);
		my @item = split/\t/;
		$total++;
		$rank++ if($item[0]<=$tmb);
	}close I;
	my $r = sprintf("%.4f",$rank/$total);
	return $r;
}
####################################################
sub read_title {
	my $txt = shift;

	my $key;
	chomp $txt;
	$txt =~s/\r//g;
	$txt =~ s/ +//g;
	$txt =~ s/^\#//;
	@$key = split("\t", $txt);
	return $key;
}

#####################################################
sub read_hash {
	my ($iput, $hash) = @_;
	open IPUT, "<$iput" or die;
	my $txt;
	while ($txt=<IPUT>) {
		$txt =~ /^#/ or next;
		my $key = read_title ($txt);
		$txt = <IPUT>;
		my $tab = read_value ($key, $txt);
		$$hash{$_} = $$tab{$_} foreach (keys %$tab);
	}
	close IPUT;
}
#####################################################
sub read_value {
	my ($key, $txt) = @_;
	my $hash;
	chomp $txt;
	$txt=~s/\r//g;
#	$txt =~ s/ //g;
	my @val = split "\t", $txt;
	for (my $i=0; $i<@{$key}; $i++) {
		$$hash{$$key[$i]} = (defined $val[$i])? $val[$i]:'';
	}
	return $hash;
}
####################################################
sub parse_depth{
	my $file = shift @_;
	my @depth = (); my $base_num = 1;
	open I,$file or die $!;
	while (<I>) {
		chomp;next if(/^$/||/^\#/);
		my @item = split/\t/;
		push @depth,$item[2];
		next if($item[2] < $minus_depth);
		$pos{$item[0]}{$item[1]} = $item[2] ;
		$base_num++;
	}close I;
	my $mean = mean(@depth);
	$mean=~s/,//g;
	my $percentage20 = grep $_>0.2*$mean,@depth;
	$percentage20 = sprintf("%.4f",$percentage20/@depth);
	return($mean,$percentage20,$base_num);
}

####################################################
sub usage{
	my $info=<<INFO;
USAGE:  $0 <option>
	-config-file <file> the cupcake conf file
	-snv <file> Commercial snv final file.
	-bed <file>  Bed region file name, for exom and target region capture.
	-driver <file> Driver genesymbol to skip file 
	-bam <file> BAM file 
	-type <STR>  Sample tissue type ,tumor or plasma default tumor
	-minux_af <float> minus AF ,tumor default 7%,plasma default 1%     
	-max_af <float max AF ,default 0.4                                 
	-d <int> minus cov  depth ,tumor default 100x ,plasma default 500x 
	-rank <file> the TMB database for mub ranking
	-black <file>  SNV Hit List file
	-p <STR> output file prefix
	-od     output directory
INFO
        #die $info;
}

sub load_config_file {
        my $config_file = shift;
        open IN, $config_file or die $!;
        my %config_hash;
        while (<IN>){
                chomp;
                next if (m/^#/);
                next if (m/^\s*$/);
                my ($config_name,$config_path) = split /=/,$_;
                $config_hash{$config_name} = $config_path;
        }
        close IN;
        return \%config_hash;
}


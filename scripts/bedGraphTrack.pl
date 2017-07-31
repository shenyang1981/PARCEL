#!/usr/bin/env perl

#Generate bedgraph track from bam file. Strand information is maintained and visualized in different colors (red, blue) using IGV or UCSC genome browser.

use strict;
use warnings;

my $cpuNumber = $ENV{"CPUNUM"} || 8;
my $inputBed = $ARGV[0];

if(@ARGV<3){
	die "perl bedGraphTrack.pl inputbed1,inputbed2 outputprefix genomesizefile [strand] [normlized] [sort]\n";
}

$inputBed =~ s/,/ /g;

my $outputPrefix = $ARGV[1];

my $genomeSizeFile = $ARGV[2];
my $strand = $ARGV[3] || "no";
my $normalized = $ARGV[4]||"no";
my $sort = $ARGV[5]||"no";

my $only5p = $ENV{"FIVEPOS"} || "no";
my $cutoff = $ENV{"MINCOV"} || 0;
my $region = $ENV{"INTERVAL"} || "";
my $pairend = $ENV{"PAIREND"} || "no";

my $selectInt="";
my $cat="cat $inputBed";

if($inputBed =~ m/gz$/){
	$cat = "zcat $inputBed";
}

my $if5p = $only5p eq "no" ? "": "-5";
my $issplit = $only5p eq "no" ? "-split" : "";

if($inputBed =~ m/bam$/){
	if($region ne ""){
		$cat = "samtools view -h -b $inputBed $region | bedtools bamtobed $issplit -i -";
	} else {
	$cat = "bedtools bamtobed $issplit -i $inputBed";
	}
}

$pairend = $pairend eq "no" ? "" : "-pc";
my $inputfile = $inputBed;
my $fileFlag=0;
my $totalReads = 0;
foreach my $file (split(/ /,$inputfile)){
	$fileFlag++ if ! -s $file;
	if($normalized ne "no"){
		if($file=~ m/bam$/){
			$totalReads += `samtools idxstats $file | mawk '{sum+=\$3}END{print sum}'`;
		} else {
			$totalReads += `wc -l $file | mawk '{print \$1}'`;
		}
	}	
}
my $scaleFactor=defined($ENV{"SCALE"}) ? 1e6/$ENV{"SCALE"}:1;

if($normalized ne "no"){
	$scaleFactor=1e6/$totalReads;
}

my $sortcmd = $sort eq "no" ? "": qq{ sort -k1,1 -k2,2n --parallel $cpuNumber -S 4G | };
my $cmd;
my $cmd1;
my $filterCmd = $cutoff ==0 ? "":qq{ | mawk '{OFS="\\t"}{if(\$1~ /^track/ || \$4>=$cutoff || \$4<= $cutoff*-1){print;}}'};

if($strand eq "no"){
	if($pairend ne ""){
		$cmd = qq{ bedtools genomecov $if5p $issplit -scale $scaleFactor -ibam $inputBed -g $genomeSizeFile -bg -trackline -trackopts 'name=\"$outputPrefix\" visibility=2 color=255,30,30' }.$filterCmd;
	} else {
		$cmd = qq{$cat | } .$sortcmd.qq{ bedtools genomecov $if5p $issplit -scale $scaleFactor -i - -g $genomeSizeFile -bg -trackline -trackopts 'name=\"$outputPrefix\" visibility=2 color=255,30,30' }.$filterCmd;

	}
	print STDERR $cmd."\n";
	system($cmd) if ! -f $outputPrefix."bedgraph.gz" && $fileFlag==0;
} else {
	my $pairendstrand = $pairend eq "" ? "" : "-du";
	if($pairendstrand eq "-du"){
		$cmd = qq{ bedtools genomecov $pairendstrand $if5p $issplit -scale $scaleFactor -strand "+" -ibam $inputBed -g $genomeSizeFile -bg -trackline -trackopts 'name=\"$outputPrefix\" visibility=2 color=255,30,30'}.$filterCmd;
		$cmd1 = qq{ bedtools genomecov $pairendstrand $if5p $issplit -scale $scaleFactor -strand "-" -ibam $inputBed -g $genomeSizeFile -bg  | mawk '{OFS="\\t"}{print \$1,\$2,\$3,\$4*-1;}'}.$filterCmd;
	} else {
	$cmd = qq{$cat | } .$sortcmd.qq{ bedtools genomecov $if5p $issplit -scale $scaleFactor -strand "+" -i - -g $genomeSizeFile -bg -trackline -trackopts 'name=\"$outputPrefix\" visibility=2 color=255,30,30'}.$filterCmd;
	$cmd1 = qq{$cat | } .$sortcmd.qq{ bedtools genomecov $if5p $issplit -scale $scaleFactor -strand "-" -i - -g $genomeSizeFile -bg  | mawk '{OFS="\\t"}{print \$1,\$2,\$3,\$4*-1;}'}.$filterCmd;
	}
	print STDERR $cmd."\n";
	print STDERR $cmd1."\n";
	system($cmd) if $fileFlag==0;
	system($cmd1) if $fileFlag==0;
}


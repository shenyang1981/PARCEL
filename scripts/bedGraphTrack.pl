use strict;
use warnings;

my $cpuNumber = $ENV{"CPUNUM"} || 8;
my $inputBed = $ARGV[0];

if(@ARGV<3){
	die "perl bedGraphTrack.pl inputbed1,inputbed2 outputprefix speciescode [strand] [normlized] [sort]\n";
}

$inputBed =~ s/,/ /g;

my $outputPrefix = $ARGV[1];

my $genome = $ARGV[2];
my $strand = $ARGV[3] || "no";
my $normalized = $ARGV[4]||"no";
my $sort = $ARGV[5]||"no";

my $only5p = $ENV{"FIVEPOS"} || "no";
my $cutoff = $ENV{"MINCOV"} || 0;
my $region = $ENV{"INTERVAL"} || "";
my $pairend = $ENV{"PAIREND"} || "no";
my $genomeList={
	"HMPREF_Gut" => "/home/sheny/prog/database/Genome/UCSC/HMPREF/gastrointestinal_tract/allchr.sizes",
	"AllBac" => "/home/sheny/prog/database/Genome/UCSC/AllBacteria/allchr.chrom.sizes",
	"Eightbac" => "/home/sheny/prog/database/Genome/UCSC/EightBacteriaPools/allchr.chrom.sizes",
	"Yeast" => "~/prog/database/Genome/UCSC/Yeast/sacCer3/chrom.genomes",
	"YeastEnsembl" => "~/prog/database/Genome/Ensembl/S.cer/R64/chrom.sizes",
	"Dm6" => "/home/sheny/prog/database/Genome/Flybase/Dmel/FB2016_02/chrom.size",
	"Dm6TrL" => "/home/sheny/prog/database/Genome/Flybase/Dmel/FB2016_02/transcriptome_LongestTr.size",
	"Dm6Tr" => "/home/sheny/prog/database/Genome/Flybase/Dmel/FB2016_02/transcriptome.size",
	"Pombe" => "/home/sheny/prog/database/Genome/UCSC/Yeast/asm294/chrom.sizes",
	"Bsubtilis"=>"/home/sheny/prog/database/Genome/UCSC/EightBacteriaPools/allchr.chrom.sizes",
	"Paer"=>"/home/sheny/prog/database/Genome/UCSC/EightBacteriaPools/allchr.chrom.sizes",
	"Candida"=>"~/prog/database/Genome/UCSC/C.albicans/A21/allchr.sizes",
	"CandidaTr"=>"/home/sheny/prog/sunm/others/candida2/ref/C_albicans_genes_Snyder_UTRs.fa.size",
	"mm10"=> "~/prog/database/Genome/UCSC/Mouse/mm10/mm10.chrom.sizes",
	"mm10Rep" => "~/prog/database/repBase/repeatmaskerlibraries-20140131/mm10Rep/mm10_rmsk_repGenome.chrom.sizes",
	"mm10CDS"=> "~/prog/database/Genome/Ensembl/Mouse/GRCm38/Mus_musculus.GRCm38.cds.all.longest.Length.txt",
	"rn5"=> "~/prog/database/Genome/UCSC/Rat/rn5/rn5.chrom.sizes",
	"hg19"=> "~/prog/database/Genome/UCSC/Human/hg19/hg19.chrom.sizes",
	"g1k"=> "/mnt/AnalysisPool/libraries/genomes/human_g1k_v37/human_g1k_v37.chrom.sizes",
	"repeatGenome" => "~/prog/database/repBase/repeatmaskerlibraries-20140131/chrRep.sizes",
	"repeatMouse" =>"~/prog/database/repBase/repeatmaskerlibraries-20140131/RepeatMaskerLib_Mouse.length.txt",
	"subZfp281Con" => "~/gseq/prog/Chengqi/metaAnalysis/L1evolution/fullL1/subZfp281.con.cls.length.txt",
	"repeat" => "~/prog/database/repBase/repeatmaskerlibraries-20140131/RepeatMaskerLib.Length.txt"};


if(!defined($genomeList->{$genome})){
	die "$genome is not found in genome DBs\n";
}

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
	#$cmd = qq{$cat | } .$sortcmd.qq{ bedtools genomecov -split -scale $scaleFactor -i - -g $genomeList->{$genome} -bg -trackline -trackopts 'name=\"$outputPrefix\" visibility=2 color=255,30,30' }.$filterCmd.qq{ | pigz -p $cpuNumber > $outputPrefix.bedgraph.gz};
	if($pairend ne ""){
		$cmd = qq{ bedtools genomecov $if5p $issplit -scale $scaleFactor -ibam $inputBed -g $genomeList->{$genome} -bg -trackline -trackopts 'name=\"$outputPrefix\" visibility=2 color=255,30,30' }.$filterCmd;
	} else {
		$cmd = qq{$cat | } .$sortcmd.qq{ bedtools genomecov $if5p $issplit -scale $scaleFactor -i - -g $genomeList->{$genome} -bg -trackline -trackopts 'name=\"$outputPrefix\" visibility=2 color=255,30,30' }.$filterCmd;

	}
	print STDERR $cmd."\n";
	system($cmd) if ! -f $outputPrefix."bedgraph.gz" && $fileFlag==0;
} else {
	my $pairendstrand = $pairend eq "" ? "" : "-du";
	if($pairendstrand eq "-du"){
		$cmd = qq{ bedtools genomecov $pairendstrand $if5p $issplit -scale $scaleFactor -strand "+" -ibam $inputBed -g $genomeList->{$genome} -bg -trackline -trackopts 'name=\"$outputPrefix\" visibility=2 color=255,30,30'}.$filterCmd;
		$cmd1 = qq{ bedtools genomecov $pairendstrand $if5p $issplit -scale $scaleFactor -strand "-" -ibam $inputBed -g $genomeList->{$genome} -bg  | mawk '{OFS="\\t"}{print \$1,\$2,\$3,\$4*-1;}'}.$filterCmd;
	} else {
	$cmd = qq{$cat | } .$sortcmd.qq{ bedtools genomecov $if5p $issplit -scale $scaleFactor -strand "+" -i - -g $genomeList->{$genome} -bg -trackline -trackopts 'name=\"$outputPrefix\" visibility=2 color=255,30,30'}.$filterCmd;
	$cmd1 = qq{$cat | } .$sortcmd.qq{ bedtools genomecov $if5p $issplit -scale $scaleFactor -strand "-" -i - -g $genomeList->{$genome} -bg  | mawk '{OFS="\\t"}{print \$1,\$2,\$3,\$4*-1;}'}.$filterCmd;
	}
	print STDERR $cmd."\n";
	print STDERR $cmd1."\n";
	system($cmd) if $fileFlag==0;
	system($cmd1) if $fileFlag==0;
}


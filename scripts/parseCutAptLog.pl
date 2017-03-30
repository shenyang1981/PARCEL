use strict;
use warnings;

use IO::Handle;
use IO::File;

my $fh;

if(@ARGV==0){
	$fh = IO::Handle->new();
	$fh->fopen(fileno(STDIN),"r");
} else {
	$fh = IO::File->new();
	$fh->open($ARGV[0],"r");
}

my $id = $ENV{"ID"} || "All";
my $totalReads=0;
my $writtenReads=0;
my $adapterReads=0;
while(my $line=$fh->getline()){
	if($line=~ m/Total reads processed:(.*)/){
		$totalReads=$1;
		$totalReads=~ s/\s//sg;
		$totalReads=~ s/,//sg;
	}
	if($line=~ m/Reads with adapters:(.*)\(/){
		$adapterReads = $1;
		$adapterReads =~ s/\s//sg;
		$adapterReads =~ s/,//sg;
	}
	if($line=~ m/Reads written.*:(.*)\(/){
		$writtenReads = $1;
		$writtenReads =~ s/\s//sg;
		$writtenReads =~ s/,//sg;
	}
}
print $id."\t".$totalReads."\t".$writtenReads."\t".$adapterReads."\n";



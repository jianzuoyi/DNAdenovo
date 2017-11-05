#!/usr/bin/perl
use strict;
use Getopt::Long;
my $k;
GetOptions(
	"k:i"=>\$k,
);
my $histo=$ARGV[0].".histo";
my $stats=$ARGV[0].".stats";

print STDERR "The following computation only can be used for the normal k-mer distribution of genome reads with depth>20!\n";
my $peak=`sed '1d;2d;3d;4d;5d;6d;7d;8d;9d;10d;11d;12d;13d;14d;15d;16d;17d;18d;19d;\$D' $histo | sort -nr -k 2 | head -1 | cut -d' ' -f1`+1;
my $tmp=`cat $stats`;
$tmp =~ m/\D*(\d+)\D*(\d+)\D*(\d+)\D*(\d+)/;
my ($unique,$distinct,$total,$max_count)=($1,$2,$3,$4);
my $genome_size=int(($total-$unique)/$peak);
print "kmer:        $k\nmer_num:     $total\nnum_1:       $unique\npeak:        $peak\ngenome_size: $genome_size\nnode_num:    $distinct\n";

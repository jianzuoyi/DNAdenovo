#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ( $low_quality_p, $N_p, $Q_cut, $short ) = ( 40, 20, 20, 0 );
my $help;

GetOptions(
    "Q:i" => \$Q_cut,
    "B:i" => \$low_quality_p,
    "w:i" => \$N_p,
    "s:i" => \$short,
    "h"   => \$help
);

my $usage = "
NAME:
	FAS - Filter_After_Sffinfo

USAGE:
	perl $0 [Options] <qual_file> <seq_file>

OPTIONS:
	-Q <int>  Base with less value than this is defined as low quality bases, default 20
	-B <int>  filter reads with >X percent base are low quality bases,set a cutoff, default 40
	-w <int>  filter reads with >X percent base are Ns,set a cutoff, default 20
	-s <int>  filter reads with length shorter than this threshold value, default 0 nucleotides.
	-h	  print this help information

NOTE:
	The -B and -w parameters is same as Solexa fastaq data filter pipeline

EXAMPLE:
	perl $0 raw.1.qual raw.1.fa >filter.1.fa 2>filter.1.id

AUTHOR:
	DENG Cao
	dengcao\@bgitechsolutions.com
	brentcaodeng\@hotmail.com
	2012-08-13
";
die print $usage if $help;
( @ARGV && ( @ARGV % 2 == 0 ) ) || die print $usage;

my $qual  = shift @ARGV;
my $fasta = shift @ARGV;
my %hash;

# filter the low_quality and too short reads
open( QUAL, $qual ) or die "Cannot Open File: $qual\n";
$/ = "\n>";
while (<QUAL>) {
    chomp;
    s/^>//;
    next if /^\s+$/;
    s/(^\s+)|(\s+$)//;
    my ( $ID, $q_values ) = split /\n/, $_, 2;
    ($ID) = $ID =~ /^(\S+)/;
    my @values = split( /\s+/, $q_values );
    next if @values < $short; ## filter the too short reads
    my $low_qual = 0;

    for ( my $i = 0 ; $i < @values ; $i++ ) {
        $low_qual++ if ( $values[$i] < $Q_cut );
    }
    my $ratio = $low_qual / @values * 100;
    if ( $ratio < $low_quality_p ) { $hash{$ID} = 1; }
}
close QUAL;

# filter the Ns
open( FA, $fasta ) or die "Cannot Open File: $fasta\n";
while (<FA>) {
    chomp;
    s/^>//;
    next if /^\s+$/;
    s/\s+$//;
    my ( $head, $seq_raw ) = split /\n/, $_, 2;
    my ($id) = $head =~ /^(\S+)/;
    next unless exists $hash{$id};
    my $seq = $seq_raw;
    $seq =~ s/\s+//g;
    my $N = $seq =~ tr/N/n/;
    my $ratio = $N / ( length $seq ) * 100;
    next if $ratio > $N_p;
    # print filtered data
    print STDERR "$id\n";
    print ">$head\n$seq_raw\n";
}
close FA;
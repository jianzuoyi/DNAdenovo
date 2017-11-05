#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $len_col = 2;
GetOptions( "len_col:i" => \$len_col );
@ARGV == 1 || die "
Usage: perl $0 [Options] length_ file > length_distribution.xls
-len_col <int> default 2;
";

my %len;
open IN, $ARGV[0];
while (<IN>) {
    next if /^#/ || /^\s+$/;
    my ($len) = ( split /\s+/ )[ $len_col - 1 ];
    $len{$len}++;
}
close IN;

foreach my $len ( sort { $a <=> $b } keys %len ) {
    print "$len\t$len{$len}\n";
}

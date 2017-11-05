#!/usr/bin/perl -w
use strict;
@ARGV == 1 || die "Usage: perl $0  quality_ file >quality_distribution.xls";

my %qual;
open IN, $ARGV[0];
while (<IN>) {
    next if /^#/ || /^\s+$/;
    next if (/^>/);
    s/^\s+//;
    s/\s+$//;
    my @quals = split;
    for ( my $i = 0 ; $i < @quals ; $i++ ) {
        $qual{ $quals[$i] }++;
    }
}

foreach ( sort { $a <=> $b } keys %qual ) {
    print "$_\t$qual{$_}\n";
}

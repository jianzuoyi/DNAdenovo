#!/usr/bin/perl -w
use strict;

@ARGV == 1 || die "
Usage: perl $0 <seq_length_infor> >result.xls
";

my ( %infor, $total_read, $total_base );
open IN, $ARGV[0];
while (<IN>) {
    my ( $sample, $reads, $length ) = split /\s+/;
    $infor{$sample}{'read_num'}++;
    $infor{$sample}{'read_len'} += $length;
    $infor{$sample}{'read_max'} ||= $length;
    $infor{$sample}{'read_min'} ||= $length;
    $infor{$sample}{'read_max'} = $length if ( $infor{$sample}{'read_max'} < $length );
    $infor{$sample}{'read_min'} = $length if ( $infor{$sample}{'read_min'} > $length );
    $total_read++;
    $total_base += $length;
}
close $ARGV[0];

print "Sample\tRead Number\tPercent of Reads (%)\tTotal Length(bp)\tPercent of Bases (%)";
print "\tMin Length(bp)\tMax Length(bp)\tAverage Length(bp)\n";
foreach my $sample ( sort keys %infor ) {
    print "$sample\t$infor{$sample}{'read_num'}\t";
    printf "%.2f", ( $infor{$sample}{'read_num'} / $total_read * 100 );
    print "\t$infor{$sample}{'read_len'}\t";
    printf "%.2f", ( $infor{$sample}{'read_len'} / $total_base * 100 );
    print "\t$infor{$sample}{'read_min'}\t$infor{$sample}{'read_max'}\t";
    print int( $infor{$sample}{'read_len'} / $infor{$sample}{'read_num'} );
    print "\n";
}

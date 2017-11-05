#!/usr/bin/perl -w
use strict;
@ARGV == 3 || die "
Usage:
	perl $0 <gff> <quality_in> <quality_out> 2> running.log
Author:
	DENG Cao
	2012-09-03
	brentcaodeng\@hotmail.com
	dengcao\@bgitechsolutions.com
";
##=============##
## Processing  ##
##=============##
my ( $gff, $in, $out ) = @ARGV;
my $infor;
open GFF, $gff;
while (<GFF>) {
    next if /^#/;
    next if /^\s+$/;
    my ( $id, $full_length, $start, $end ) = ( split /\s+/ )[ 0, 1, 2, 3 ];
    $infor->{$id}->{'full_length'} = $full_length;
    $infor->{$id}->{'start'}       = $start - 1;     # change it to start from 0
    $infor->{$id}->{'end'}         = $end - 1;       # change it to start from 0
}
close GFF;
open IN,  $in;
open OUT, ">$out";
$/ = "\n>";
while (<IN>) {
    chomp;
    s/^>//;
    next if /^\s+$/;
    s/(^\s+)|(\s+$)//;
    my ( $head, $quals ) = split /\n/, $_, 2;
    my ($id) = $head =~ /^(\S+)/;
    next unless exists $infor->{$id}; ## escape the filtered id;
    my @quals = split /\s+/, $quals;
    my $real = scalar @quals;

    if ( $infor->{$id}->{'full_length'} == @quals ) {
        my ( $start, $end ) = ( $infor->{$id}->{'start'}, $infor->{$id}->{'end'} );
        print OUT ">$id\n", join( " ", @quals[ $start .. $end ] ), "\n";
    }
    else {
        Message( "The length of $id in $gff is $infor->{$id}->{'full_length'}", \*STDERR );
        Message( "The length of $id in $in is $real", \*STDERR );
        Message( "Please check the input files",      \*STDERR );
    }
}
$/ = "\n";
close OUT;
close IN;
##=============##
## Subroutines ##
##=============##
sub Message {
    ## Message($infor,\*HANDLE)
    my ( $infor, $fh ) = @_;
    my @new = split /\n+/,$infor;
    if ( defined fileno $fh ) { ## STDIN = 0, STDOUT = 1, STDERR = 2, others >=3
        foreach (@new){
            print $fh Date() . "\t\t$_\n";
	}
	print $fh "\n";
    }
    else {
        $fh = \*STDERR;
        print $fh Date() . "\t\tThe parameter in Message() is not a filehandle\n";
	foreach (@new){
            print $fh Date() . "\t\t$_\n";
	}
	print $fh "\n";
    }
}

sub Date {
    my $time = `date`;
    $time =~ s/^\s+|\s+$//;
    return '[' . $time . ']';
}

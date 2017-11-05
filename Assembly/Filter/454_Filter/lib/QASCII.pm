package QASCII;
use strict;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(Sanger_2QA Sanger_2ASCII Solexa_2QA Solexa_2ASCII);

sub Sanger_2QA {
    my ($ascii) = @_;
    my @quals;
    my @ascii = split //, $ascii;
    foreach my $chr (@ascii) {
        push @quals, ( ord($chr) - 33 );    ## QA = ASCII - 33;
    }
    return @quals;
}

sub Sanger_2ASCII {
    my (@qual) = @_;
    my $qual;
    foreach my $num (@qual) {
        $num += 33;
        $qual .= chr $num;                  ## ASCII = 33 + QA;
    }
    return $qual;
}

sub Solexa_2QA {
    my ($ascii) = @_;
    my @quals;
    my @ascii = split //, $ascii;
    foreach my $chr (@ascii) {
        push @quals, ( ord($chr) - 64 );    ## QA = ASCII - 33;
    }
    return @quals;
}

sub Solexa_2ASCII {
    my (@qual) = @_;
    my $qual;
    foreach my $num (@qual) {
        $num += 64;
        $qual .= chr $num;                  ## ASCII = 33 + QA;
    }
    return $qual;
}

1;

#!/usr/bin/perl -w
use strict;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use QASCII;
use Getopt::Long;

my $usage = "
Usage:
	perl $0 [Options] <file_prefix.fq|file_prefix.fa>

Author:
	DENG Cao
	2012-09-03
	brentcaodeng\@hotmail.com
	dengcao\@bgitechsolutions.com

Options:
	-manner <Sanger|Solexa>	Choose the transcoding manner between QA and ASCII, default Solexa
				if you choose Sanger, ASCII = QA + 33; if you choose Solexa, ASCII = QA + 64;
	-fq2fa			Change FQ file to FASTA and QUAL files, and the input file must be FQ format, this is default one
	-fa2fq			Change FASTA and QUAL files to FQ file, and the input file must be FASTA format and -qual must be setted
	-quality <str>		if you choose -fa2fq, this parameter must be setted to give this script QUAL file
	-out_pre <str>		default 'file_prefix.fq' or 'file_prefix.fa'
	-out_dir <str>		default .
	-help|man|?		print this help information
Note:
	Once you use -fq2fa, you cannot use -fa2fq and -quality and the input file muse be FQ format;
	Once you use -fa2fq, you must use -quality and the input file muse be FASTA format. Moreover, -fq2fa is forbidden.

Example:
	perl $0 file_prefix.fq 2>log_file
	perl $0 -fa2fq -qual file_prefix.qa file_prefix.fa 2>log_file
	perl $0 -fa2fq -qual file_prefix.qa -manner Sanger -out_pre test file_prefix.fa 2>log_file

";
##====================##
##   get  parameters  ##
##====================##
my ( $manner, $fq2fa, $fa2fq, $qual, $out_pre, $out_dir, $help );
GetOptions(
    "manner:s"   => \$manner,
    "fq2fa"      => \$fq2fa,
    "fa2fq"      => \$fa2fq,
    "quality:s"  => \$qual,
    "out_pre:s"  => \$out_pre,
    "out_dir:s"  => \$out_dir,
    "help|man|?" => \$help
);
die print $usage if @ARGV != 1 || $help;
##====================##
## default parameters ##
##====================##
$out_dir ||= '.';
die unless MakeDir( $out_dir, 0 );

unless ($out_pre) {
    $out_pre = basename( $ARGV[0] );
    if ($out_pre =~ /\.(fa|fas|fasta|fq|fastaq)(.{0}|\.gz)$/i){
        $out_pre =~ s/\.(fa|fas|fasta|fq|fastaq)(.{0}|\.gz)$//i;
    }
}
$manner ||= 'Solexa';
$fq2fa = 1 unless $fq2fa || $fa2fq;

##====================##
##  check  parameters ##
##====================##
my $in = shift @ARGV;

if ($fq2fa) {
#    if ( $in !~ /\.(fq|fastaq)(.{0}|\.gz)$/i ) {
#        Message( "Your infile is not FQ format", \*STDERR );
#    }
    if ( $fa2fq || $qual ) {
        Message( "You defined -fq2fa & -fa2fq/-qual synchronously", \*STDERR );
        die;
    }
}

if ($fa2fq) {
#    if ( $in !~ /\.(fa|fas|fasta)(.{0}|\.gz)$/i ) {
#        Message( "Your infile is not FASTA format", \*STDERR );
#    }
    unless ($qual) {
        Message( "You defined -fa2fq so you must define -qual", \*STDERR );
        die;
    }
}

unless ( $manner =~ /^Solexa$/i or $manner =~ /^Sanger$/i ) {
    Message( "Unknown transcoding manner: $manner", \*STDERR );
    die;
}

##=============##
##    fa2fq    ##
##=============##
if ($fa2fq) {

    # open filehandle
    if ( $in =~ /\.gz$/i ) {
        open FA, "gunzip -dc $in |" || die print STDERR "Cannot open $in\n";
    }
    else {
        open FA, $in || die print STDERR "Cannot open $in\n";
    }

    if ( $qual =~ /\.gz$/i ) {
        open QUAL,
          "gunzip -dc $qual|" || die print STDERR "Cannot open $qual\n";
    }
    else {
        open QUAL, $qual || die print STDERR "Cannot open $qual\n";
    }

    open OUTFQ, ">$out_dir/$out_pre.fq";
    $/ = "\n>";

    # reading quality file
    my $quality;
    while (<QUAL>) {
        chomp;
        s/^>//;
        next if /^\s+$/;
        s/(^\s+)|(\s+$)//;
        my ( $head, $quals ) = split /\n/, $_, 2;
        my ($id) = $head =~ /^(\S+)/;
        unless ( $id && $quals ) {
            Message( "Your file is not exactly FASTA format", \*STDERR );
            die;
        }
        my @quals = split /\s+/, $quals;
        my $q;
        $q = Sanger_2ASCII(@quals) if $manner =~ /^Sanger$/i;
        $q = Solexa_2ASCII(@quals) if $manner =~ /^Solexa$/i;
        $quality->{$id} = $q;
    }

    # reading fasta file and output result
    while (<FA>) {
        chomp;
        s/^>//;
        next if /^\s+$/;
        s/(^\s+)|(\s+$)//;
        my ( $head, $seq ) = split /\n/, $_, 2;
        my ($id) = $head =~ /^(\S+)/;
        unless ( $id && $seq ) {
            Message( "Your file is not exactly FASTA format", \*STDERR );
            die;
        }
        $seq =~ s/\s+//g;
        if ( !exists $quality->{$id} ) {
            Message( "The $id donot exist quality records", \*STDERR );
            die;
        }
        if ( ( length $seq ) != ( length $quality->{$id} ) ) {
            Message( "The sequence length of $id donot match to $quality->{$id}", \*STDERR );
	    my $length_seq = length $seq;
	    my $length_qual = length $quality->{$id};
	    Message("Length of seqience is $length_seq, and the sequence is:\n$seq\nLength of quality is $length_qual",\*STDERR);
            die;
        }
        print OUTFQ "\@$id\n$seq\n+\n$quality->{$id}\n";
    }
    $/ = "\n";
    close FA;
    close QUAL;
    close OUTFQ;
}
##=============##
##    fq2fa    ##
##=============##
if ($fq2fa) {

    # open filehandle
    if ( $in =~ /\.gz$/i ) {
        open FQ, "gunzip -dc $in |" || die print STDERR "Cannot open $in\n";
    }
    else {
        open FQ, $in || die print STDERR "Cannot open $in\n";
    }
    open OUTFA, ">$out_dir/$out_pre.fa";
    open OUTQA, ">$out_dir/$out_pre.qual";

    # reading fq file and output result

    $/ = "\n@";
    while (<FQ>) {
        chomp;
        s/^@//;
        next if /^\s+$/;
        s/(^\s+)|(\s+$)//;
        my ( $head, $seq, $dull, $chr ) = split /\n/, $_, 4;
        unless ( $head || $seq || $dull || $chr ) {
            Message( "Your file is not exactly FQ format", \*STDERR );
            die;
        }
        my ($id) = $head =~ /^(\S+)/;
        my @qa;
	@qa = Sanger_2QA($chr) if $manner =~ /^Sanger$/i;
        @qa = Solexa_2QA($chr) if $manner =~ /^Solexa$/i;
        if ( ( length $seq ) != @qa ) {
            Message( "The sequence length of $id donot match to quality length", \*STDERR );
            die;
        }
        print OUTFA ">$id\n$seq\n";
        print OUTQA ">$id\n", join( " ", @qa );
    }
    $/ = "\n";
    close FQ;
    close OUTFA;
    close OUTQA;
}
##=============##
## Subroutines ##
##=============##
sub MakeDir {
    ## $clean = 0/1, 1 means delete files if exists the DIR.
    ## it returns 0/1, 1 means Success while 0 fail.
    my ( $dir, $clean ) = @_;
    if ( $dir eq '.' ) {
        Message( "The dir is '.', no need to creat!", \*STDERR );
        return 1;
    }
    if ( -d $dir ) {
        Message( "Allready exists $dir", \*STDERR );
        if ($clean) {
            if ( system("rm -rf $dir/*") ) {
                Message( "Failed to delete file(s) in $dir", \*STDERR );
                return 0;
            }
            else {
                Message( "Successfully delete file(s) in $dir", \*STDERR );
                return 1;
            }
        }
	return 1;
    }
    else {
        if ( mkdir $dir ) {
            Message( "Successfully creat $dir", \*STDERR );
            return 1;
        }
        else {
            Message( "Failed to creat $dir", \*STDERR );
            return 0;
        }
    }
}

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
        print $fh Date() . "\t\t$fh is not a filehandle\n";
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

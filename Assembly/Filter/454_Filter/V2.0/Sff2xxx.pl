#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($RealBin);

my ( $seq, $qual, $flow, $tab, $notrim, $mft, $plain, $all, $lst, $dir, $help );

GetOptions(
    "seq"        => \$seq,
    "qual"       => \$qual,
    "flow"       => \$flow,
    "tab"        => \$tab,
    "notrim"     => \$notrim,
    "mft"        => \$mft,
    "plain"      => \$plain,
    "all"        => \$all,
    "lst:s"      => \$lst,
    "dir:s"      => \$dir,
    "help|?|mam" => \$help,
);

my $usage = "
NAME:
	Sff2xxx - Using sffinfo change sff to other format and get the qualilty files

USAGE:
	perl $0 [Oprionts] <sample_name> <sff_file> 2>running.log

OPTIONS:
	-s or -seq      Output just the sequences, default yes
	-q or -qual     Output just the quality scores, default default yes
	-f or -flow     Output just the flowgrams, default no
	-t or -tab      Output the seq/qual/flow as tab-delimited lines, default no
	-n or -notrim   Output the untrimmed sequence or quality scores, default no
	-m or -mft      Output the manifest text, default no
	-p or -plain	Output the plain text, default no
	-a or -all	Output the all above files, default no
	-l or -lst <s>	*.sff files list,if define this, you need not to use <sff_file> <sample_name> at commond line
	-d or -dir	Outdir of the output files. defualt .
	-h or -help	Print this information

NOTE:
	The format of *.sff files list used by -l or -lst parameter:
	# If defined this parameter, the output files will use the raw.Sample_Name as prefix
	# If the different sff file have same Sample_Name,they'll be included in one result file
	# Sample_Name	FIle_Path
	A1_1		../00.raw_data/454Reads.1.sff
	A1_2		../00.raw_data/454Reads.2.sff
	...
					  
AUTHOR:
	DENG Cao
	dengcao\@bgitechsolutions.com
	brentcaodeng\@hotmail.com
	2012-08-14
";

### default and check the parrameters
$dir ||= '.';
if ( $dir ne '.' ) {
    unless ( -d $dir ) {
        die Message ("Cannot creat dir <$dir> : $!") if system("mkdir -p $dir");
    }
}
( $seq, $qual, $flow, $tab, $notrim, $mft, $plain ) = ( 1, 1, 1, 1, 1, 1, 1 ) if $all;
( $seq, $qual ) = ( 1, 1 ) if ( !$seq && !$qual && !$flow && !$tab && !$notrim && !$mft && !$plain );
die print $usage if $help			# User's requir
  or ( defined $lst  && ( @ARGV % 2 ) != 0 )    # <sff_file> and <sample_name> not defined and -lst not defined too
  or ( !defined $lst && @ARGV != 2 );   	# Wrong commond line parameter number
### get and check the sff files
my %files;

if ( @ARGV == 2 ) {
    my ( $sample_name, $sff_file ) = @ARGV;
    unless ( -f $sff_file ) {
        Message("There's no $sff_file, Please check!");
        die;
    }
    $files{$sample_name}{$sff_file} = 1;
}
if ( defined $lst ) {
    open LST, $lst;
    while (<LST>) {
        next if /^\s+$/ or /^#/;
        my ( $sample, $path ) = /^(\S+)\s+(\S+)/;
        unless ( -f $path ) {
            print STDERR Date() , "\t\tThere's no $path file, Please check your sff files lst.\n";
            die;
        }
        if ( exists $files{$sample}{$path} ) {
            print STDERR Date() , "\t\tYou have two same lines in your sff files lst:\t$sample\t$path\n";
        }
        $files{$sample}{$path} = 1;
    }
    close LST;
}

### get result
foreach my $sample ( sort keys %files ) {
    if ( scalar keys %{ $files{$sample} } > 1 ) {
        print STDERR Date() , "\t\tYou have same Sample name <$sample> in your sff files list, and the path is different.\n";
        print STDERR Date() , "\t\tthry are: \t", join( " <-> ", keys %{ $files{$sample} } ), "\n";
        print STDERR Date() , "\t\tIf they belong to one sample, please IGNORE this paragraph, or please change\n";
        print STDERR Date() , "\t\ttheir sample name in your sff files list to be different.\n";
    }
    foreach my $path ( sort keys %{ $files{$sample} } ) {
        OneRun( $path, "$dir/raw.$sample" );
    }
}

### Subroutines
sub OneRun {
    my ( $sff, $out_pre ) = @_;
    `$RealBin/../opt/sffinfo -s $sff >> $out_pre.fa`       if $seq;
    `$RealBin/../opt/sffinfo -q $sff >> $out_pre.qual`     if $qual;
    `$RealBin/../opt/sffinfo -f $sff >> $out_pre.flow`     if $flow;
    `$RealBin/../opt/sffinfo -t $sff >> $out_pre.tab`      if $tab;
    `$RealBin/../opt/sffinfo -n $sff >> $out_pre.notrim`   if $notrim;
    `$RealBin/../opt/sffinfo -m $sff >> $out_pre.mainfest` if $mft;
    `$RealBin/../opt/sffinfo $sff >> $out_pre.txt`         if $plain;
}

sub Message {
    ## Message($infor)
    my ($infor) = @_;
    my @new = split /\n+/,$infor;
    foreach (@new){
        print STDERR Date() . "\t\t$infor\n";
    }
    print STDERR "\n";
}

sub Date {
    my $time = `date`;
    $time =~ s/(^\s+)|(\s+$)//;
    return '[' . $time . ']';
}

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($RealBin);

my ( $sff2xxx, $fas, $lst, $dir, $out_pre, $clean, $help );

GetOptions(
    "sff2xxx:s"  => \$sff2xxx,
    "fas:s"      => \$fas,
    "lst:s"      => \$lst,
    "dir:s"      => \$dir,
    "out_pre:s"  => \$out_pre,
    "clean"      => \$clean,
    "help|?|mam" => \$help,
);

my $usage = "
NAME:
	SSH - Super Sff Handler

VERSION:
	1.0 - 2012-08-14

USAGE:
	perl $0 [Oprionts] [<sample_name> <sff_file>]

OPTIONS:
	## Commom Parameters
	-sff2xxx <str>		parameters for $RealBin/Sff2xxx.pl, default -sff2xxx=\"-s -q\";
	-fas <str>		parameters for $RealBin/FAS.pl, default -fas=\"-Q 20 -B 40 -w 2\";
	-l or -lst <str>	*.sff files list, if define this, you need not to use <sff_file> <sample_name> at commond line
	-d or -dir		Outdir of the output files. defualt .
	-o or -out_pre <str>	prefix of out files,default 'filter'
	-h or -help		Print this information

	# parameters for $RealBin/Sff2xxx.pl, 
	# just list them in -sff2xxx as -sff2xxx=\"-m -t -q\";
	#	-s or -seq      Output just the sequences, default yes
	#	-q or -qual     Output just the quality scores, default default yes
	#	-f or -flow     Output just the flowgrams, default no
	#	-t or -tab      Output the seq/qual/flow as tab-delimited lines, default no
	#	-n or -notrim   Output the untrimmed sequence or quality scores, default no
	#	-m or -mft      Output the manifest text, default no
	#	-p or -plain	Output the plain text, default no
	#	-a or -all	Output the all above files, default no
	# parameters for $RealBin/FAS.pl,
	# just list them in -fas as -fas=\"-Q 30 -B 30\"
	#	-Q <int>  Base with less value than this is defined as low quality bases, default 20
	#	-B <int>  filter reads with >X percent base are low quality bases,set a cutoff, default 40
	#	-w <int>  filter reads with >X percent base are Ns,set a cutoff, default 2

NOTE:
	(1).The format of *.sff files list used by -l or -lst parameter:
		# If defined this parameter, the output files will use the (raw/filter).Sample_Name as prefix
		# If the different sff file have same Sample_Name,they'll be included in one result file
		# Sample_Name	File_Path
		A1_1		../00.raw_data/454Reads.1.sff
		A1_2		../00.raw_data/454Reads.2.sff
		...
					  
EXAMPLE:
	perl $0 Sample_A ./Sample_A.sff
	perl $0 -l sff_lst -d A1
	perl $0 -l sff_lst -d A1 -sff2xxx=\"-s -q -f\"
	perl $0 -l sff_lst -d A1 -sff2xxx=\"-s -q -f\" -fas=\"-Q 30 -B 30\";

AUTHOR:
	DENG Cao
	dengcao\@bgitechsolutions.com
	brentcaodeng\@hotmail.com
	2012-08-14
";

### default and check the parrameters
$sff2xxx ||= "-s -q";
$fas     ||= "-Q 20 -B 40 -w 2";
$dir     ||= '.';
$out_pre ||= "filter";
if ( $dir ne '.' ) {
    unless ( -d $dir ) {
        die "Cannot creat dir <$dir> : $!\n" if system("mkdir $dir");
    }
}
die print $usage if $help
  or ( defined $lst && ( @ARGV % 2 ) != 0 )    # <sff_file> and <sample_name> not defined and -lst not defined too
  or ( !defined $lst && @ARGV != 2 );
open LOG, ">$dir/SSH.log";

### get and check the sff files
my %files;

if ( @ARGV == 2 ) {
    my ( $sample_name, $sff_file ) = @ARGV;
    unless ( -f $sff_file ) {
        print LOG "There's no $sff_file, Please check!\n\n\n";
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
            print LOG "There's no $path file, Please check your sff files lst.\n\n\n";
            die;
        }
        if ( exists $files{$sample}{$path} ) {
            print LOG "You Have two same lines in your sff files lst:\n$sample\t$path\n\n\n";
        }
        $files{$sample}{$path} = 1;
    }
    close LST;
}

### get change sff file to fasta and quality files and then filter
foreach my $sample ( sort keys %files ) {
    if (scalar keys %{ $files{$sample} } > 1){
        print LOG "You have same Sample name <$sample> in your sff files list, and the path is different.\n";
        print LOG "thry are: \n\t",join( "\n\t", keys %{ $files{$sample} } ), "\nif they ";
        print LOG "belong to one sample, please IGNORE this paragraph, or please change their sample name in ";
        print LOG "your sff files list to be different.\n\n\n";
    }
    
    foreach my $path ( sort keys %{ $files{$sample} } ) {
	die print LOG "Canot execute the $RealBin/Sff2xxx.pl: $!\n\n\n" if system( "perl $RealBin/Sff2xxx.pl $sample $path -dir $dir $sff2xxx" );
    }
    die print LOG "Canot execute the $RealBin/FAS.pl: $!\n\n\n" if system( "perl $RealBin/FAS.pl $fas $dir/raw.$sample.qual $dir/raw.$sample.fa >$dir/filter.$sample.fa 2>$dir/filter.$sample.id");

    # get reads length
    die print LOG "Canot execute the $RealBin/../opt/fastaDeal.pl: $!\n\n\n" if system( "perl $RealBin/../opt/fastaDeal.pl -attribute id:len $dir/raw.$sample.fa > $dir/raw.$sample.len");

    die print LOG "Canot execute the $RealBin/../opt/fastaDeal.pl: $!\n\n\n" if system( "perl $RealBin/../opt/fastaDeal.pl -attribute id:len $dir/filter.$sample.fa > $dir/filter.$sample.len");

    die print LOG "Canot execute the $RealBin/../opt/fishInWinter.pl: $!\n\n\n" if system( "perl $RealBin/../opt/fishInWinter.pl --bformat table --fformat fasta $dir/filter.$sample.id $dir/raw.$sample.qual >$dir/filter.$sample.qual");
    # change fa+qual to fq
    die print LOG "Canot execute the $RealBin/Fa2Fq.pl$!\n\n\n" if system( "perl $RealBin/Fa2Fq.pl -fa2fq -qual $dir/filter.$sample.qual -manner Sanger -out_pre filter.$sample -out_dir $dir $dir/filter.$sample.fa 2>$dir/Fa2Fq.log");
    die print LOG "Canot execute the $RealBin/Fa2Fq.pl$!\n\n\n" if system( "perl $RealBin/Fa2Fq.pl -fa2fq -qual $dir/raw.$sample.qual -manner Sanger -out_pre raw.$sample -out_dir $dir $dir/raw.$sample.fa 2>$dir/Fa2Fq.log");
    # get all reads length
    die print LOG "Canot execute the awk: $!\n\n\n" if system( "awk '{print \"$sample\t\"\$1\"\t\"\$2}' $dir/raw.$sample.len >> $dir/raw.len");

    die print LOG "Canot execute the awk: $!\n\n\n" if system( "awk '{print \"$sample\t\"\$1\"\t\"\$2}' $dir/filter.$sample.len >> $dir/filter.len");
}
# get lenth distribution and draw
die print LOG "$!\n" if system( "perl $RealBin/length_distribution.pl -len_col 3 $dir/raw.len > $dir/raw.length_distribution.txt");
die print LOG "$!\n" if system( "perl $RealBin/length_distribution.pl -len_col 3 $dir/filter.len > $dir/filter.length_distribution.txt");
die print LOG "$!\n" if system( "perl $RealBin/../opt/line_diagram.pl $dir/raw.length_distribution.txt -x_title \"Length(bp)\" -y_title \"Contents\" > $dir/raw.length_distribution.svg");
die print LOG "$!\n" if system( "perl $RealBin/../opt/line_diagram.pl $dir/filter.length_distribution.txt -x_title \"Length(bp)\" -y_title \"Contents\" > $dir/filter.length_distribution.svg");
# get quality distribution and draw
die print LOG "$!\n" if system( "cat $dir/raw.*.qual > $dir/raw.qual");
die print LOG "$!\n" if system( "cat $dir/filter.*.qual > $dir/filter.qual");
die print LOG "$!\n" if system( "perl $RealBin/quality_distribution.pl $dir/raw.qual > $dir/raw.quality_distribution.txt");
die print LOG "$!\n" if system( "perl $RealBin/quality_distribution.pl $dir/filter.qual > $dir/filter.quality_distribution.txt");
die print LOG "$!\n" if system( "rm $dir/raw.qual $dir/filter.qual");
die print LOG "$!\n" if system( "perl $RealBin/../opt/line_diagram.pl $dir/raw.quality_distribution.txt -x_title \"Length(bp)\" -y_title \"Contents\" > $dir/raw.quality_distribution.svg");
die print LOG "$!\n" if system( "perl $RealBin/../opt/line_diagram.pl $dir/filter.quality_distribution.txt -x_title \"Length(bp)\" -y_title \"Contents\" > $dir/filter.quality_distribution.svg");
# get statistic result
die print LOG "$!\n" if system( "perl $RealBin/statistic.pl $dir/raw.len > $dir/raw.statistic.xls");
die print LOG "$!\n" if system( "perl $RealBin/statistic.pl $dir/filter.len > $dir/filter.statistic.xls");
# change file name
die print LOG "$!\n" if system( "tar -zcf $dir/$out_pre.raw.fa.tar.gz $dir/raw.*.fa " );
die print LOG "$!\n" if system( "tar -zcf $dir/$out_pre.raw.qual.tar.gz $dir/raw.*.qual " );

die print LOG "$!\n" if system( "tar -zcf $dir/$out_pre.filter.fa.tar.gz $dir/filter.*.fa " );
die print LOG "$!\n" if system( "tar -zcf $dir/$out_pre.filter.qual.tar.gz $dir/filter.*.qual " );

die print LOG "$!\n" if system( "tar -zcf $dir/$out_pre.raw.fq.tar.gz $dir/raw.*.fq " );

die print LOG "$!\n" if system( "tar -zcf $dir/$out_pre.filter.fq.tar.gz $dir/filter.*.fq " );
# Clean the data
if ($clean){
    die print LOG "$!\n" if system( "rm $dir/*.fa $dir/*.len $dir/*.qual $dir/*.id" );
}else{
    unless (-d "$dir/origin_data"){
	die print LOG "$!\n" if system( "mkdir $dir/origin_data" );
	die print LOG "$!\n" if system( "mv $dir/*.fa $dir/*.len $dir/*.qual $dir/*.id $dir/*.fq $dir/origin_data" );
    }
}
close LOG;


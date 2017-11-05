#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($RealBin);

my $usage = "
NAME:
	SSH - Super Sff Handler

VERSION:
	1.0 - 2012-08-14
	2.0 - 2012-09-04

AUTHOR:
	DENG Cao
	dengcao\@bgitechsolutions.com
	brentcaodeng\@hotmail.com
	2012-08-14

USAGE:
	perl $0 [Oprionts] [<sample_name> <sff_file>]

DESCRIPTION:
	This pipeline is used for 4-Primer Amplicon Tagging data filter;

	There are 2 pairs of primers that go into the PCR reaction, which consist of:
	1. the TS(Target-specific primers) and the CS(Consensus Sequence)
	2. 454 Adapter sequence (key), Barcode (specific to each DNA template) sequence and the CS.
	So, the PCR production will be:
	Adapter(key)->Barcode->CS1->TS1--> SEQUENCE [-->RC_TS2->RC_CS2->RC_Barcode->RC_Adapter(key)]
	Adapter(key)->Barcode->CS2->TS2--> SEQUENCE [-->RC_TS1->RC_CS1->RC_Barcode->RC_Adapter(key)]
	RC: Reverse Complete

	This pipeline will remove Adapter(key)->Barcode->CS (and/or RC_CS->RC_Barcode->RC_Adapter(key))
	from each sample.


OPTIONS:
	## Commom Parameters
	-sff2xxx <str>		parameters for HOME/bin/Sff2xxx.pl, default -sff2xxx=\"-s -q\";
	-remove_CS <str>	parameters for HOME/bin/remove_CS.pl, default -remove_CS=\"-cat yes -remove yes\";
	-fas <str>		parameters for HOME/bin/FAS.pl, default -fas=\"-Q 20 -B 40 -w 2 -s 0\";
	-fa2fq <str>		parameters for HOME/bin/Fa2Fq.pl,default -fa2fq=\"-manner Sanger\";
	-l or -lst <str>	*.sff files list, if define this, you need not to use <sff_file> <sample_name> at commond line;
	-p or -project_id <str>	Project ID, such as ASPxfbD, defualt SSH
	-h or -help		Print this information

	# parameters for HOME/bin/Sff2xxx.pl, 
	# just list them in -sff2xxx as -sff2xxx=\"-m -t -q\";
	#	-s or -seq      Output just the sequences, default yes
	#	-q or -qual     Output just the quality scores, default default yes
	#	-f or -flow     Output just the flowgrams, default no
	#	-t or -tab      Output the seq/qual/flow as tab-delimited lines, default no
	#	-n or -notrim   Output the untrimmed sequence or quality scores, default no
	#	-m or -mft      Output the manifest text, default no
	#	-p or -plain	Output the plain text, default no
	#	-a or -all	Output the all above files, default no

	# parameters for HOME/bin/remove_CS.pl,
	# just list them in -remove_CS as -remove_CS=\"-cat yes -remove yes\"
	#	-CS1 <str>                      Consensus sequence 1
	#	-CS2 <str>                      Consensus sequence 2
	#	    				The default CS pairs are 'ACACTGACGACATGGTTCTACA' (CS1)
	#							       & 'TACGGTAGCAGAGACTTGGTCT' (CS2)
	#	-max_mismatch <int>             default 2 nucletides
	#	-max_gap <int>                  default 2 nucletides
	#	-cat <yes|no>                   cat the *.A/B.*.* to *.deCS, default no
	#	-remove <yes|no>                remove the sequnce without CS, default no

	# parameters for HOME/bin/Fa2Fq.pl,
	# just list them in -fa2fq as -fa2fq=\"-manner Sanger\"
	#	-manner <Sanger|Solexa> Choose the transcoding manner between QA and ASCII, default Solexa
	#				if you choose Sanger, ASCII = QA + 33;
	#				if you choose Solexa, ASCII = QA + 64;


	# parameters for HOME/bin/FAS.pl,
	# just list them in -fas as -fas=\"-Q 30 -B 30\"
	#	-Q <int>  Base with less value than this is defined as low quality bases, default 20
	#	-B <int>  filter reads with >X percent base are low quality bases,set a cutoff, default 40
	#	-w <int>  filter reads with >X percent base are Ns,set a cutoff, default 2
	#	-s <int>  filter reads with length shorter than this threshold value, default 30 nucleotides.

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

";
##====================##
##    get  options    ##
##====================##
my ( $sff2xxx, $remove_CS, $fas, $fa2fq, $lst, $Project, $help );

GetOptions(
    "sff2xxx:s"   => \$sff2xxx,
    "remove_CS:s" => \$remove_CS,
    "fas:s"       => \$fas,
    "fa2fq:s"     => \$fa2fq,
    "lst:s"       => \$lst,
    "project_id:s"=> \$Project,
    "help|?|mam"  => \$help,
);
##====================##
## default parameters ##
##====================##
die print $usage if $help or ( defined $lst && ( @ARGV % 2 ) != 0 ) or ( !defined $lst && @ARGV != 2 );
# <sff_file> and <sample_name> not defined and -lst not defined too

$sff2xxx   ||= "-s -q";
$remove_CS ||= "-cat yes -remove yes";
$fas       ||= "-Q 20 -B 40 -w 2 -s 0";
$fa2fq     ||= "-manner Sanger";
$Project   ||= 'SSH';
##====================##
## preparing working  ##
##====================##
die unless MakeDir($Project, 1);
my $data = "$Project/DATA";
die unless MakeDir($data, 0);
open LOG, ">$Project/SSH.log";
my %files;

if ( @ARGV == 2 ) {
    my ( $sample_name, $sff_file ) = @ARGV;
    unless ( -f $sff_file ) {
        Message("There's no $sff_file, Please check!", \*LOG);
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
            Message("There's no $path file, Please check your sff files lst.", \*LOG);
            die;
        }
        if ( exists $files{$sample}{$path} ) {
            Message("You have two same lines in your sff files lst:\n$sample-->$path", \*LOG);
	    die;
        }
        $files{$sample}{$path} = 1;
    }
    close LST;
}

### get change sff file to fasta and quality files and then filter
foreach my $sample ( sort keys %files ) {
    if ( scalar keys %{ $files{$sample} } > 1 ) {
        Message("You have same Sample name <$sample> in your sff files list, and the path is different.They are:", \*LOG);
	Message(join( "\n", keys %{ $files{$sample} } ), \*LOG);
        Message("If they belong to one sample, please IGNORE this paragraph,", \*LOG);
        Message("or please change their sample name in your sff files list to be different.", \*LOG);
    }

    foreach my $path ( sort keys %{ $files{$sample} } ) {
        ## change sff file to FA and QUAL files
	die Message("Canot execute the $RealBin/Sff2xxx.pl: $!", \*LOG) 
	if system("perl $RealBin/Sff2xxx.pl $sample $path -dir $data $sff2xxx 2>>$Project/SSH.log");
    }

    ## remove CS to generate raw data
    die Message("Canot execute the $RealBin/remove_CS.pl: $!", \*LOG) 
    if system("perl $RealBin/remove_CS.pl -out_dir $data $remove_CS $data/raw.$sample.fa 2>>$Project/SSH.log");
    ## fetch raw quality files
    die Message("Canot execute the $RealBin/fetch_quality_file_from_gff.pl: $!", \*LOG) 
    if system("perl $RealBin/fetch_quality_file_from_gff.pl $data/raw.$sample.fa.gff $data/raw.$sample.qual $data/raw.$sample.qual.deCS 2>>$Project/SSH.log");

    ## filter FA and quality files
    die Message("Canot execute the $RealBin/FAS.pl: $!", \*LOG)
    if system("perl $RealBin/FAS.pl $fas $data/raw.$sample.qual.deCS $data/raw.$sample.fa.deCS >$data/filter.$sample.fa 2>$data/filter.$sample.id");

    die Message("Canot execute the $RealBin/../opt/fishInWinter.pl: $!", \*LOG) 
    if system("perl $RealBin/../opt/fishInWinter.pl --bformat table --fformat fasta $data/filter.$sample.id $data/raw.$sample.qual.deCS >$data/filter.$sample.qual");

    ## change FA and quality files to FQ file
    die Message("Canot execute the $RealBin/Fa2Fq.pl: $!", \*LOG) 
    if system("perl $RealBin/Fa2Fq.pl -fa2fq $fa2fq -quality $data/raw.$sample.qual.deCS -out_pre $data/raw.$sample.fa.deCS $data/raw.$sample.fa.deCS  2>>$Project/SSH.log");
    die Message("Canot execute the $RealBin/Fa2Fq.pl: $!", \*LOG) 
    if system("perl $RealBin/Fa2Fq.pl -fa2fq $fa2fq -quality $data/filter.$sample.qual -out_pre $data/filter.$sample $data/filter.$sample.fa  2>>$Project/SSH.log");
    
    # get reads length
    die Message("Canot execute the $RealBin/../opt/fastaDeal.pl: $!", \*LOG) 
    if system("perl $RealBin/../opt/fastaDeal.pl -attribute id:len $data/raw.$sample.fa.deCS > $data/raw.$sample.fa.deCS.len");

    die Message("Canot execute the $RealBin/../opt/fastaDeal.pl: $!", \*LOG) 
    if system("perl $RealBin/../opt/fastaDeal.pl -attribute id:len $data/filter.$sample.fa > $data/filter.$sample.len");


    # get all reads length
    die Message("Canot execute the awk: $!", \*LOG) 
    if system("awk '{print \"$sample\t\"\$1\"\t\"\$2}' $data/raw.$sample.fa.deCS.len >> $Project/raw.len");

    die Message("Canot execute the awk: $!", \*LOG) 
    if system("awk '{print \"$sample\t\"\$1\"\t\"\$2}' $data/filter.$sample.len >> $Project/filter.len");
}

# get lenth distribution and draw
die Message("$!", \*LOG) 
if system("perl $RealBin/length_distribution.pl -len_col 3 $Project/raw.len > $Project/raw.length_distribution.txt");

die Message("$!", \*LOG) 
if system("perl $RealBin/length_distribution.pl -len_col 3 $Project/filter.len > $Project/filter.length_distribution.txt");

die Message("$!", \*LOG) 
if system("perl $RealBin/../opt/line_diagram.pl $Project/raw.length_distribution.txt -x_title \"Length(bp)\" -y_title \"Contents\" > $Project/raw.length_distribution.svg");

die Message("$!", \*LOG) 
if system("perl $RealBin/../opt/line_diagram.pl $Project/filter.length_distribution.txt -x_title \"Length(bp)\" -y_title \"Contents\" > $Project/filter.length_distribution.svg");

# get quality distribution and draw
die Message("$!", \*LOG) if system("cat $data/raw.*.qual.deCS > $Project/raw.qual");
die Message("$!", \*LOG) if system("cat $data/filter.*.qual > $Project/filter.qual");
die Message("$!", \*LOG) if system("perl $RealBin/quality_distribution.pl $Project/raw.qual > $Project/raw.quality_distribution.txt");
die Message("$!", \*LOG) 
if system("perl $RealBin/quality_distribution.pl $Project/filter.qual > $Project/filter.quality_distribution.txt");
die Message("$!", \*LOG) if system("rm $Project/raw.qual $Project/filter.qual");
die Message("$!", \*LOG) 
if system("perl $RealBin/../opt/line_diagram.pl $Project/raw.quality_distribution.txt -x_title \"Length(bp)\" -y_title \"Contents\" > $Project/raw.quality_distribution.svg");
die Message("$!", \*LOG) if system("perl $RealBin/../opt/line_diagram.pl $Project/filter.quality_distribution.txt -x_title \"Length(bp)\" -y_title \"Contents\" > $Project/filter.quality_distribution.svg");

# get statistic result
die Message("$!", \*LOG) if system("perl $RealBin/statistic.pl $Project/raw.len > $Project/raw.statistic.xls");
die Message("$!", \*LOG) if system("perl $RealBin/statistic.pl $Project/filter.len > $Project/filter.statistic.xls");

# change file name
die Message("$!", \*LOG) if system("rm -f $data/raw.*.fa $data/raw.*.qual $data/*.len $data/*.id");

close LOG;

##==================##
##    Subroutines   ##
##==================##
sub MakeDir {
    ## $clean = 0/1, 1 means delete files if exists the DIR.
    ## it returns 0/1, 1 means Success while 0 fail.
    my ( $new_dir, $clean ) = @_;
    if ( $new_dir eq '.' ) {
        Message( "The dir is current working dir, no need to creat!", \*STDERR);
        return 1;
    }
    if ( -d $new_dir ) {
        Message( "Allready exists $new_dir", \*STDERR);
        if ($clean) {
            if ( system("rm -rf $new_dir/*") ) {
                Message( "Failed to delete file(s) in $new_dir", \*STDERR);
                return 0;
            }
            else {
                Message( "Successfully delete file(s) in $new_dir", \*STDERR);
                return 1;
            }
        }
    }
    else {
        if ( mkdir $new_dir ) {
            Message( "Successfully creat $new_dir", \*STDERR);
            return 1;
        }
        else {
            Message( "Failed to creat $new_dir", \*STDERR);
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
    }else{
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

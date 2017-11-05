#!/usr/bin/perl -w
use strict;
use FindBin qw($RealBin);
use File::Basename qw(basename);
use Getopt::Long;
use lib "$RealBin/../lib";
use DPalign;
use ForkManager;
use Data::Dumper;

my $usage = "
Usage:
	perl $0 [Options] <seqence_files> 2>remove_CS.log

Author:
	DENG Cao
	2012-08-30
	brentcaodeng\@hotmail.com
	dengcao\@bgitechsolutions.com

Options:
	-CS1 <str>			Consensus sequence 1
	-CS2 <str>			Consensus sequence 2
					The default CS pairs are 'ACACTGACGACATGGTTCTACA' (CS1) 
							       & 'TACGGTAGCAGAGACTTGGTCT' (CS2)
	-max_mismatch <int>		default 2 nucletides
	-max_gap <int>			default 2 nucletides

	-type <fasta|fastaq|sff>	the format of input sequence file, default fasta
	-cat <yes|no>			cat the *.A/B.*.* to *.deCS, default no
	-remove	<yes|no>		remove the sequnce without CS, default no
	-parallel <int>			parallel running jobs number, default 500
	-out_pre <str>			default basename(\$ARGV[0])
	-out_dir <str>			default .
	-help|man|?			print this help information
Note:
	1. <seqence_files>: gz readable
	2. '-CS1' & '-CS2': both or none
	3. The output file will be the same format as input. Therefore, if the sequence 
	   file is fasta foramt, you have to use the *.gff and quality file to generate
	   fastaq file if you need this format.

Desciption:
	Please note there might be a small portion of mismatch due to sequencing error
	(just as with the barcodes), so please take that into account. So the sequences
	are at the beginning and end of the reads (maybe not at the end of all reads if
	the molecule was not sequenced in total). 
	
	There are two possible situations:
		1. You should find at one end:
			Consensus sequence 1: ACACTGACGACATGGTTCTACA 
		   And at the other end:
			Consensus sequence 2: TACGGTAGCAGAGACTTGGTCT
			       (Reverse comp: AGACCAAGTCTCTGCTACCGTA)

		2. You should find at one end:
			Consensus sequence 2: TACGGTAGCAGAGACTTGGTCT
		   And at the other end:
			Consensus sequence 1: ACACTGACGACATGGTTCTACA
			       (Reverse comp: TGTAGAACCATGTCGTCAGTGT)
Example:
	perl remove_CS.pl example.fasta
	perl remove_CS.pl -type fastaq example.fq
	perl remove_CS.pl -CS1 AGCTGCTACGATCGAT -CS2 GCTAGCAGTTCAT example.fasta.gz
";

##====================##
##    get  options    ##
##====================##
my (
    $CS1,    $CS2, $max_mismatch, $max_gap, $type,
    $remove, $cat, $out_pre,      $out_dir, $help,
    $parallel
);
GetOptions(
    "CS1:s"          => \$CS1,
    "CS2:s"          => \$CS2,
    "max_mismatch:i" => \$max_mismatch,
    "max_gap:i"      => \$max_gap,
    "type:s"         => \$type,
    "remove:s"       => \$remove,
    "cat:s"          => \$cat,
    "out_pre:s"      => \$out_pre,
    "out_dir:s"      => \$out_dir,
    "help|man|?"     => \$help,
    "parallel:i"     => \$parallel
);

##====================##
## default parameters ##
##====================##
die print $usage if @ARGV != 1;

if ( !defined $CS1 && !defined $CS2 ) {
    $CS1 = 'ACACTGACGACATGGTTCTACA';
    $CS2 = 'TACGGTAGCAGAGACTTGGTCT';
}
elsif ( defined $CS1 || defined $CS2 ) {
    Message("If you defined one CS, you must define another one!\n********\n");
    die print $usage;
}
$type         ||= 'fasta';
$remove       ||= 'no';
$cat          ||= 'no';
$out_pre      ||= basename( $ARGV[0] );
$out_dir      ||= '.';
$max_mismatch ||= 1;
$max_gap      ||= 1;
$parallel     ||= 500;

##====================##
## preparing working  ##
##====================##
die unless MakeDir( $out_dir, 0 );
my ($seq_file) = @ARGV;
if ( $out_pre =~ /(fa$|fasta$|fa\.gz$|fasta\.gz$)/i && $type ne 'fasta' ) {
    Message( "Your file type is '$type', which is not consistent with the file name: '$1'" );
    die;
}
if ( $out_pre =~ /(fastaq$|fq$|fq\.gz$|fastaq\.gz$)/i && $type ne 'fastaq' ) {
    Message( "Your file type is '$type', which is not consistent with the file name: '$1'" );
    die;
}
if ( $out_pre =~ /(sff$|sff\.gz$)/i && $type ne 'sff' ) {
    Message( "Your file type is '$type', which is not consistent with the file name: '$1'" );
    die;
}

# reverse sequence and caculate the Max & Min, rc means reverse complement
my $RC_CS1 = Reverse($CS1);
die unless $RC_CS1;
my $RC_CS2 = Reverse($CS2);
die unless $RC_CS2;

my $PCS1 = substr($CS1,-11);
my $PCS2 = substr($CS2,-11);

my $RC_PCS1 = Reverse($PCS1);
die unless $RC_PCS1;
my $RC_PCS2 = Reverse($PCS2);
die unless $RC_PCS2;

my $full_MinMaxScore = MinMaxScore($CS1,$CS2,$max_gap,$max_mismatch);
my $part_MinMaxScore = MinMaxScore($PCS1,$PCS2,1,1);
## open filehandles: SEQ
my $sffinfo = "$RealBin/../opt/sffinfo";
if ( $type =~ /sff/i ) {
    if ( $seq_file =~ /\.sff$/i ) {
        open SEQ, "$sffinfo $seq_file |" || die Message("Failed to open sequence file --> $seq_file");
    }
    else {
        Message( "Your file type is 'sff', which is not consistent with the file name:$seq_file" );
        die;
    }
}
else {
    if ( $out_pre =~ /\.gz$/ ) {
        open SEQ, "gunzip -dc $seq_file |" || die Message("Failed to open $seq_file");
    }
    else {
        open SEQ, $seq_file || die Message("Failed to open $seq_file");
    }
}
## A terminal output files: AINTACT & ACS1 & ACS2
my $handle;
open $handle->{AINTACT}, ">$out_dir/$out_pre.A.intact.deCS" || die Message("Failed to open filehandle AINTACT");
open $handle->{ACS1}, ">$out_dir/$out_pre.A.fragment.deCS1" || die Message("Failed to open filehandle ACS1");
open $handle->{ACS2}, ">$out_dir/$out_pre.A.fragment.deRcCS2" || die Message("Failed to open filehandle ACS2");

## B terminal output files: BINTACT & BCS1 & BCS2
open $handle->{BINTACT}, ">$out_dir/$out_pre.B.intact.deCS" || die Message("Failed to open filehandle BINTACT");
open $handle->{BCS2}, ">$out_dir/$out_pre.B.fragment.deCS2" || die Message("Failed to open filehandle BCS2");
open $handle->{BCS1}, ">$out_dir/$out_pre.B.fragment.deRcCS1" || die Message("Failed to open filehandle BCS1");

open $handle->{NON}, ">$out_dir/$out_pre.NON.fragment.noCS" || die Message("Failed to open filehandle NON");

## Cannot distiguish A or B terminal: NON
unless ( $remove =~ /y/i ) {
    open $handle->{CSONLY}, ">$out_dir/$out_pre.CsOnly" || die Message("Failed to open filehandle CSONLY");
    open $handle->{CSWRONG}, ">$out_dir/$out_pre.CsWrong" || die Message("Failed to open filehandle CSWRONG");
}
## GFF file for generating qualtity files for fasta format sequence file: GFF
open  $handle->{GFF}, ">$out_dir/$out_pre.gff" || die Message("Failed to open filehandle GFF");
	     
##=====================##
## process fork result ##
##=====================##
my $pm = new ForkManager($parallel);
$pm->run_on_finish(
  sub{
      my ($no_result,$data) = ( $_[1], $_[-1] );
      unless ($no_result){
          foreach (keys %{$data}){
	      my $fh = $handle->{$_};
	     print $fh $data->{$_};
	  }
      }
  }
);
##====================##
## process  sff  file ##
##====================##
if ( $type =~ /^sff$/i ) {
    $/ = "\n>";
    ## waiting for source code
    $/ = "\n";
}

##====================##
## process fasta file ##
##====================##
elsif ( $type =~ /^fasta$/i ) {
    $/ = "\n>";
    while (<SEQ>) {
        chomp;
        s/^>//;
        next if /^\s+$/;
        s/^\s+|\s+$//;
	$pm->start and next; # do the fork
	my $no_result = 1;
	my $data;
        my ( $head, $seq ) = split /\n/, $_, 2;
        my ($id) = $head =~ /^(\S+)/;
        $seq =~ s/\s+//g;
        my $trimCS = TrimCS( $seq, $CS1, $CS2 );
        my $full_length = length $seq;

        ## print CS_WRONG adn CS_ONLY
	unless ( $remove =~ /y/i ) {
            if (exists $trimCS->{cs_only}){
	        $data->{CSONLY} .= ">$id\n$seq\n";
		$no_result = 0;
	    }
	    if (exists $trimCS->{cs_wrong}){
	        $data->{CSWRONG} .= ">$id\n$seq\n";
		$no_result = 0;
	    }
        }
        
	## print NON
	if ( exists $trimCS->{non} ) {
            $data->{NON} .= ">$id\n$seq\n";
            $data->{GFF} .= "$id\t$full_length\t1\t$full_length\t-\t-\t-\t-\n";
	    $no_result = 0;
        }
	## print A
        elsif ( exists $trimCS->{A} ) {
            foreach my $file_type ( keys %{ $trimCS->{A} } ) {
	        next if (!defined $trimCS->{A}->{$file_type}->{new_seq} || $trimCS->{A}->{$file_type}->{new_seq} eq '');
                $data->{GFF} .= "$id\t$full_length\t$trimCS->{A}->{$file_type}->{start}\t$trimCS->{A}->{$file_type}->{end}\t";
                $data->{GFF} .= "$trimCS->{A}->{$file_type}->{cs1}\t$trimCS->{A}->{$file_type}->{rc_cs2}\n";
                $data->{AINTACT} .= ">$id\n$trimCS->{A}->{$file_type}->{new_seq}\n" if $file_type eq 'intact';
                $data->{ACS1} .= ">$id\n$trimCS->{A}->{$file_type}->{new_seq}\n" if $file_type eq 'cs1';
                $data->{ACS2} .= ">$id\n$trimCS->{A}->{$file_type}->{new_seq}\n" if $file_type eq 'rc_cs2';
		$no_result = 0;
            }
        }
	## print B
        elsif ( exists $trimCS->{B} ) {
            foreach my $file_type ( keys %{ $trimCS->{B} } ) {
	        next if (!defined $trimCS->{B}->{$file_type}->{new_seq} || $trimCS->{B}->{$file_type}->{new_seq} eq '');
                $data->{GFF} .= "$id\t$full_length\t$trimCS->{B}->{$file_type}->{start}\t$trimCS->{B}->{$file_type}->{end}\t";
                $data->{GFF} .= "$trimCS->{B}->{$file_type}->{cs2}\t$trimCS->{B}->{$file_type}->{rc_cs1}\n";
                $data->{BINTACT} .= ">$id\n$trimCS->{B}->{$file_type}->{new_seq}\n" if $file_type eq 'intact';
                $data->{BCS1} .= ">$id\n$trimCS->{B}->{$file_type}->{new_seq}\n" if $file_type eq 'rc_cs1';
                $data->{BCS2} .= ">$id\n$trimCS->{B}->{$file_type}->{new_seq}\n" if $file_type eq 'cs2';
		$no_result = 0;
            }
        }
	$pm->finish($no_result,$data); # do the exit in the child process
    }
    $pm->wait_all_children;
    $/ = "\n";
}

##=====================##
## process fastaq file ##
##=====================##
elsif ( $type =~ /^fastaq$/i ) {
    $/ = "\n@";
    while (<SEQ>) {
        chomp;
        s/^@//;
        next if /^\s+$/;
        s/^\s+|\s+$//;
	$pm->start and next; # do the fork
	my $no_result = 1;
	my $data;
        my ( $head, $seq, $dull, $qual ) = split /\n/;
        my ($id) = $head =~ /^(\S+)/;
        my $trimCS = TrimCS( $seq, $CS1, $CS2 );
        my $full_length = length $seq;

        ## print CS_OLY and CS_WRONG
	unless ( $remove =~ /y/i ) {
            if ( exists $trimCS->{cs_only} ) {
	        $data->{CSONLY} .= "\@$id\n$seq\n$dull\n$qual\n";
		$no_result = 0;
	    }
	    if ( exists $trimCS->{cs_wrong} ) {
	        $data->{CSWRONG} .= "\@$id\n$seq\n$dull\n$qual\n";
		$no_result = 0;
	    }
	}

	## print NON
	if ( exists $trimCS->{non} ) {
	    unless ( $remove =~ /y/i ) {
	        $data->{NON} .= "\@$id\n$seq\n$dull\n$qual\n";
	        $data->{GFF} .= "$id\t$full_length\t1\t$full_length\t-\t-\t-\t-\n";
		$no_result = 0;
	    }
	}
	## print A
        elsif ( exists $trimCS->{A} ) {
            foreach my $file_type ( keys %{ $trimCS->{A} } ) {
	        next if (!defined $trimCS->{A}->{$file_type}->{new_seq} || $trimCS->{A}->{$file_type}->{new_seq} eq '');
                my $start = $trimCS->{A}->{$file_type}->{start} - 1;
                my $length = $trimCS->{A}->{$file_type}->{end} - $trimCS->{A}->{$file_type}->{start} + 1;
                my $new_qual = substr( $qual, $start, $length );
                $data->{GFF} .= "$id\t$full_length\t$trimCS->{A}->{$file_type}->{start}\t$trimCS->{A}->{$file_type}->{end}\t";
                $data->{GFF} .= "$trimCS->{A}->{$file_type}->{cs1}\t$trimCS->{A}->{$file_type}->{rc_cs2}\n";
                $data->{AINTACT} .= "\@$id\n$trimCS->{A}->{$file_type}->{new_seq}\n$dull\n$new_qual\n" if $file_type eq 'intact';
                $data->{ACS1} .= "\@$id\n$trimCS->{A}->{$file_type}->{new_seq}\n$dull\n$new_qual\n" if $file_type eq 'cs1';
                $data->{ACS2} .= "\@$id\n$trimCS->{A}->{$file_type}->{new_seq}\n$dull\n$new_qual\n" if $file_type eq 'rc_cs2';
		$no_result = 0;
            }
        }
	## print B
	elsif ( exists $trimCS->{B} ) {
            foreach my $file_type ( keys %{ $trimCS->{B} } ) {
	        next if (!defined $trimCS->{B}->{$file_type}->{new_seq} || $trimCS->{B}->{$file_type}->{new_seq} eq '');
                my $start = $trimCS->{B}->{$file_type}->{start} - 1;
                my $length = $trimCS->{B}->{$file_type}->{end} - $trimCS->{B}->{$file_type}->{start} + 1;
                my $new_qual = substr( $qual, $start, $length );
                $data->{GFF} .= "$id\t$full_length\t$trimCS->{B}->{$file_type}->{start}\t$trimCS->{B}->{$file_type}->{end}\t";
                $data->{GFF} .= "$trimCS->{B}->{$file_type}->{cs2}\t$trimCS->{B}->{$file_type}->{rc_cs1}\n";
                $data->{BINTACT} .= "\@$id\n$trimCS->{B}->{$file_type}->{new_seq}\n$dull\n$new_qual\n" if $file_type eq 'intact';
                $data->{BCS1} .= "\@$id\n$trimCS->{B}->{$file_type}->{new_seq}\n$dull\n$new_qual\n" if $file_type eq 'rc_cs1';
                $data->{BCS2} .= "\@$id\n$trimCS->{B}->{$file_type}->{new_seq}\n$dull\n$new_qual\n" if $file_type eq 'cs2';
		$no_result = 0;
	    }
        }
	$pm->finish($no_result,$data); # do the exit in the child process
    }
    $pm->wait_all_children;
    $/ = "\n";
}
else {
    Message("Unknown sequence file format!");
}

##====================##
## process files over ##
##====================##
if ($cat =~ /y/i){
    die Message("$!") 
    if system("cat $out_dir/$out_pre.A.* $out_dir/$out_pre.B.* $out_dir/$out_pre.NON.fragment.noCS > $out_dir/$out_pre.deCS");
    die Message("$!") if system("rm -f $out_dir/$out_pre.A.* $out_dir/$out_pre.B.* $out_dir/$out_pre.NON.fragment.noCS");
}

##==================##
##    Subroutines   ##
##==================##
sub TrimCS {
    my ( $seq, $cs1, $cs2 ) = @_;
    
    my $trimCS;
    ## $trimCS->{A}
    ## $trimCS->{B}
    ## $trimCS->{non}
    ## $trimCS->{cs_only}
    ## $trimCS->{cs_wrong}

    my $score;
    my $judge = 'full';
    my $do_REDO = 1;
    my ($rc_cs1, $rc_cs2) = ($RC_CS1, $RC_CS2);
    my $cs1_align    = local_affine( $cs1,    $seq );
    my $cs2_align    = local_affine( $cs2,    $seq );

    my ($rc_cs1_align, $rc_cs2_align);
    
    $score->{cs1} = $cs1_align->{score} if $cs1_align->{score} >= $full_MinMaxScore->{cs1}->{min};
    $score->{cs2} = $cs2_align->{score} if $cs2_align->{score} >= $full_MinMaxScore->{cs2}->{min};
    
    ## Need recheck
    REDO:{
        $rc_cs1_align = local_affine( $rc_cs1, $seq );
        $rc_cs2_align = local_affine( $rc_cs2, $seq );
    
        # score standardization
        if ($judge eq 'full'){
            $score->{rc_cs2} = $rc_cs2_align->{score} if $rc_cs2_align->{score} >= $full_MinMaxScore->{cs2}->{min};
            $score->{rc_cs1} = $rc_cs1_align->{score} if $rc_cs1_align->{score} >= $full_MinMaxScore->{cs1}->{min};
        }elsif($judge eq 'part'){
            $score->{rc_cs2} = $rc_cs2_align->{score} if $rc_cs2_align->{score} >= $part_MinMaxScore->{cs2}->{min};
            $score->{rc_cs1} = $rc_cs1_align->{score} if $rc_cs1_align->{score} >= $part_MinMaxScore->{cs1}->{min};
        }
    }
    ## return NON
    if ( scalar keys %{$score} == 0) {
        if ($do_REDO){
	    $do_REDO--;
            ## reset value to run REDO
	    ($rc_cs1,$rc_cs2) = ($RC_PCS1,$RC_PCS2);
	    $judge = 'part';
	    goto REDO;
	}else{
	    $trimCS->{non} = 1;
	    return $trimCS;
	}
    }

    ## return
    foreach my $cs_type ( sort { $score->{$b} <=> $score->{$a} } keys %{$score} ){
        if ( $cs_type eq 'cs1' ) {
            if ( exists $score->{rc_cs2} ) {
                my $RealCS1   = IsRealCS( $cs1,    $cs1_align,  $seq );
	        my $RealRcCS2 = IsRealCS( $rc_cs2, $rc_cs2_align, $seq );
                $trimCS = IntactWrong($seq, $cs1, $rc_cs2, $RealCS1, $RealRcCS2, "A", "cs1", $cs1_align, $rc_cs2_align);
		return $trimCS;
            }
            else {
                my $RealCS1 = IsRealCS( $cs1, $cs1_align, $seq );
                $trimCS->{A}->{cs1}->{start}  = $cs1_align->{ye} + 1;
                $trimCS->{A}->{cs1}->{end}    = length $seq;
                $trimCS->{A}->{cs1}->{cs1}    = "$cs1\t$RealCS1";
                $trimCS->{A}->{cs1}->{rc_cs2} = "$rc_cs2\t-";
                $trimCS->{A}->{cs1}->{new_seq} = substr( $seq, $cs1_align->{ye} );
		$trimCS = CsOnly($seq, $trimCS,"A","cs1");
		return $trimCS;
            }
        }
        elsif ( $cs_type eq 'rc_cs2' ) {
            if ( exists $score->{cs1} ) {
                my $RealCS1   = IsRealCS( $cs1,    $cs1_align,  $seq );
	        my $RealRcCS2 = IsRealCS( $rc_cs2, $rc_cs2_align, $seq );
                $trimCS = IntactWrong($seq, $cs1, $rc_cs2, $RealCS1, $RealRcCS2, "A", "cs1", $cs1_align, $rc_cs2_align);
		return $trimCS;
            }
            else {
                my $RealRcCS2 = IsRealCS( $rc_cs2, $rc_cs2_align, $seq );
                $trimCS->{A}->{rc_cs2}->{start}  = 1;
                $trimCS->{A}->{rc_cs2}->{end}    = $rc_cs2_align->{ys} - 1;
                $trimCS->{A}->{rc_cs2}->{cs1}    = "$cs1\t-";
                $trimCS->{A}->{rc_cs2}->{rc_cs2} = "$rc_cs2\t$RealRcCS2";
                $trimCS->{A}->{rc_cs2}->{new_seq} = substr( $seq, 0, $rc_cs2_align->{ys} - 1 );
		$trimCS = CsOnly($seq, $trimCS,"A","rc_cs2");
		return $trimCS;
            }
        }
        elsif ( $cs_type eq 'cs2' ) {
            if ( exists $score->{rc_cs1} ) {
                my $RealCS2   = IsRealCS( $cs2,    $cs2_align,  $seq );
	        my $RealRcCS1 = IsRealCS( $rc_cs1, $rc_cs1_align, $seq );
                $trimCS = IntactWrong($seq, $cs2, $rc_cs1, $RealCS2, $RealRcCS1,  "B", "cs2", $cs2_align, $rc_cs1_align);
		return $trimCS;
            }
            else {
                my $RealCS2 = IsRealCS( $cs2, $cs2_align, $seq );
                $trimCS->{B}->{cs2}->{start}  = $cs2_align->{ye} + 1;
                $trimCS->{B}->{cs2}->{end}    = length $seq;
                $trimCS->{B}->{cs2}->{cs2}    = "$cs2\t$RealCS2";
                $trimCS->{B}->{cs2}->{rc_cs1} = "$rc_cs1\t-";
                $trimCS->{B}->{cs2}->{new_seq} = substr( $seq, $cs2_align->{ye} );
		$trimCS = CsOnly($seq, $trimCS,"B","cs2");
		return $trimCS;
            }
        }
        elsif ( $cs_type eq 'rc_cs1' ) {
            if ( exists $score->{cs2} ) {
                my $RealCS2   = IsRealCS( $cs2,    $cs2_align,  $seq );
	        my $RealRcCS1 = IsRealCS( $rc_cs1, $rc_cs1_align, $seq );
                $trimCS = IntactWrong($seq, $cs2, $rc_cs1, $RealCS2, $RealRcCS1,  "B", "cs2", $cs2_align, $rc_cs1_align);
		return $trimCS;
            }
            else {
                my $RealRcCS1 = IsRealCS( $rc_cs1, $rc_cs1_align, $seq );
                $trimCS->{B}->{rc_cs1}->{start}  = 1;
                $trimCS->{B}->{rc_cs1}->{end}    = $rc_cs1_align->{ys} - 1;
                $trimCS->{B}->{rc_cs1}->{cs2}    = "$cs2\t-";
                $trimCS->{B}->{rc_cs1}->{rc_cs1} = "$rc_cs1\t$RealRcCS1";
                $trimCS->{B}->{rc_cs1}->{new_seq} = substr( $seq, 0, $rc_cs1_align->{ys} - 1 );
		$trimCS = CsOnly($seq, $trimCS,"B","rc_cs1");
		return $trimCS;
            }
        }
        last;
    }
    return $trimCS;
}

sub IntactWrong{
    my ( $seq, $cs, $rc_cs, $RealCS, $RealRcCS, $type, $new_class, $first_align, $second_align ) = @_;
    my $trimCS;

    ##                          2----------
    ## 1----------------------
    if ($first_align->{ye} < $second_align->{ys}){
        $trimCS->{$type}->{intact}->{start}  = $first_align->{ye} + 1;
	$trimCS->{$type}->{intact}->{end}    = $second_align->{ys} - 1;
	$trimCS->{$type}->{intact}->{$new_class} = "$cs\t$RealCS";
	if ($new_class eq 'cs1'){
           $trimCS->{$type}->{intact}->{rc_cs2} = "$rc_cs\t$RealRcCS";
	}
	if ($new_class eq 'cs2'){
           $trimCS->{$type}->{intact}->{rc_cs1} = "$rc_cs\t$RealRcCS";
	}
	$trimCS->{$type}->{intact}->{new_seq} = substr( $seq, $first_align->{ye}, $second_align->{ys} - $first_align->{ye} - 1 );
        $trimCS = CsOnly($seq, $trimCS, $type, "intact");
	return $trimCS;
    }

    ##                    2-----------
    ## 1-----------------------
    elsif ($first_align->{ye} >= $second_align->{ys} && $first_align->{ye} <= $second_align->{ye}){
        $trimCS->{cs_only} = 1;
	return $trimCS;
    }
    ##            2---------  |           2-----    | 2------------------
    ## 1--------------------  | 1-----------------  |                     1-------------------
    elsif ($first_align->{ye} >= $second_align->{ye}){
       if ((length $cs ) == (length $rc_cs)){
	   $trimCS->{cs_wrong} = 1;
	   return $trimCS;
       }else{
           $trimCS->{$type}->{$new_class}->{start}  = $first_align->{ye} + 1;
           $trimCS->{$type}->{$new_class}->{end}    = length $seq;
           $trimCS->{$type}->{$new_class}->{$new_class}    = "$cs\t$RealCS";
	   if ($new_class eq 'cs1'){
               $trimCS->{$type}->{$new_class}->{rc_cs2} = "$rc_cs\t-";
	   }
	   if ($new_class eq 'cs2'){
               $trimCS->{$type}->{$new_class}->{rc_cs1} = "$rc_cs\t-";
	   }
           $trimCS->{$type}->{$new_class}->{new_seq} = substr( $seq, $second_align->{ye} );
           $trimCS = CsOnly($seq, $trimCS, $type, $new_class);
	   return $trimCS;
       }
    }
}

sub CsOnly{
    my ( $seq, $trimCS, $type, $class ) = @_;
    if ($trimCS->{$type}->{$class}->{start} - (length $seq) > 0 ){
#        print "$trimCS->{$type}->{$class}->{start},$trimCS->{$type}->{$class}->{end}\n",length $seq,"\n";
#        print "$trimCS->{$type}->{$class}->{cs1},$trimCS->{$type}->{$class}->{rc_cs2}\n" if $type eq 'A';
#        print "$trimCS->{$type}->{$class}->{cs2},$trimCS->{$type}->{$class}->{rc_cs1}\n" if $type eq 'B';
#	print "$seq\n\n\n\n";
        my $trimCS_new;
	$trimCS_new->{cs_only} = 1;
	return $trimCS_new;
    }else{
        return $trimCS;
    }
}

sub MinMaxScore {
    my ( $cs1, $cs2, $gap, $mismatch ) = @_;
    my $ref; # score scheme used by DPalign Module: match = 2, mismatch = 1; gap = -3; gap_extension = -1
    $ref->{cs1}->{min} = ( ( length $cs1 ) - $gap - $mismatch ) * 2 + $mismatch * 1 - 3 * $gap;
    $ref->{cs2}->{min} = ( ( length $cs2 ) - $gap - $mismatch ) * 2 + $mismatch * 1 - 3 * $gap;
    $ref->{cs1}->{max} = ( length $cs1 ) * 2;
    $ref->{cs2}->{max} = ( length $cs2 ) * 2;
    return $ref;
}

sub IsRealCS {
    my ( $cs, $align, $seq ) = @_;
    my $real = $align->{y};
    unless ( ( length $cs ) == ( length $align->{x} ) ) {
        my $left  = $align->{xs} - 1;
        my $right = ( length $cs ) - $align->{xe};
        if ( ( $align->{ys} - $left ) < 1 ) {
            $align->{ys} = 1;
        }
        else {
            $align->{ys} = $align->{ys} - $left;
        }
        if ( ( $align->{ye} + $right ) > ( length $seq ) ) {
            $align->{ye} = length $seq;
        }
        else {
            $align->{ye} = $align->{ye} + $right;
        }
        $real = substr( $seq, $align->{ys} - 1, $align->{ye} - $align->{ys} + 1 );
    }
    $real =~ s/-//g;
    return $real;
}

sub Reverse {
    my ($sense) = @_;
    my $antisense = '';
    if ( $sense =~ /[^ATGCatgc]+/g ) {
        Message( "Your Consensus Sequence seems wrong --> $sense, please check!" );
        return $antisense;
    }
    $antisense = reverse $sense;
    $antisense =~ tr/ATGCatgc/TACGtacg/;
    return $antisense;
}

sub MakeDir {
    ## $clean = 0/1, 1 means delete files if exists the DIR.
    ## it returns 0/1, 1 means Success while 0 fail.
    my ( $dir, $clean ) = @_;
    if ( $dir eq '.' ) {
        Message("The dir is current working dir, no need to creat!");
        return 1;
    }
    if ( -d $dir ) {
        Message("Allready exists $dir");
        if ($clean) {
            if ( system("rm -rf $dir/*") ) {
                Message("Failed to delete file(s) in $dir");
                return 0;
            }
            else {
                Message("Successfully delete file(s) in $dir");
                return 1;
            }
        }else{
	    return 1;
	}
    }
    else {
        if ( mkdir $dir ) {
            Message("Successfully creat $dir");
            return 1;
        }
        else {
            Message("Failed to creat $dir");
            return 0;
        }
    }
}

sub Message {
    my ($infor) = @_;
    my @new = split /\n+/, $infor;
    foreach (@new) {
        print STDERR Date() . "\t\t$_\n";
    }
    print STDERR "\n";
}

sub Date {
    my $time = `date`;
    $time =~ s/^\s+|\s+$//;
    return '[' . $time . ']';
}

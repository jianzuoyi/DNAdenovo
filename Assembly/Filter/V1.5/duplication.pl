#!/usr/bin/perl
use strict;

die "Usage:$0 <read1.fq> <read2.fq> <clean1.fq> <clean2.fq> <dup.stat>\n" if @ARGV<5;

my $fq1 = shift;
my $fq2 = shift;
my $outfile1 = shift;
my $outfile2 = shift;
my $statfile = shift;

my %Seq;
my $Duplicate = 0;

open OUT,">$statfile" or die "$!";
print OUT "Start time:".time."\n";

open IN1,$fq1 or die "$!";
open IN2,$fq2 or die "$!";
open OUT1,">$outfile1" or die "$!";
open OUT2,">$outfile2" or die "$!";
my $num = 0;
my $clean = 0;
while(my $id1 = <IN1> ){
	$num++;
	my $seq1 = <IN1>;
	chomp $seq1;
	my $qid1 = <IN1>;
	my $q1 = <IN1>;
	my $id2 = <IN2>;
	my $seq2 = <IN2>;
	chomp $seq2;
	my $qid2 = <IN2>;
	my $q2 = <IN2>;
	
	if ( not defined $Seq{"$seq1$seq2"} ){
			print OUT1 "$id1$seq1\n$qid1$q1";
			print OUT2 "$id2$seq2\n$qid2$q2";
			$clean ++;;
			$Seq{"$seq1$seq2"} ++;
			
	}else{
		$Seq{"$seq1$seq2"}++;
	}
}
close IN1;
close IN2;
close OUT1;
close OUT2;

foreach ( sort values %Seq ){
	$Duplicate += $_ if $_ > 1;
}

print OUT "Total_reads:".($num)."\n"."Duplicate_reads:".($Duplicate)."\n";
print OUT "Clean_reads:".($clean)."\n";
print OUT "END time:".time."\n";
close OUT;


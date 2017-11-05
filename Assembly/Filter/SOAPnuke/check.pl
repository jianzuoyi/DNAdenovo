#!/usr/bin/perl
=head1 Version
Author: tangqiqin@genomics.org.cn
Date: 2012.7.18
=head1 Example
        perl check.pl
        -o:	old Fq1.gz file
	-n:	new Fq1.gz file
=cut

use strict;
use warnings;
use Getopt::Long;

my ($old,$new);
GetOptions(
        "o:s" =>\$old,
        "n:s"   =>\$new
);

die `pod2text  $0` if( !$old or !$new);

if($new=~/.*\.gz/){
	open N,"gzip -dc $new |";
}else{
	open N,"$new";
}

my %new_id;
my %seq;
my %qul;
my $flag=0;
while(<N>){
	chomp;
	$flag++;
	$new_id{$_}++;
	$seq{$_}=<N>;<N>;$qul{$_}=<N>;
}
close N;
print "Number of New $flag\n";
open O3,">new_repeat.txt";
foreach(keys %new_id){
	if($new_id{$_}>1){
		print O3 $_."\t".$new_id{$_}."\n";
	}
}
close O3;
if($old=~/.*\.gz/){
	open O,"gzip -dc $old |";
}else{
	open O,"$old";
}
open O1,">add_new_file.fq";
open O2,">add_old_file.fq";
$flag=0;
while(<O>){
	$flag++;
	chomp;
	my $seq=<O>;
	my $jia=<O>;
	my $qul=<O>;
	if( ! exists $new_id{$_} ){  #在新个过滤结果中不存在，也就是
		print O2 $_."\n";
		print O2 $seq.$jia.$qul;
	}else{
		$new_id{$_}=0;
	}
}
print "Number of Old $flag\n";

foreach(keys %new_id){
        if($new_id{$_}>0){
                print O1 $_."\n".$seq{$_}."+\n".$qul{$_};
        }
}

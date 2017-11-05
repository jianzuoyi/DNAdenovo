#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

=head1 Description

	It is a program for simulating solexa PE reads,with a series of problems generate by solexa 
	sequencing machine in thought,such as insertsize distribution,error rate,heterozygosis SNP
	and heterozygosis Indel in diploid. User should set the value of insertsize_mean and 
	insertsize_sd ,they are the mean value and standard deviation of the normal distribution
	that used as the model function when simulating insertsize distribution,usually we set the 
	insertsize_sd to be 1/20 of the insertsize_mean.The normal distribution function model we
	used in this program is f(x)=1/��/sqrt(2*pi)/exp((x-��)**2 / (2*��**2)) ,and the insertsize 
	distribution range is limited in (��-5�ң���+5��), which will cover almost all the data.
	To simulate illumina error rates on different cycles,we use the function f(x)=0.00001*x**3 
	as a model function,because the error rate of a number of bases at the end is much larger 
	than other area in a read. User should set the heterozygosis SNP rate and heterozygosis Indel
	rate as also ,if you want to use those function,but remember that heterozygosis SNP rate 
	and heterozygosis Indel rate is only exists in diploid. At last ,you should set another
	several parameters ,read length ,coverage of reads ,input sequence and output prefix,the 
	input sequence must be set ,because there is not default value.

=head1 Contact and Version

	Contact:	lujianliang@genomics.org.cn yuezhen@genomics.org.cn
	Version:	1.0		Data:	2010-2-24

=head1 Option

	-input			<string>	input reference genome sequence
	-read_len		<int>		set read length, read1 and read2 have the same length,default:100
	-coverage		<int>		set the sequencing coverage(sometimes called depth),default:40
	-insertsize_mean	<int>		set the average value of insert size,default:500
	-insertsize_sd		<int>		set the standard deviation of insert sizes, default:25
	-error_rate		<float>		set the average error rate over all cycles,default:0.01
	-heterSNP_rate		<float>		set the heterozygous SNP rate of the diploid genome,default:0
	-heterIndel_rate	<float>		set the heterozygous indel rate of the diploid genome,default:0
	-output			<string>	output file prefix default:solexa
	-help					output help infomation

=head1 Example

	1. perl simulate_solexa_reads.pl -input ref_sequence.fa 
	every parameter use the default one. 
	2. perl simulate_solexa_reads.pl -input ref_sequence.fa -read_len 150 -coverage 20 -o humen_500_100
	just set read length and coverage you needed.
	3. perl simulate_solexa_reads.pl -input ref_sequence.fa -o humen -insertsize_mean 600 -insertsize_sd 30 
	   -error_rate 0.01
	set insertsize distribution and error rate.
	4. perl simulate_solexa_reads.pl -input ref_sequence.fa -o humen -heterSNP_rate 0.001 -heterIndel_rate 0.001
	the genome is diploid and you want to produce heterozygosis SNPs  heterozygosis Indels in reads.

=cut
################################################################################################
my ($inputref,$read_len,$coverage,$insertsize_mean,$insertsize_sd,$error_rate,$heterSNP_rate);
my ($heterIndel_rate,$ref_ploid,$n_ploid,$o,$help);

GetOptions(
	"input:s"=>\$inputref,
	"read_len:i"=>\$read_len,
	"coverage:i"=>\$coverage,
	"insertsize_mean:i"=>\$insertsize_mean,
	"insertsize_sd:i"=>\$insertsize_sd,
	"error_rate:f"=>\$error_rate,
	"heterSNP_rate:f"=>\$heterSNP_rate,
	"heterIndel_rate:f"=>\$heterIndel_rate,
	"output:s"=>\$o,
	"help"=>\$help
);
die `pod2text $0` if ($help);
die `pod2text $0` if (!$inputref);
#############################Set default value#############################################
$read_len=100 if (!defined $read_len);
$coverage=40 if (!defined $coverage);
$insertsize_mean=500 if (!defined $insertsize_mean);
$insertsize_sd=25 if (!defined $insertsize_sd);
$error_rate=0.01 if (!defined $error_rate);
$heterSNP_rate=0 if (!defined $heterSNP_rate);
$heterIndel_rate=0 if (!defined $heterIndel_rate);
$o="solexa" if (!defined $o);

die ("reads length is no bigger than 0\n") if ($read_len<=0);
die ("coverage is no bigger than 0\n") if ($coverage<=0);
die ("the mean insertsize is no bigger than 0\n") if ($insertsize_mean<=0);
die ("the insertsize_sd is no bigger than 0\n") if ($insertsize_sd<=0);
die ("the heterSNP_rate is smaller than 0\n") if ($heterSNP_rate<0);
die ("the heterIndel_rate is smaller than 0\n") if ($heterIndel_rate<0);
#############################Main Program##################################################
my ($O1,$O2);
open (REF,"$inputref") || die ("can't open file $inputref\n");
open ($O1,">$o\_1.fa") || die ("can't open file $o\_1.fa\n");
open ($O2,">$o\_2.fa") || die ("can't open file $o\_2.fa\n");
$/=">";

while (<REF>) {
	chomp;
	next if ($_ eq "");
	my $line=$_;
	my $id=(split (/\s+/,$line))[0];
	my $refseq=(split (/\n/,$line,2))[1];
	$refseq=~s/\n//g;
	$refseq=~s/[BDEFHIJKLMNOPQRSUVWXYZ]+//g;
	my $length=length $refseq;
	my ($read_num_pair,%insertsize,%read_info);
	if ($heterSNP_rate>0 || $heterIndel_rate>0) {
		$read_num_pair=int ($length*$coverage/(2*2*$read_len));
	}else{
		$read_num_pair=int ($length*$coverage/(2*$read_len));
	}
	&insertsize_distribution($insertsize_mean,$insertsize_sd,$read_num_pair,\%insertsize);
	if ($error_rate>0){
		if ($error_rate>0.001){
			&error_distribution($length,$read_len,$error_rate,$read_num_pair,\%read_info);
		}else{
			print "The error rate is smaller than basic error rate 0.001,you'd better set the value either bigger or 0\n";
			exit;
		}
	}
	&getreads($refseq,\%insertsize,\%read_info);
	
	%insertsize=();
	%read_info=();
	if ($heterSNP_rate>0 || $heterIndel_rate>0) {	###heterozygous SNP and heterozygous indel exists is diploid
		my $refseqsnp=&heter_SNP_Indel($refseq,$heterSNP_rate,$heterIndel_rate,$id);
		$length=length $refseqsnp;
		$read_num_pair=int ($length*$coverage/(2*2*$read_len));
		&insertsize_distribution($insertsize_mean,$insertsize_sd,$read_num_pair,\%insertsize);
		if ($error_rate>0){
                	if ($error_rate>0.001){
                        	&error_distribution($length,$read_len,$error_rate,$read_num_pair,\%read_info);
                	}else{
                        	print "The error rate is smaller than basic error rate 0.001\n";
				exit;
                	}
        	}
		&getreads($refseqsnp,\%insertsize,\%read_info);
	}
}
$/="\n";
close REF;
close $O1;
close $O2;
#############################Subroutines###################################################
########simulate the insertsize distribution with the model of normal distribution function
########The insertsize range is limited in (��-5�ң���+5��), which covers almost all the data.
sub insertsize_distribution{		
	(my $mean,my $sd,my $rd_num,my $isize)=@_;
	my ($total,$total1);
	my $pi=3.1415926535;
	for (my $i=$mean-5*$sd;$i<=$mean+5*$sd;$i++) {
		$isize->{$i}=1/$sd/sqrt(2*$pi)/exp(($i-$mean)**2 / (2*$sd**2));
		$total+=$isize->{$i};
	}
	for (my $i=$mean-5*$sd;$i<=$mean+5*$sd;$i++) {
		$isize->{$i}=0.5+$isize->{$i}/$total*$rd_num;
		$isize->{$i}=int $isize->{$i};
		$total1+=$isize->{$i};
	}
	$isize->{$mean}+=$rd_num-$total1 if ($rd_num>$total1);
}

#######Simulate illumina error distribution on different cycles,
#######with the model function f(x)=0.00001*x**3
sub error_distribution{		
	(my $seq_len,my $rd_len,my $er_rate,my $rd_num,my $rd_info)=@_;
	my (%hash,$total);
	my $basic_err=0.001;
	my $total_err=($er_rate-$basic_err)*$rd_len;
	for (my $i=1;$i<=$rd_len;$i++) {
		$hash{$i}=0.00001*$i**3;
		$total+=$hash{$i};
	}
	for (my $circle=1;$circle<=$rd_len;$circle++) {
		$hash{$circle}=$basic_err+$hash{$circle}/$total*$total_err;
		$hash{$circle}=$hash{$circle}*$rd_num*2;
		$hash{$circle}=int $hash{$circle};
	}
	
	for(my $circle=1;$circle<=$rd_len;$circle++) {
		while ($hash{$circle}>0) {
			my $j=int (rand $rd_num);
			push @{$rd_info->{$j}},$circle;
			$hash{$circle}--;
		}
	}
	undef %hash;
}
####Rrealization of error sequencing

sub match{
	(my $seq,my $position,my $r_l)=@_;
	$position=int $position;
	my %get_match=(
                "A"=>["T","G","C"],
                "T"=>["A","G","C"],
                "G"=>["A","T","C"],
                "C"=>["A","T","G"],
                "a"=>["t","g","c"],
                "t"=>["a","g","c"],
                "g"=>["a","t","c"],
                "c"=>["a","t","g"],
        );
	my $s1=substr $seq,0,$position-1;
	my $s2=substr $seq,$position-1,1;
	my $s3=substr $seq,-($r_l-$position),$r_l-$position;
	my $k=int ((rand 100)%3);
	if ($s2 eq "A") {
		$s2=$get_match{A}[$k];
	}elsif($s2 eq "T"){
		$s2=$get_match{T}[$k];
	}elsif($s2 eq "G"){
		$s2=$get_match{G}[$k];
	}elsif($s2 eq "C"){
		$s2=$get_match{C}[$k];
	}elsif($s2 eq "a"){
		$s2=$get_match{a}[$k];
	}elsif($s2 eq "t"){
		$s2=$get_match{t}[$k];
	}elsif($s2 eq "g"){
		$s2=$get_match{g}[$k];
	}else{
		$s2=$get_match{c}[$k];
	}
	my $s;
	if ($position==1) {
		$s=$s2.$s3;
	}elsif($position==$r_l){
		$s=$s1.$s2;
	}else{
		$s=$s1.$s2.$s3;
	}
	return $s;
}

#####Getting insert sequence when heterozygous indel rate is bigger than 0
sub getinsertion{
	my $num=shift @_;
	my $base;
	my @baseA=("A","T","G","C");
	for (my $i=1;$i<=$num;$i++) {
		my $j=int ((rand 100)%4);
		$base.=$baseA[$j];
	}
	return $base;
}
#####Produce heterozygous SNPs and indels in multiploid
sub heter_SNP_Indel{			
	(my $sequence,my $hsnp_rate,my $hindel_rate,my $ID)=@_;
#	print "It's begin to get snp sequence !\n";
	my %get_match=(
                "A"=>["T","G","C"],
                "T"=>["A","G","C"],
                "G"=>["A","T","C"],
                "C"=>["A","T","G"],
                "a"=>["t","g","c"],
                "t"=>["a","g","c"],
                "g"=>["a","t","c"],
                "c"=>["a","t","g"],
        );
	my $seqlen=length $sequence;
	my @list;
	@list=split /\s*/,$sequence;
	$sequence="";
	if ($hsnp_rate>0) {
		my $snp_num=$seqlen*$hsnp_rate;
		$snp_num=int $snp_num;
		for (my $i=1;$i<=$snp_num;$i++) {
#			print "$ID\t";
			my $j=int (rand $seqlen);
#			my $q=$j+1;
#			print "$q\t$list[$j]\t";
			my $k=$j%3;
			$list[$j]= $get_match{$list[$j]}[$k];
#			print "$list[$j]\t$snp_num\t$seqlen\n";
		}
		if (($hsnp_rate>0) && ($hindel_rate<=0)){
			foreach my $g (@list){
				$sequence.=$g;
			}
		}
	}
	
	if ($hindel_rate>0) {
		my %insertion;
		my $indel=$seqlen*$hindel_rate;
		my $insertion=$indel/2;
		my $deletion=$indel/2;
		my @list1=(2,3,6);
		my $p=1;
		foreach my $k (@list1) {
			for (my $i=1;$i<=$insertion/$k;$i++) {
				my $j=int (rand $seqlen);
				my $insert=&getinsertion($p);
				$insertion{$j}=$insert;		
			}
			$p++;
		}
		$p=1;
		foreach my $i (@list1) {
			for (my $j=1;$j<=$deletion/$i;$j++) {
				my $k=int (rand $seqlen);
				redo if ($k+$p>$seqlen);
				for (my $t=$k;$t<$k+$p;$t++){
					$list[$t]="";
				}
			}
			$p++;
		}
		my $cs=0;
		foreach my $ce (sort {$a<=>$b} keys %insertion){
			for (my $i=$cs;$i<=$ce;$i++){
                               $sequence.=$list[$i];

                	}
			$sequence.=$insertion{$ce};
			$cs=$ce+1;
		}
		for (my $i=$cs;$i<$seqlen;$i++){
			$sequence.=$list[$i];

        }
		undef %insertion;
		undef @list;
#		foreach my $i (keys %insertion){
#			delete $insertion{$i};
#		}
#	print "Get deletion has been finished!\n";
	}
#	print "Finish Get snp and indel sequence\n";
	return $sequence;
}
#########Produce solexa reads
sub getreads{
	my ($refseq,$ISize,$Rd_Info)=@_;
	my $length=length $refseq;
	my $read_count=0;
	for (my $i=$insertsize_mean-5*$insertsize_sd;$i<=$insertsize_mean+5*$insertsize_sd;$i++) {
                while ($ISize->{$i}>0) {
                        my $pos=rand $length;
                        $pos=int $pos;
                        my $substr=substr $refseq,$pos,$i;
			redo if (length $substr<$i);
                        my $a=substr $substr,0,$read_len;
                        my $b=substr $substr,-$read_len,$read_len;
			if ($insertsize_mean <=1000){
                        	$b=reverse $b;
                        	$b=~tr/ATGCatgc/TACGtacg/;
			}elsif($insertsize_mean >1000){
				$a=reverse $a;
				$a=~tr/ATGCatgc/TACGtacg/;
			}else{
				print "Insertsize can't be smaller than 0\n";
				exit;
			}
                        $read_count++;
                        if (exists $Rd_Info->{$read_count}) {
				 foreach my $j (@{$Rd_Info->{$read_count}}) {
				 	my $k=int ((rand 10)%2);
				 	if ($k) {
				 		$a=&match($a,$j,$read_len);
					 }else{
				 		$b=&match($b,$j,$read_len);
					 }
				 }
                        }
                        my $d=int ((rand 10)%2);
                        if ($d) {
                                print $O1 ">reads_$read_count\_$read_len\_1\n$a\n";
                                print $O2 ">reads_$read_count\_$read_len\_2\n$b\n";
                        }else{
                                print $O1 ">reads_$read_count\_$read_len\_1\n$b\n";
                                print $O2 ">reads_$read_count\_$read_len\_2\n$a\n";
                        }
                        $ISize->{$i}--;
                }
        }
}


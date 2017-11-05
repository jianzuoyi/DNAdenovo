﻿#!/usr/bin/perl
=head1 Version
  Author: Wang Chongzhi, wangchongzhi@genomics.org.cn
  Version: 0.1,  Date: 2011-9-6
  Version: 0.2,  Date: 2011-9-8
  Version: 0.5,  Date: 2011-9-9
  Version: 1.0,  Date: 2011-9-19 

=head1 Usage
  perl kmer_counter.pl <list_file_of_sequence_data_file(s)>
  --k <num>          Length of mers(default=17, must in [11,31])
  --memory <num>     Memory usage(default=10,in Gbytes,i.e. 32th power of 2 in bytes, 1G at least)
  --prefix <str>     Output prefix(default=count[k]mer_by_[current job id])
  --help             Output help information to screen and exit
  You must confirm that your input file is the list of sequence data file(s) in the absolute path(s) and the sequence data file(s) should be gziped.
Example
  1.work with default options (the most simplest way)
  perl kmer_counter.pl my_sequence_data_file.list
  2.work with user specifed options: (to select --k,--memory,--prefix)
  perl kmer_counter.pl --k 23 --memory 10 --prefix k17mer my_sequence_data_file.list
=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);

##get options from command line into variables and set default values
my ($K, $Memory, $Prefix, $Help, $T, $Time_file, $Hash_size);#$T represents Thread_number.
GetOptions(
        "k:i"=>\$K,
        "memory:f"=>\$Memory,
        "prefix:s"=>\$Prefix,
        "help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);
$K ||= 17;
die "Invaid value of --k :  must be in [11,31]!" if $K<11||$K>31;
$Memory ||= 10;
$Memory=1 if $Memory<1;
$T ||= 1;

#read path from config file
open FI,"$Bin/../config" or die "can't read config file!\n";
my %hash_path;
while(<FI>){
	chomp;
	my($key,$value) = split;
	$hash_path{$key} = $value;
}
close FI;
my $jellyfish_path = $hash_path{"jellyfish_path"};
my $gcc_path = $hash_path{"gcc_path"};

#obtain the list of data file(s) which contain(s) the sequnece in ".gz" format of fa or fq file(s).
my $sources_list_file=pop @ARGV;
open IN,"<$sources_list_file";
my @input_files = <IN>;
$/="";
chomp @input_files;
$/="\n";
close IN;
$Prefix ||= "count".$K."mer_by_".$$;
$Time_file=$Prefix.".time";

&compute_hash_size($K,$Memory);

################################################################################
#####    generate the shellscript for work    ##################################
my $work_file=$Prefix.".sh";
my $log_file=$Prefix.".log";
my $error_file=$Prefix.".error";
open OUT,">$work_file";

#set Environment Variables
my $PATH = $hash_path{"PATH"};
print OUT "export PATH=\"$PATH\"\n";
my $LD_LIBRARY_PATH = $hash_path{"LD_LIBRARY_PATH"};
print OUT "export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH\"\n";
my $LIBRARY_PATH = $hash_path{"LIBRARY_PATH"};
print OUT "export LIBRARY_PATH=\"$LIBRARY_PATH\"\n";

#generate the command line to count kmers
print OUT "cat  @input_files | jellyfish count -m $K -o $Prefix --timing $Time_file -s $Hash_size -t $T -c 8 -C /dev/fd/0 1>$log_file 2>$error_file\n";
#generate the command line to merge the binary files of hash table generated by "jellyfish count"
my $premerged=$Prefix."_*";
my $merged_file=$Prefix.".jf";
print OUT "jellyfish merge -v -o $merged_file $premerged 1>>$log_file 2>>$error_file\n";
#generate the command line to dump the the binary files of hash table,usually useless.
my $dump_file=$Prefix.".dump";
print OUT "jellyfish dump -c -t -o $dump_file $merged_file 1>>$log_file 2>>$error_file\n";
#generate the command line to obtain statistics of kmers
my $stats_file=$Prefix.".stats";
print OUT "jellyfish stats -o $stats_file $merged_file 2>>$error_file\n";
#generate the command line to histogram kmers
#my $high_frequency=511;
my $histo_file=$Prefix.".histo";
print OUT "jellyfish histo -t $T $merged_file | sed 's/ /\t/g' >$histo_file\n";
#print OUT "jellyfish histo -h $high_frequency -t $Threads_number $merged_file | sed 's/ /\t/g' >$histo_file\n";
print OUT "perl /home/wangchongzhi/bin/kmeranalysis.pl -k $K $Prefix";
close OUT;

#####      finish the shellscript for work    ##################################
################################################################################

#   print the informations suggested for qsub this kmers-counting job
$Memory=$Memory."G";
print "That's OK! The shell script is finished in $work_file!\n\n";
print "You can qsub your job as following with some appropriate modifications:\n";
print "\t qsub -cwd -l vf=$Memory -l p=$T -q bc_mem.q -P bc_mem $work_file\n";
print "Here,\"vf=$Memory\" and \"p=$T\" is arranged previously by yourself;\"$work_file\" is generated just now.\n";
print "\"bc_mem.q\" can be replaced by the queue you want and \"bc_mem\" by the name of the Project your job belongs to.\n";



##the subroutine which convert memory to hash size
sub compute_hash_size{
	my ($k,$memory,$t) = @_;
	my $m = $memory * 1073741824;
	$Hash_size=4**$k;
	for (my $l=2*$k;$l>16;$l--){
		#2k-l+r+1, default r=11,log_2_1953.
		my $h=((2*$k-$l+12)/8+0.5)*$Hash_size;#h represent hash_and_counter_usage_modified, 0.5 is the coefficient.
		if($h<$m){
			for ($t=2;$t<=32;$t++){
				last if $h>$m-&basic_usage($t);
			}
			$T=$t-1;
			$Memory=sprintf "%.1f",($h+&basic_usage($T))/1073741824;
			return;
		}
		$Hash_size/=2;
        }
}

sub basic_usage{
	my ($t)=@_;
        #The default value of reprobe is 1953 in jellyfish-1.1.2, which is different from that its help information shows.
        my $max_reprobe ||= 1953;#obtained by adding -v option when using the subroutine MERGE of Jellyfish.
        my $out_buffer_size ||= 20000000;#obtained by "jellyfish count --full-help".
        my $buffer_size ||= 8192;#obtained by "jellyfish count --full-help".
        my $buffers ||= 24;#obtained from much uninteresting test and 'top'.

        my $heap_size = $max_reprobe*($max_reprobe+1)/2;#from the article.
        my $total_heap_size = $heap_size*$t;
        my $total_buffer_size = $buffer_size*$buffers*$t;
        my $total_out_buffer_size = $out_buffer_size*$t;
	return $total_heap_size+$total_buffer_size+$total_out_buffer_size+27000000*$t+10000000;#The size of Jellyfish's code is 5.4M. 
}


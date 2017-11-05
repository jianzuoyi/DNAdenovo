#!user/bin/perl

## dorp the N in the fasta file 

$file=shift; ## file in fasta
$n_num=0;
$num=0;
$gc=0;
$write="genome.fa";
open WR,">$write";
print WR ">test\n";
open IN,$file;
while(<IN>){
	if(/>/){
#		print WR $_;
		next;
	}else{
		chomp;
		@t=split //;
		$print =0;
		for($i=0;$i<@t;$i=$i+1){
			$num=$num+1;
			if($t[$i] eq "a" || $t[$i] eq "A" || $t[$i] eq "g" || $t[$i] eq "G" || $t[$i] eq "t" || $t[$i] eq "T" || $t[$i] eq "c" || $t[$i] eq "C"){
				print WR $t[$i];$print=1;
				if($t[$i] eq "G" ||$t[$i] eq "C"){$gc=$gc+1;}
			}else{
				$n_num=$n_num+1;
			}
		}
		if($print eq 0){next;}
		print WR "\n";
	}
}
close IN;
print "the total base num: $num ; the N base num: $n_num ;\n";
print "the total GATC is: ",$num-$n_num,"\n";
print "GC percent: ",$gc/($num-$n_num),"\n";

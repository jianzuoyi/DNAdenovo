
#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";

my $usage =<<END;

Usage: perl $0 dogma.result > gene.func.lst\n
  This is used to get the function information of chloroplast's coding gene,
  based on dogma results.

END

die $usage if (@ARGV < 1);
my $dogma_file=$ARGV[0];
my $funtion_tab="$Bin/../lib/chloroplast.genefunc.txt";

my (%hash, %hash_2);
open IN1,$funtion_tab or die $!;
while (<IN1>){
	chomp;
	my @line = split(/[\t\n\r]+/,$_);
	if ($line[2] !~ /^NON/){
		$line[2] =~ s/;//g;
		my @gene_id= split(/\s/,$line[2]);
		foreach my $key (@gene_id){
			$hash{$key} = $line[1];		
		}
	}
	$hash{$line[0]} = $line[1];

	$hash_2{$line[0]} = $line[3] if ($line[3] ne "");
}
close IN1;

open IN2,$dogma_file or die $!;
<IN2>;
while (<IN2>){
	chomp;
	my @row = split(/[\t\r\s]+/,$_);
	my $start = $row[0];
	my $end = $row[1];
	my $direction;
	if ($row[-1] == 0){$direction = "-";}
	else{$direction = "+";}
	my $gene_name = $row[2];
	if (exists $hash{$gene_name}) {
		my $h2 = exists $hash_2{$gene_name} ? $hash_2{$gene_name} : ".";
		print "$gene_name\t$start\t$end\t$direction\t$hash{$gene_name}\t$h2\n";
	}
	else {print "$gene_name\t$start\t$end\t$direction\tunknown\t.\n";}
}
close IN2;

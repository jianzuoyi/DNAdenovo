Version: simulate_solexa_reads_0414.cpp
output file:	*_100_170_1.fa ,*_100_170_2.fa;("100" is reads length,"170" is the mean insertsize); 
		snp_solexa_100_170.lis,indel_solexa_100_170.lis((170 is the mean insertsize)if -s and -d have been set);
		*_100_170.log 

example:
*_100_170_1.fa
>read_500_1_1 I 82473 100
TAGAAAAACCAGAGTGGTTGCTGTGGTAGTACGTTGGGGGCTCGTTTGTCGAGGTTTTTAAGCAGGTTTGAATGGAAGAAGAAAAAAACGAGCACTTACT
>read_500_2_1 I 126149 100 95,A;
GCCATTCTAACTCTCCAATTCACGTCCTTGGCTAATTCTGTTATGGCAGGTAACAAAGAGTCTGATAGCAGCTCAATTCCTATTACGTCATTGACTACCT
>read_500_3_1 I 9704 100 74,G;90,T;
CCTTTCATAACTAAACCAATTACAAAAAGTGCAGAAATGTCATGGATACCATTGGCCTTAGATTTTTTAGAAGAAACTGGTAACCGCGAGGAGAGCAACC

in the first line,"500" is the mean insertsize,"1" behind the "500" is the num of reads ,the last "1"
is to show that this read is read 1,"I" is reference id,"82473" is the start position of this read in 
reference,"100" is the length of this read. 
in the third line,the meanning of "read_500_2_1 I 126149 100" is the same as the first line;"95" is 
the error position of this reads and the right base of this position is "A";

snp_solexa_100_170.lis
I       39540   G       A
I       45342   C       T
I       104775  C       T
I       220818  C       G

the first column is the reference id,the second column is the position of this snp,the third column 
is the base in reference and the last column is the snp base;

indel_solexa_100_170.lis
I       -       4280    1       C
I       -       104003  1       G
I       -       81587   1       G
I       +       206841  2       TC

the first column is the reference id,int the second column, "-" show that this is a deletion, "+" thow
that this is a insertion, the third is the position of this indel in reference,the 4th column is the 
bases number of this indel,and the last column is bases that inserted or deleted.

 

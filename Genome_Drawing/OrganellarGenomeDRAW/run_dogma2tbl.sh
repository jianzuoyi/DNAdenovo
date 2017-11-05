#!/bin/bash

DOGMA_RESULT=dogma.result

echo ">Feature scaffold1_pilon
1	152305	REFERENCE
			PubMed	000000F"
grep rrn $DOGMA_RESULT | awk '{if (($4=="-")) {tmp=$1;$1=$2;$2=tmp} {print $1"\t"$2"\tgene";print "\t\t\tgene\t"$3} {print $1"\t"$2"\trRNA";print "\t\t\tproduct\t"$3}}'
grep trn $DOGMA_RESULT | awk '{if (($4=="-")) {tmp=$1;$1=$2;$2=tmp} {print $1"\t"$2"\tgene";print "\t\t\tgene\t"$3} {print $1"\t"$2"\ttRNA";print "\t\t\tproduct\t"$3}}'
grep -v rrn $DOGMA_RESULT | grep -v trn $DOGMA_RESULT | awk '{if (($4=="-")) {tmp=$1;$1=$2;$2=tmp} {print $1"\t"$2"\tgene";print "\t\t\tgene\t"$3} {print $1"\t"$2"\tgene";print "\t\t\tgene\t"$3}}'

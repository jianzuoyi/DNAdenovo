#!/bin/bash

RGI=/its1/GB_BT2/jianzuoyi/biosoft/anaconda/bin/rgi
THREADS=24
OUTDIR=/its1/GB_BT2/jianzuoyi/projects/baumannii/08.CARD

while read genome
do
    name=$(basename $genome .fa)
    echo $name
    OUT_DIR=${OUTDIR}
    echo "$RGI -i $genome -o ${OUTDIR}/$name -n $THREADS -d chromosome"
done < input.fofn

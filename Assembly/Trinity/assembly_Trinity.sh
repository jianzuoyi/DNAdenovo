echo start assembly at `date +%Y-%m-%d' '%H:%M:%S` && \
export PATH=/opt/bio/perl/bin:/its1/GB_BT1/pipeline/RNAdenovo-cloud/Modules/assembly_Trinity/jdk1.7.0_51/bin:/its1/GB_BT1/pipeline/RNAdenovo-cloud/Modules/assembly_Trinity/samtools1/bin:/its1/GB_BT1/pipeline/RNAdenovo-cloud/Modules/assembly_Trinity/bowtie2:/its1/GB_BT1/pipeline/RNAdenovo-cloud/Modules/assembly_Trinity/bowtie:$PATH && \
/its1/GB_BT1/pipeline/RNAdenovo-cloud/Modules/assembly_Trinity/trinityrnaseq-2.0.6/Trinity  --group_pairs_distance 500  --path_reinforcement_distance 75  --min_glue 3  --min_kmer_cov 3  --min_contig_length 100  --seqType fq --left /its1/GB_BT1/pipeline/RNAdenovo-cloud/test.dir/output/filter/B3/B3_1.fq.gz,/its1/GB_BT1/pipeline/RNAdenovo-cloud/test.dir/output/filter/B1/B1_1.fq.gz,/its1/GB_BT1/pipeline/RNAdenovo-cloud/test.dir/output/filter/B2/B2_1.fq.gz --right /its1/GB_BT1/pipeline/RNAdenovo-cloud/test.dir/output/filter/B3/B3_2.fq.gz,/its1/GB_BT1/pipeline/RNAdenovo-cloud/test.dir/output/filter/B1/B1_2.fq.gz,/its1/GB_BT1/pipeline/RNAdenovo-cloud/test.dir/output/filter/B2/B2_2.fq.gz \
--max_memory 50G --CPU 24 --inchworm_cpu 24 \
--bflyHeapSpaceInit 1G --bflyHeapSpaceMax 4G --bfly_opts "-V 5 --edge-thr=0.1 --stderr"  \
--output /its1/GB_BT1/pipeline/RNAdenovo-cloud/test.dir/output/assembly_Trinity_trinity && \
rm -rf /its1/GB_BT1/pipeline/RNAdenovo-cloud/test.dir/output/assembly_Trinity_trinity/both.fa* && \
echo finish assembly at `date +%Y-%m-%d' '%H:%M:%S` && \
echo finish &> /its1/GB_BT1/pipeline/RNAdenovo-cloud/test.dir/shell/assembly_Trinity/assembly_Trinity.sh.ok

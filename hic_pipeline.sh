#!/bin/bash
### input data

#######################
genome="/data/home/ruanlab/huangjiaxiang/pipeline/shell_version/mm10.fa" #genome path that have bwa index
genome_size="/data/home/ruanlab/huangjiaxiang/pipeline/shell_version/mm10.chrom.size"
out_path=$2 #output path
prefix=$1 # prefix of file path
post1fix=_R1.fastq.gz # postfix of R1 file
post2fix=_R2.fastq.gz # postfix of R2 file
NTHREAD=8 # thread number
########################

mkdir ${out_path}
cd ${out_path}

ChIApipeV2="/data/home/ruanlab/huangxingyu/chiaV2.sif"
scom="singularity run ${ChIApipeV2}"
name=`basename ${out_path}`
pbs="$name.hic.jobslurm"
fq1=${prefix}${post1fix}
fq2=${prefix}${post2fix}
bwa="${scom} bwa"
pairtools="${scom} pairtools"
samtools="${scom} samtools"
get_qc="${scom} get_qc"
preseq="${scom} preseq"
bgzip="${scom} bgzip"
juicertool="$scom juicer_tools"
path=${out_path}/result
data=${out_path}/mapped.pairs

echo ${genome}
mkdir result
### fastq to valid pairs
$bwa mem -5SP -T0 -t ${NTHREAD} ${genome} ${fq1} ${fq2} -o aligned.sam
$pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in ${NTHREAD} --nproc-out ${NTHREAD} --chroms-path ${genome} aligned.sam >  parsed.pairsam
$pairtools sort --nproc ${NTHREAD} parsed.pairsam > sorted.pairsam
$pairtools dedup --nproc-in ${NTHREAD} --nproc-out ${NTHREAD} --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam
$pairtools split --nproc-in ${NTHREAD} --nproc-out ${NTHREAD} --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam
$samtools sort -m 2G -@ ${NTHREAD} -T temp.bam -o mapped.PT.bam unsorted.bam
$samtools index mapped.PT.bam
### qc
$get_qc -p stats.txt > result/${name}.out.qc
$preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output result/${name}.out.preseq mapped.PT.bam
## juicetools
awk '!/^#/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}/^#/{print}' mapped.pairs > mapped.pairs.cut
$bgzip mapped.pairs.cut
mv mapped.pairs.cut.gz result/${name}.mapped.pairs.cut.gz -f 
$juicertool pre result/${name}.mapped.pairs.cut.gz result/${name}.hic ${genome_size} -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000
### rm sam
rm aligned.sam parsed.pairsam sorted.pairsam dedup.pairsam unsorted.bam

### Generate heatmap filter reads <1000bp
zcat ${path}/${name}.mapped.pairs.cut.gz|grep "#" > tmppair1
zcat ${path}/${name}.mapped.pairs.cut.gz|grep -v "#" |awk '($2==$4 && $5-$3>1000)||($2!=$4)'> tmppair2
cat tmppair1 tmppair2 | $bgzip > result/${name}.flt1000.mapped.pairs.cut.gz
rm tmppair1 tmppair2
$juicertool pre result/${name}.flt1000.mapped.pairs.cut.gz result/${name}.flt1000.hic ${genome_size} -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000

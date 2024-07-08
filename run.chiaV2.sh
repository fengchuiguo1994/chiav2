#!/usr/bin/bash

set -e 

configfile=$1 #configuration general to a set of libraries
fqfile=$2 # fqpath
RUN=$3 # RUN name

source $configfile
scom="singularity run ${ChIApipeV2}"
maindir=$PWD
pigz="$scom pigz"
samtools="$scom samtools"
bedtools="$scom bedtools"
cpuprog="$scom cpu"
chiasigprog="$scom ChiaSig"
juicertool="$scom juicer_tools"
fasta=$genome
#defined names
pairlabel='singlelinker.paired'
singlabel='singlelinker.single'
nonelabel='none'

#------------------------------------

# read1=${fqfile}${R1suf}
# read2=${fqfile}${R2suf}
nf=$( ls ${fqfile}*${R1suf}| wc -l )
r1=$( ls ${fqfile}*${R1suf} )
r2=$( ls ${fqfile}*${R2suf} )

if [ "$nf" -gt 1 ]; then
   cat ${fqfile}*${R1suf} > ${RUN}/$read1
   cat ${fqfile}*${R2suf} > ${RUN}/$read2
else
    ln -s $r1 ${RUN}/$read1
    ln -s $r2 ${RUN}/$read2
fi

pbs="$RUN.insituChiapet.jobslurm"
peakfile="${RUN}_peaks.narrowPeak"
cis_prefix="${RUN}.e${extbp}.clusters.cis.chiasig"
trans_prefix="${RUN}.e${extbp}.clusters.trans.chiasig"
chiasigoutput="${cis_prefix}.BE2.sigf.interactions"
statout=${RUN}.final_stats.tsv

shortname=${RUN}.chiapet
#Create directory for each RUN


if [ ! -d $RUN ]; then
    echo "Create $RUN"
    mkdir $RUN
fi

cd $RUN

echo "Workdir: $RUN"

LOGFILE=${shortname}.log

echo "Reference genome:"
/bin/ls ${fasta}*

# perform linker detection and generation of different category of fastq files
echo "Linker detection on: ${RUN} " 2>${LOGFILE}

$cpuprog stag -A ${linker} -W -T 18 -t ${NTHREAD} -O ${RUN} ${RUN}/$read1 ${RUN}/$read2  2>>${LOGFILE}

echo "--- linker detection completed ----" >>${LOGFILE}

#Get the stat
$cpuprog stat -s -p -T 18 -t 1 ${RUN}.cpu 2>>${LOGFILE} 1>${RUN}.stat
echo "--- statistics done  ----" >>${LOGFILE}

echo "--- pigziiping   ----" >>${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.singlelinker.paired.fastq 2>>${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.singlelinker.single.fastq 2>>${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.none.fastq 2>>${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.conflict.fastq 2>>${LOGFILE}
$pigz -p ${NTHREAD} ${RUN}.tied.fastq 2>>${LOGFILE}

#------------------------------------------------------------------------------

#Pair

# mapping
echo  START  ${RUN} cpu memaln .. | tee -a ${LOGFILE}
echo  Mapping paired tags .. | tee -a ${LOGFILE}
$cpuprog memaln -T $minMapScore -t ${NTHREAD} $fasta ${RUN}.$pairlabel.fastq.gz 2>> ${LOGFILE} 1>${RUN}.$pairlabel.sam

$pigz -p ${NTHREAD} ${RUN}.$pairlabel.sam | tee -a ${LOGFILE}
echo  ENDED pair mapping | tee -a ${LOGFILE}

#pairing
echo  STARTED ${RUN} cpu pair .. | tee -a ${LOGFILE}
echo  Pairing paired tags .. | tee -a ${LOGFILE}
$cpuprog pair -s $selfbp -S -t ${NTHREAD} ${RUN}.$pairlabel.sam.gz 1>${RUN}.$pairlabel.stat.xls 2>> ${LOGFILE}
echo  ENDED ${RUN} cpu pair .. | tee -a ${LOGFILE}

# span
echo  STARTED ${RUN} cpu span .. | tee -a ${LOGFILE}
echo  Computing span of paired tags .. | tee -a ${LOGFILE}
$cpuprog span -s $selfbp -g -t ${NTHREAD} ${RUN}.$pairlabel.UU.bam 2>> ${LOGFILE} 1>${RUN}.$pairlabel.UU.span.xls
echo  ENDED ${RUN} span pair .. | tee -a ${LOGFILE}

# deduplication
echo  STARTED ${RUN} cpu dedup .. | tee -a ${LOGFILE}
echo  De-duplicating paired tags UU .. | tee -a ${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} ${RUN}.$pairlabel.UU.bam 1>${RUN}.$pairlabel.UU.dedup.lc 2>> ${LOGFILE}
echo  ENDED ${RUN} cpu dedup .. | tee -a ${LOGFILE}

# deduplicated span
echo  STARTED ${RUN} cpu dedup span.. | tee -a ${LOGFILE}
echo  Computing span of paired tags UU nr .. | tee -a ${LOGFILE}
$cpuprog span -s $selfbp -t ${NTHREAD} ${RUN}.$pairlabel.UU.nr.bam 2>> ${LOGFILE} 1>${RUN}.$pairlabel.UU.nr.span.xls
echo  ENDED ${RUN} cpu dedup span.. | tee -a ${LOGFILE}

# cluster tags
echo  STARTED ${RUN} cpu clustering.. | tee -a ${LOGFILE}
$cpuprog cluster -s $selfbp -M -B 1000 -5 5,-20 -3 5,480 -t ${NTHREAD} -O ${RUN}.e$extbp -j -x -v 1 -g ${RUN}.$pairlabel.UU.nr.bam 2>&1 | tee -a ${LOGFILE}
echo  ENDED ${RUN} cpu clustering.. | tee -a ${LOGFILE}

#------------------------------------------------------------------------------
#None tag

# mapping
echo  STARTED ${RUN}.$nonelabel cpu memaln .. | tee -a ${LOGFILE}
$cpuprog memaln -T 15 -t ${NTHREAD} $fasta ${RUN}.$nonelabel.fastq.gz 2>> ${LOGFILE} 1>${RUN}.$nonelabel.sam
$pigz -p ${NTHREAD} ${RUN}.$nonelabel.sam 2>&1  | tee -a ${LOGFILE}
echo  ENDED ${RUN} cpu memaln .. | tee -a ${LOGFILE}

# pairing
echo Pairing $nonelabel tags .. | tee -a ${LOGFILE}
$cpuprog pair -S -t ${NTHREAD} ${RUN}.$nonelabel.sam.gz 1>${RUN}.$nonelabel.stat.xls 2>> ${LOGFILE}
echo  ENDED ${RUN} cpu pair .. | tee -a ${LOGFILE}

# span
echo  STARTED ${RUN}.$nonelabel cpu span .. | tee -a ${LOGFILE}
$cpuprog span -g -t ${NTHREAD} ${RUN}.$nonelabel.UU.bam 2>> ${LOGFILE} 1>${RUN}.$nonelabel.UU.span.xls
echo  ENDED ${RUN}.$nonelabel span pair .. | tee -a ${LOGFILE}

# deduplication
echo  STARTED ${RUN}.$nonelabel cpu dedup .. | tee -a ${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} ${RUN}.$nonelabel.UU.bam 1>${RUN}.$nonelabel.UU.dedup.lc 2>> ${LOGFILE}

#------------------------------------------------------------------------------

#1tag

# mapping
echo STARTED ${RUN}.$singlabel cpu memaln .. | tee -a ${LOGFILE}
$cpuprog memaln -T 15 -t ${NTHREAD} $fasta ${RUN}.$singlabel.fastq.gz 2>> ${LOGFILE} 1>${RUN}.$singlabel.sam
$pigz -p ${NTHREAD} ${RUN}.$singlabel.sam 2>&1  | tee -a ${LOGFILE}

# pairing to get bam
echo Pairing $singlabel tags .. | tee -a ${LOGFILE}
$cpuprog pair -S -t ${NTHREAD} ${RUN}.$singlabel.sam.gz 1>${RUN}.$singlabel.stat.xls 2>> ${LOGFILE}
$scom samtools view -q ${minMapScore} -o ${RUN}.$singlabel.UxxU.q30.bam ${RUN}.$singlabel.UxxU.bam
mv ${RUN}.$singlabel.UxxU.q30.bam ${RUN}.$singlabel.UxxU.bam

# deduplication
echo STARTED ${RUN}.$singlabel cpu dedup .. | tee -a ${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} ${RUN}.$singlabel.UxxU.bam 1>${RUN}.$singlabel.UxxU.dedup.lc 2>> ${LOGFILE}
echo  ENDED ${RUN}.$singlabel cpu dedup .. | tee -a ${LOGFILE}

#------------------------------------------------------------------------------
#Filter out non primary reads
$samtools view -F 2048 -@ $NTHREAD -h ${RUN}.$pairlabel.UU.nr.bam  | awk 'length($10) > 30 || $1 ~ /^@/' \
        | $samtools sort -@ $NTHREAD - -o ${RUN}.$pairlabel.F2048.bam
$samtools view -F 2048 -@ $NTHREAD -h ${RUN}.$singlabel.UxxU.nr.bam | awk 'length($10) > 30 || $1 ~ /^@/' \
        | $samtools sort -@ $NTHREAD - -o ${RUN}.$singlabel.F2048.bam
$samtools view -F 2048 -@ $NTHREAD -h ${RUN}.$nonelabel.UU.nr.bam  | awk 'length($10) > 30 || $1 ~ /^@/' \
        | $samtools sort -@ $NTHREAD - -o ${RUN}.$nonelabel.F2048.bam

echo -e "Converting file formats..\n" >> ${LOGFILE}
$samtools sort -@ $NTHREAD -o ${RUN}.$pairlabel.nr.sorted.bam ${RUN}.$pairlabel.F2048.bam
$samtools sort -@ $NTHREAD -o ${RUN}.$singlabel.nr.sorted.bam ${RUN}.$singlabel.F2048.bam    
$samtools sort -@ $NTHREAD -o ${RUN}.$nonelabel.nr.sorted.bam ${RUN}.$nonelabel.F2048.bam  

#forBASIC still have blacklist.  Removed secondary alignment with -F2048; for clean view we have NDP ${RUN}_treat_pileup.clip.sorted.bdg
$samtools merge -@ $NTHREAD ${RUN}.forBASIC.bam ${RUN}.$pairlabel.nr.sorted.bam ${RUN}.$singlabel.nr.sorted.bam ${RUN}.$nonelabel.nr.sorted.bam

echo 'Convert bam to bed for macs2 callpeak and clean from blacklist... '
if [ ${blackbed} != "null" ];then
        $bedtools bamtobed -i ${RUN}.$pairlabel.nr.sorted.bam | $bedtools subtract -a stdin -b $blackbed > ${RUN}.$pairlabel.nr.sorted.bed
        $bedtools bamtobed -i ${RUN}.$singlabel.nr.sorted.bam | $bedtools subtract -a stdin -b $blackbed > ${RUN}.$singlabel.nr.sorted.bed
        $bedtools bamtobed -i ${RUN}.$nonelabel.nr.sorted.bam | $bedtools subtract -a stdin -b $blackbed > ${RUN}.$nonelabel.nr.sorted.bed
else
        $bedtools bamtobed -i ${RUN}.$pairlabel.nr.sorted.bam > ${RUN}.$pairlabel.nr.sorted.bed
        $bedtools bamtobed -i ${RUN}.$singlabel.nr.sorted.bam > ${RUN}.$singlabel.nr.sorted.bed
        $bedtools bamtobed -i ${RUN}.$nonelabel.nr.sorted.bam > ${RUN}.$nonelabel.nr.sorted.bed
fi
ls $RUN.*.nr.sorted.bed


$scom macs2 callpeak --keep-dup all --nomodel --shift $shiftsize  --extsize $peakext  -B --SPMR -t ${RUN}.$pairlabel.nr.sorted.bed ${RUN}.$singlabel.nr.sorted.bed ${RUN}.$nonelabel.nr.sorted.bed  -f BED -g $organism -n $RUN  --qvalue $macsq >> ${LOGFILE}

echo  Generating coverage density.. | tee -a ${LOGFILE}
#$samtools view -H ${RUN}.$pairlabel.UU.nr.bam | grep '^@SQ' | cut -f 2-3 | sed s?SN:?? | sed s?LN:?? > ${RUN}.genome.length
$samtools view -H ${RUN}.forBASIC.bam | grep '^@SQ' | cut -f 2-3 | sed s?SN:?? | sed s?LN:?? > ${RUN}.genome.length

$scom bedClip ${RUN}_treat_pileup.bdg ${RUN}.genome.length ${RUN}_treat_pileup.clip.bdg
sort -k1,1 -k2,2n ${RUN}_treat_pileup.clip.bdg > ${RUN}_treat_pileup.clip.sorted.bdg

$scom bedGraphToBigWig ${RUN}_treat_pileup.clip.sorted.bdg ${RUN}.genome.length ${RUN}.treat_pileup.NDP.bw
echo ENDED ${RUN} coverage density generated. | tee -a ${LOGFILE}

# Make bedgraph
$bedtools genomecov -ibam ${RUN}.forBASIC.bam -bg > ${RUN}.forBASIC.bedgraph

# Sort bedgraph
$bedtools sort -i ${RUN}.forBASIC.bedgraph > ${RUN}.forBASIC.sorted.bedgraph

# Make bigwig
$scom bedGraphToBigWig ${RUN}.forBASIC.sorted.bedgraph  ${RUN}.genome.length ${RUN}.coverage.forBASIC.bw
echo -e  "Done converting file formats, ${RUN}.coverage.forBASIC.bw generated.\n" >> ${LOGFILE}

#------------------------------------------------------------------------------
#Significant interaction
#get ipet > 1
echo  Creating BE2 file and clean blacklist using $blackbed  >> ${LOGFILE}
zcat ${cis_prefix}.gz | awk '{ if ( $7 > 1 ) print }' > ${cis_prefix}.BE2

#run chiasig
$chiasigprog -c $PET -t $NTHREAD -p -m $selfbp ${cis_prefix}.BE2

# Heatmap for Juicebox

sed 's/chr//' ${RUN}.genome.length | sed 's/^M/ZZZZ/' | sort -k1,1 -V | sed 's/^ZZZZ/M/' >  ${RUN}.nochr.genome.length
$juicertool pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000 -k KR,GW_KR ${RUN}.e${extbp}.juice.gz $RUN.hic  ${RUN}.nochr.genome.length  2>> ${LOGFILE}
rm  ${RUN}.nochr.genome.length

$pigz -p ${NTHREAD} $RUN.*.bedgraph
$pigz -p ${NTHREAD} $RUN.*.bdg #keep the clip.sorted.bdg

# Get the summary
echo Starting to sumarize $PWD
echo "Output file: $statout" 

cat ${RUN}.stat | grep "Filename" | awk -F'[ \t.]' '{print "Library_ID\t"$2}' > ${statout}
echo -e "Reference_genome\t"${genome} >> ${statout}

## PET count
# Get PET count
n_pet=$( cat ${RUN}.stat | grep "Total pairs" | awk -F'[ \t]' '{print $3}' )
n_pet=$( printf "%'.f" ${n_pet} )
echo -e "Total_PE_reads\t"${n_pet} >> ${statout}
## Get PET count with linker
read_link=$( cat ${RUN}.stat | grep "Linker detected" | awk -F '[ \t]' '{print $3}' | xargs printf "%'.f")
echo -e "Linker_detected_in_pair_reads\t"${read_link} >> ${statout}
## Get Pai count with linker
pet_link=$( cat ${RUN}.stat | grep "Single Linker 2 tags" | awk -F '[ \t]' '{print $6}' | xargs printf "%'.f")
echo -e "PET_with_linker\t"${pet_link} >> ${statout}

## Mapping
#cd $cpuoutdir
# Get mapped PET count 
# Get uniquely mapped PET count 
unique=$( cat ${RUN}.$pairlabel.UU.span.xls | grep "Total pairs" | awk -F '[\t]' '{print $2}' )
# Get uniquely mapped and non-redundant PET count 
nr=$( cat ${RUN}.$pairlabel.UU.nr.span.xls | grep "Total pairs" | awk -F '[\t]' '{print $2}' )
# Compute redundancy
redun=$( echo "(${unique} - ${nr}) / ${unique}" | bc -l )
# Write uniquely mapped PET count
unique=$( printf "%'.f" ${unique} )
echo -e "Uniquely_mapped_PET\t"${unique} >> ${statout}
# Write unique mapped and non-redundant PET count
nr=$( printf "%'.f" ${nr} )
echo -e "Non-redundant_PET\t"${nr} >> ${statout}
# Write redundancy
redun=$( printf %.2f ${redun} )
echo -e "Redundancy\t"${redun} >> ${statout}

## Get number of peaks
#cd $peakdir
np=$( cat ${RUN}_peaks.narrowPeak | wc -l )
n_peak=$( printf "%'.f" ${np} )
echo -e "Peak\t"$n_peak >> ${statout}

## Interaction types
# Get self-ligation PET count
self_lig=$( cat ${RUN}.$pairlabel.UU.nr.span.xls | grep "second/best<0.95" -A5 | awk -F '[\t]' '{if(NR==4)print $2}' )
self_lig=$( printf "%'.f" ${self_lig} )
echo -e "Self-ligation_PET\t"${self_lig} >> ${statout}
# Get inter-ligation PET count (intra-chr)
intra_chr_pet=$( cat ${RUN}.$pairlabel.UU.nr.span.xls | grep "second/best<0.95" -A5 | awk -F '[\t]' '{if(NR==5)print $2}' )
# Get inter-ligation PET count (inter-chr)
inter_chr_pet=$( cat ${RUN}.$pairlabel.UU.nr.span.xls | grep "second/best<0.95" -A5 | awk -F '[\t]' '{if(NR==2)print $2}' )
# Compute ratio of intra-chr to inter-chr inter-ligation PETs
pet_ratio=$( echo "${intra_chr_pet} / ${inter_chr_pet}" | bc -l )
# Compute inter-ligation PET count (all)
inter_lig_all=$( echo "${intra_chr_pet} + ${inter_chr_pet}" | bc )
# Write inter-ligation PET count (all)
inter_lig_all=$( printf "%'.f" ${inter_lig_all} )
echo -e "Inter-ligation_PET\t"${inter_lig_all} >> ${statout}
# Write inter-ligation PET count (intra-chr)
intra_chr_pet=$( printf "%'.f" ${intra_chr_pet} )
echo -e "Intra-chr_PET\t"${intra_chr_pet} >> ${statout}
# Write inter-ligation PET count (inter-chr)
inter_chr_pet=$( printf "%'.f" ${inter_chr_pet} )
echo -e "Inter-chr_PET\t"${inter_chr_pet} >> ${statout}
# Write ratio of intra-chr to inter-chr inter-ligation PETs
pet_ratio=$( printf %.2f ${pet_ratio} )
echo -e "ratio_of_intra/inter_PET\t"${pet_ratio} >> ${statout}

## Singleton
# Get singleton PET count (all)
singleton=$(zcat *chiasig.gz | awk '$7==1{print}' | wc -l)
singleton=$( printf "%'.f" ${singleton} )
echo -e "Singleton\t"$singleton >> ${statout}
# Get singleton PET count (intra-chr)
intra_singleton=$(zcat *cis.chiasig.gz | awk '$7==1{print}' | wc -l)
intra_singleton=$( printf "%'.f" ${intra_singleton} )
echo -e "Intra-chr_singleton\t"$intra_singleton >> ${statout}
# Get singleton PET count (inter-chr)
inter_singleton=$(zcat *trans.chiasig.gz | awk '$7==1{print}' | wc -l)
inter_singleton=$( printf "%'.f" ${inter_singleton} )
echo -e "Inter-chr_singleton\t"$inter_singleton >> ${statout}

## Clusters (overall)
# Get cluster count
total_cluster_number=$(zcat *chiasig.gz | awk '$7 != 1{print}' | wc -l)
total_cluster_number=$( printf "%'.f" ${total_cluster_number} )
echo -e "PET_cluster\t"${total_cluster_number} >> ${statout}
# Get intra-chr cluster count
intra_cluster=$( zcat *cis.chiasig.gz | awk '$7 >=2 {print}' | wc -l )
# Get inter-chr cluster count
inter_cluster=$( zcat *trans.chiasig.gz | awk '$7 >=2 {print}' | wc -l)

# Compute ratio of intra-chr to inter-chr clusters
cluster_ratio=$( echo "${intra_cluster} / ${inter_cluster}" | bc -l )
cluster_ratio=$( printf %.2f ${cluster_ratio} )

# Write cluster ratio
echo -e "ratio_of_intra/inter_cluster\t"${cluster_ratio} >> ${statout}
## Clusters (intra-chr)

# Write intra-chr cluster count
intra_cluster=$( printf "%'.f" ${intra_cluster} )
echo -e "Intra-chr_PET_cluster\t"${intra_cluster} >> ${statout}

# Get intra-chr cluster count by number of PETs (1 - 10)
for i in $(seq 2 10)
do
    intra_pets_number=$(zcat *cis.chiasig.gz | awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | xargs printf "%'.f")
    echo -e "pets_number_"${i}"\t"${intra_pets_number} >> ${statout}
done

# Get intra-chr cluster count with > 10 PETs
echo -e "pets_number>10\t"$(zcat *cis.chiasig.gz | awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${statout}

## Clusters (inter-chr)
# Write inter-chr cluster count
inter_cluster=$( printf "%'.f" ${inter_cluster} )
echo -e "Inter-chr_PET_cluster\t"${inter_cluster} >> ${statout}

# Get inter-chr cluster count by number of PETs (1 - 10)
for i in $(seq 2 10)
do
    inter_pets_number=$(zcat *trans.chiasig.gz | awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | xargs printf "%'.f")
    echo -e "pets_number_"${i}"\t"${inter_pets_number} >> ${statout}
done

# Get inter-chr cluster count with > 10 PETs
echo -e "pets_number>10\t"$(zcat *trans.chiasig.gz | awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${statout}

#significant interaction 
cis_prefix="${RUN}.e${extbp}.clusters.cis.chiasig"
chiasigoutput="${cis_prefix}.BE2.sigf.interactions"
ns=$( grep -v : $chiasigoutput | wc -l )
ns=$( printf "%'.f" ${ns} )
echo -e "Significant_interaction>=${PET} \t"$ns >> ${statout}

# Cleaning files
# /bin/rm ${RUN}.*.fastq.gz 2>> ${LOGFILE}
# /bin/rm ${RUN}.*.cluster 2>> ${LOGFILE}
# /bin/rm ${RUN}.*.dedup 2>> ${LOGFILE}
# /bin/rm ${RUN}.*.nr.bed 2>> ${LOGFILE}
# /bin/rm ${RUN}_treat_pileup.clip.bdg ${RUN}_treat_pileup.bdg  2>> ${LOGFILE}
# /bin/rm ${RUN}.*dup.bam ${RUN}.*.self.*.bam ${RUN}.*.discordant.bam ${RUN}.*.nn.bam ${RUN}.*.xx.bam 2>> ${LOGFILE}
# /bin/rm $RUN.*.sam.gz $RUN.*.intra.bam $RUN.*.inter.bam $RUN.*.failed.bam 
# /bin/rm $RUN.bedpe.anchors.bed $RUN.cis.sigf.anchors.intersect.bed $RUN.bedpe.anchors.gene.bed
# /bin/rm ${RUN}.forBASIC.bedgraph*
# /bin/rm ${RUN}_control_lambda.bdg
# rm/bin/ $RUN.pooled.clusters.gz

echo "CPU pipeline run $RUN  completed!"

cd $maindir 
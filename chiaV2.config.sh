#!/bin/bash

## Config file for CPU pipeline

#Define directories & filenames
R1suf="_R1.fastq.gz" # suffix of the fastq files, modify if different
R2suf="_R2.fastq.gz" # read2 suffix

#Chiapet params
extbp=500
selfbp=8000
minMapScore=30 #this is the default CPU memaln set

#genome='mm10' #genome version

genome="/data/home/ruanlab/huangjiaxiang/pipeline/shell_version/mm10.fa"
blackbed="null"
ChIApipeV2="/data/home/ruanlab/huangjiaxiang/software/chiav2.sif"
#Peak calling params
organism='mm' #used in macs: human='hs'  mouse='mm'  fly='dm'
macsq=0.000001
#macsq=0.01
peakext=150
shiftsize=$(( -$peakext/2 ))

#Chiasig params
PET=3

### linker seq
linker=ACGTGATATTTCACGACTCT ## linker 3

#cluster job setting parameters:
NTHREAD=10 #cpu usage per job
hrs=72 #number of hours allocation
GB=72 #memory in GB
seed=12639

#matrix params
binsize=100000 # matrix resolution in bp


## Not important for running, metadata only

# The type of sequencing run:
#    "miseq" - around 30 million reads
#    "nextseq" - around 300 million reads
#    "hiseq" - around 300 million reads
#    "pooled" - around 1 billion reads
run_type="HiMi"

# The factor for which the IP was performed
ip_factor="ATAC"

# Cell type
cell_type="NA"




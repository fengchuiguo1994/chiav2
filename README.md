# chiav2

ChiaV2 is a packaged version of the [ChIA-PIPE workflow](https://github.com/TheJacksonLaboratory/ChIA-PIPE), which eliminates the hassle of software installation. <br/>
**DockerHub Link: https://hub.docker.com/r/fengchuiguo1994/chiav2/tags**

##  Introduction
chiav2 processes the sequencing data of ChIATAC/ChIA-PET to generate chromatin loops, coverage signal and contact matrix. chiaV2 includes linker calling, mapping, duplicate removing, loop calling.

##  [System Requirements](https://sylabs.io/guides/3.8/admin-guide/installation.html)
###   Hardware
ChiaV2 should be run on a Linux system that meets the following requirements:
* 4-core Intel or AMD processor (>12 cores recommended)
* 64GB RAM (>256GB recommended)
* 64-bit CentOS/RedHat 7.8 or Ubuntu 20.04

###   Software
* Singularity: a container platform, >=3.8

####   Quick installation of Singularity
```
## On Red Hat Enterprise Linux or CentOS install the following dependencies:
$ sudo yum update -y && \
     sudo yum groupinstall -y 'Development Tools' && \
     sudo yum install -y \
     openssl-devel \
     libuuid-devel \
     libseccomp-devel \
     wget \
     squashfs-tools \
     cryptsetup

## On Ubuntu or Debian install the following dependencies:
$ sudo apt-get update && sudo apt-get install -y \
    build-essential \
    uuid-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup-bin

## Install Go
$ export VERSION=1.14.12 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz

$ echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

## Install singularity on CentOS without compile
$ yum install -y singularity
```

####   Quick download chiaV2 from DockerHub
You can download chiav2 by running the following command:
```
singularity build chiav2.sif docker://fengchuiguo1994/chiav2:1.0
```


### Usage
A ChIATAC demo data from mice was provided.
```
# 1. Prepare genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.genome.fa.gz
gunzip GRCm39.genome.fa.gz
singularity run chiav2.sif bwa index GRCm39.genome.fa

# 2. Modify the information in the config file.

# 3. Execute program.
bash run.chiaV2.sh chiaV2.config.sh demoData/SCG0192 result
```

###   Example result
| file name | description |
| ----------- | ----------- |
| *.narrowPeak | DNA binding region (like ChIPSeq) |
| *.treat_pileup.NDP.bw | coverage signal (like ChIPSeq) |
| *.hic | contact matrix (Visualization in [juicer-box](https://aidenlab.org/juicebox/))|
| *.BE2.sigf.interactions | chromatin loops |

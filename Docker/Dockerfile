# Use an ubuntu 20.04 LTS as OS
FROM ubuntu:20.04

# Set the working directory to /app
WORKDIR /app

RUN mkdir clinsv 

ADD clinsvScripts.tar.gz /app
ADD installKNITR.r /app

# Enviornment variables
ENV clinsv_path "/app/clinsv"

# Move all contents of clinsv scripts to the clinsv folder
RUN mv clinsvScripts/* clinsv/
RUN rm -r clinsvScripts

ENV PATH="/app/clinsv/bin:${PATH}"

# Get links to repos
RUN apt-get update

# To stop essential asking for input
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Australia/Sydney

RUN apt-get -y install build-essential

RUN apt-get -y install wget curl git nano

# Dependencies for Tabix and Big perl modules
RUN apt-get -y install zlib1g-dev 

# Samtools dependencies
RUN apt-get -y install libbz2-dev liblzma-dev libcurl4-openssl-dev libncurses5-dev

# Python 3.8 dependencies
RUN apt-get -y install python3.8 python3-pip

# Root v6 dependencies
RUN apt-get -y install dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev libssl-dev

# Lumpy dependencies
RUN apt-get -y install dh-autoreconf

# Change sh command to use bash rather than dash from bin.
RUN ln -sf ../bin/bash ../bin/sh

# Install cpanm
RUN curl -L https://cpanmin.us | perl - App::cpanminus

RUN mkdir clinsv/perl

# Install Bio::DB::Big perl module
RUN cd $clinsv_path/perl && \
    git clone https://github.com/dpryan79/libBigWig && \
    cd libBigWig/ && \
    make && \
    export LIBBIGWIG_DIR=$PWD && \
    cpanm install Bio::DB::Big

# Install samtools bcftools and htslib 
RUN cd $clinsv_path && \
    mkdir samtools && \
    cd samtools && \
    wget https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 && \
    wget https://github.com/samtools/bcftools/releases/download/1.4/bcftools-1.4.tar.bz2 && \
    wget https://github.com/samtools/htslib/releases/download/1.4/htslib-1.4.tar.bz2 && \
    tar xvfj samtools-1.4.tar.bz2 && \
    tar xvfj bcftools-1.4.tar.bz2 && \
    tar xvfj htslib-1.4.tar.bz2 && \
    cd htslib-1.4 && \
    make && \ 
    make prefix=$clinsv_path/samtools install && \
    cd ../samtools-1.4 && \
    make && \ 
    make prefix=$clinsv_path/samtools install && \
    cd ../bcftools-1.4 && \
    make && \
    make prefix=$clinsv_path/samtools install && \
    cd .. && \
    rm -r htslib* samtools* bcftools*

# Create symbolic links to all the samtools  tools
RUN cd $clinsv_path/bin && \
    ln -s ../samtools/bin/bcftools bcftools && \
    ln -s ../samtools/bin/samtools samtools && \
    ln -s ../samtools/bin/bgzip bgzip && \
    ln -s ../samtools/bin/tabix tabix 

# Install bigwigtools
RUN mkdir $clinsv_path/bigWigTools && \
    cd $clinsv_path/bigWigTools && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigToWig && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig && \
    chmod ug+x $clinsv_path/bigWigTools/* 

RUN cd $clinsv_path/bin && \
    ln -s ../bigWigTools/bigWigToWig bigWigToWig && \
    ln -s ../bigWigTools/wigToBigWig wigToBigWig

# Install python2
RUN mkdir $clinsv_path/python

RUN cd $clinsv_path/python && \
    wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz && \
    tar zxf Python-2.7.18.tgz && \
    cd Python-2.7.18 && \
    ./configure  --prefix=$clinsv_path/python && \
    make && \
    make install

RUN cd $clinsv_path/python && \
    rm -r Python-2.7.18*

RUN cd $clinsv_path/bin && \
    ln -s ../python/bin/python python && \
    ln -s ../python/bin/python python2

# Install python modules
RUN pip3 install numpy
RUN pip3 install pysam

#installing root 6 (root 5 deprecated)
RUN cd $clinsv_path && \
    git clone --branch latest-stable https://github.com/root-project/root.git root_src && \
    mkdir root_build root_install && cd root_build && \
    cmake -DCMAKE_INSTALL_PREFIX=../root_install ../root_src && \
#add -j4 perhaps, means amount of cores available to build
    cmake --build . -- install

SHELL ["/bin/bash", "-c"] 

RUN mv $clinsv_path/root_install $clinsv_path/root

RUN mkdir $clinsv_path/cnvnator-multi && \
    cd $clinsv_path/cnvnator-multi && \
    git clone --recursive https://github.com/hall-lab/speedseq 
    
RUN cd $clinsv_path/cnvnator-multi/speedseq && \
#need to remove "-std=c++11" from /app/clinsv/cnvnator-multi/speedseq/src/CNVnator/Makefile , prior to running make
    sed -i 's/\-std=c++11//g' /$clinsv_path/cnvnator-multi/speedseq/src/CNVnator/Makefile && \
    source $clinsv_path/root/bin/thisroot.sh && \
    make cnvnator && \
    cd $clinsv_path/cnvnator-multi && \
    ln -s speedseq/bin/cnvnator cnvnator-multi && \
    ln -s speedseq/bin/cnvnator_wrapper.py . 

# Install Lumpy
RUN cd $clinsv_path && \
    git clone --recursive https://github.com/arq5x/lumpy-sv.git && \
    cd lumpy-sv && \
    make

RUN cd $clinsv_path/bin && \
    ln -s ../lumpy-sv/bin/lumpy lumpy

# Install tabix perl module
RUN cd $clinsv_path/perl && \
    git clone https://github.com/MinocheAE/tabix.git && \
    cd tabix && \
    make && \
    cd perl && \
    perl Makefile.PL && \
    make && \
    make install 

# Install perl modules required for CNVnator section
RUN cpanm Compress::Zlib && \
    cpanm Sort::Naturally && \
    cpanm List::BinarySearch::XS && \
    cpanm Excel::Writer::XLSX && \
    cpanm LWP

RUN apt install -y r-base && \
    R -f installKNITR.r

# Install texlive
RUN apt -y install texlive-latex-extra

# Add entry point
ENTRYPOINT [ "perl", "/app/clinsv/bin/clinsv" ]

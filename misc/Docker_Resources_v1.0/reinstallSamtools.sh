#!/bin/bash

clinsv_path=$PWD/clinsv

#get rid of existing samtools
rm -rf clinsv/samtools

cd $clinsv_path
mkdir samtools
cd samtools
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2

tar xvfj samtools-1.10.tar.bz2
tar xvfj bcftools-1.10.2.tar.bz2
tar xvfj htslib-1.10.2.tar.bz2

cd htslib-1.10.2
make; make prefix=$clinsv_path/samtools install

cd ../samtools-1.10
make; make prefix=$clinsv_path/samtools install

cd ../bcftools-1.10.2
make; make prefix=$clinsv_path/samtools install

cd ..
rm -rf htslib* samtools* bcftools*

cd $clinsv_path/bin
ln -sf ../samtools/bin/bcftools bcftools
ln -sf ../samtools/bin/samtools samtools
ln -sf ../samtools/bin/bgzip bgzip
ln -sf ../samtools/bin/tabix tabix
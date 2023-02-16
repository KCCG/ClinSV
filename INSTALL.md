# Installation of ClinSV from source
We recommend installation from Singularity or Docker images. 

This section provides an overview of the dependencies and assistence to compile them from source.

If using local versions of R, Perl or Python is preferred, remove the corresponding symlink in `$clinsv_path/bin`. Make sure that the required packages are installed.

## Create a ClinSV installation directory
```
clinsv_path=$PWD/clinsv
mkdir -p $clinsv_path/bin
export PATH=$clinsv_path/bin:$PATH
```


## Get the latest ClinSV version
* excluding precompiled dependencies

1) pre-publication access
```
wget https://nci.space/clinsv/clinsv_b37/ClinSV_scripts_only_v0.9.tar.gz
tar zxf ClinSV_scripts_only_v0.9.tar.gz
```
2) post-publication access
```
git clone https://github.com/KCCG/ClinSV.git
```

## Install perl and dependencies
### Install perl
```
mkdir $clinsv_path/perl
cd $clinsv_path/perl

wget http://www.cpan.org/src/5.0/perl-5.24.1.tar.gz
tar -xzf perl-5.24.1.tar.gz
cd perl-5.24.1
./Configure -des -Dprefix=$clinsv_path/perl -Duserelocatableinc
make
make test
make install

cd $clinsv_path/perl
rm -r perl-5.24.1*
```

### Install cpanm
```
cd $clinsv_path/perl
curl -L https://cpanmin.us/ -o cpanm
chmod +x cpanm
```

### create links to bin
```
cd $clinsv_path/bin
ln -s ../perl/bin/perl perl
ln -s ../perl/cpanm cpanm
```

### Install required perl modules with cpanm
```
cpanm Bio::Perl
cpanm Compress::Zlib
cpanm Sort::Naturally
cpanm List::BinarySearch::XS
cpanm Excel::Writer::XLSX
cpanm LWP
```

### Install Tabix perl module
```
cd $clinsv_path/perl
git clone https://github.com/MinocheAE/tabix.git
cd tabix
make
cd perl
perl Makefile.PL
make
make test
make install
```

### Install Bio::DB::BigWig perl module
```
cd $clinsv_path/perl
wget https://github.com/ENCODE-DCC/kentUtils/archive/v286.0.0.tar.gz
tar zxf v286.0.0.tar.gz

cd kentUtils-286.0.0/src/lib
export MACHTYPE=x86_64
make CFLAGS='-fPIC' > /dev/null
cd ..
export KENT_SRC=`pwd`

cpanm Bio::DB::BigWig
```

## Install Samtools BCFtools Tabix BGZip

```
cd $clinsv_path
mkdir samtools
cd samtools
wget https://github.com/samtools/htslib/releases/download/1.4/htslib-1.4.tar.bz2
wget https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2
wget https://github.com/samtools/bcftools/releases/download/1.4/bcftools-1.4.tar.bz2
tar xvfj htslib-1.4.tar.bz2
tar xvfj samtools-1.4.tar.bz2
tar xvfj bcftools-1.4.tar.bz2

cd htslib*
make; make prefix=$clinsv_path/samtools install

cd ../samtools*
make; make prefix=$clinsv_path/samtools install

cd ../bcftools*
make; make prefix=$clinsv_path/samtools install

cd ..
rm -r htslib* samtools* bcftools*

cd $clinsv_path/bin
ln -s ../samtools/bin/bcftools bcftools
ln -s ../samtools/bin/samtools samtools
ln -s ../samtools/bin/bgzip bgzip
ln -s ../samtools/bin/tabix tabix

```

## Install BigWigTools

For source and binaries please refer to:
https://github.com/ENCODE-DCC/kentUtils

Save the binaries into: `$clinsv_path/bigWigTools`

```
cd $clinsv_path/bin
ln -s ../bigWigTools/bigWigToWig bigWigToWig
ln -s ../bigWigTools/wigToBigWig wigToBigWig
```



## Install Python

```
mkdir $clinsv_path/python
cd $clinsv_path/python

wget https://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz
tar zxf Python-2.7.11.tgz
cd Python-2.7.11
./configure  --prefix=$clinsv_path/python
make && make install

rm -r Python-2.7.11*
cd $clinsv_path/bin
ln -s ../python/bin/python python
ln -s ../python/bin/pip pip
```

## install python modules
```
pip install numpy
pip install pysam
```



## Install CNVnator

For installing CNVnator
https://github.com/hall-lab/speedseq#installation



### Install root
```
cd $clinsv_path

wget https://root.cern.ch/download/root_v5.34.38.source.tar.gz
tar -zxf root_v5.34.38.source.tar.gz
cd root
./configure
make

source $PWD/bin/thisroot.sh
```

### Install CNVnator-multi
```
mkdir $clinsv_path/cnvnator-multi
cd $clinsv_path/cnvnator-multi

git clone --recursive https://github.com/hall-lab/speedseq
cd speedseq
make cnvnator
```

* The `-std=c++11` may need to be removed from `src/CNVnator/Makefile`
```
cd $clinsv_path/cnvnator-multi
ln -s speedseq/bin/cnvnator cnvnator-multi .
ln -s speedseq/bin/cnvnator_wrapper.py .
```

## Install Lumpy
```
cd $clinsv_path
git clone --recursive https://github.com/arq5x/lumpy-sv.git
cd lumpy-sv
make

cd $clinsv_path/bin
ln -s ../lumpy-sv/bin/lumpy lumpy
```


## Install R
```
mkdir $clinsv_path/R
cd $clinsv_path/R
wget http://cran.rstudio.com/src/base/R-3/R-3.1.1.tar.gz
tar zxf R-3.1.1.tar.gz
cd R-3.1.1
./configure --prefix=$clinsv_path/R
make && make install
rm -fr $clinsv_path/R/R-3.1.1*

cd $clinsv_path/bin
ln -s ../R/bin/R R
```

### install R package knitr locally
```
cd $clinsv_path/R/lib64/R/library
R
.libPaths(getwd())
install.packages("knitr", getwd(),repos='http://cran.rstudio.org');
library("knitr");
quit(save="no")
```

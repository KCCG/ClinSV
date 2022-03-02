# Docker files

## Base Image
The docker image is based on Ubuntu:20.04 LTS.

## Dependencies 
The Docker images uses these specific depenedencies:

- Python v3.8
- Python v2.7.18
- Root v6 
- Lumpy
- Perl
- LibBigWig
- Samtools v1.4
- Bcftool v1.4
- Htslib v1.4
- Tabix perl module
- bigwigtools
- CNVnator

These dependecies also have these dependencies:

**For Tabix and Big perl modules**
- zlib1g-dev

**For Samtools v1.4**
- libbz2-dev 
- liblzma-dev 
- libcurl4-openssl-dev 
- libncurses5-dev

**For Python**
- python v2.7.18
- python v3.8

  **Python 3 modules**
  - numpy
  - pysam

**For root v6**
- dpkg-dev 
- cmake 
- g++ 
- gcc 
- binutils 
- libx11-dev 
- libxpm-dev 
- libxft-dev 
- libxext-dev 
- libssl-dev

**For Lumpy**
- dh-autoreconf

**For Perl**
- cpanm
 
  **Perl modules**
  - Compress::Zlib
  - Sort::Naturally 
  - List::BinarySearch::XS
  - Excel::Writer::XLSX
  - LWP
  - Bio::DB::Big


# Use an CentOS 8 as OS
FROM centos:centos8

# Set the working directory to /app
WORKDIR /app

# Copy the ClinSV binaries to  /app
ADD ClinSV_x86_64_v1.0.tar.gz /app

#Copy the samtools reinstalltion script to /app
ADD reinstallSamtools.sh /app

# Install root
RUN yum -y --enablerepo=extras install epel-release
RUN yum -y install gcc-c++ bzip2 git vim libpng libjpeg      python2-devel libSM libX11 libXext libXpm libXft gsl-devel make cmake3     && yum -y clean all
RUN yum -y groupinstall 'development tools'

# Installing dependencies required for samtools reinstallation
RUN yum -y install autoconf automake make gcc perl-Data-Dumper zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel ncurses-devel 

# Install wget to allow for the installation of samtools + htslib + bcftools
RUN yum -y install wget

# Install terminal text editor for debugging containers
RUN yum -y install nano

# Run the samtools reinstallation script
RUN bash reinstallSamtools.sh

# Set ROOT environment
ENV ROOTSYS         "/app/clinsv/root"
ENV PATH            "$ROOTSYS/bin:/bin/bin:$PATH"
ENV LD_LIBRARY_PATH "$ROOTSYS/lib:$LD_LIBRARY_PATH"
ENV PYTHONPATH      "$ROOTSYS/lib:$PYTHONPATH"

# Install any needed packages specified in requirements.txt
RUN yum install -y libpng-devel
RUN yum install -y mysql

# for QC reprot
RUN yum install -y libgfortran
RUN yum install -y texlive

# Define environment variable
ENV PATH="/app/clinsv/bin:${PATH}"

# Run app.py when the container launches
ENTRYPOINT ["/app/clinsv/bin/perl","/app/clinsv/bin/clinsv"]
FROM ubuntu:20.04
RUN apt-get update -y
# default is python 3.8.10
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
# Make Python3 the default
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1
RUN update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

RUN apt-get install -y libz-dev libbz2-dev liblzma-dev
# install HTSLIB - /usr/include/htslib/
RUN apt-get install libseqlib-dev
RUN apt-get install samtools

RUN pip3 install pysam==0.15.4 pyBigWig==0.3.17 pandas==1.2.4 
RUN pip3 install scipy==1.7.1 matplotlib==3.4.1 pyyaml==5.3.1 argparse==1.1
RUN pip3 install pybedtools==0.8.0 
RUN pip3 install statsmodels

WORKDIR ${BUILD_DIR}
SHELL ["/bin/bash", "-c"]

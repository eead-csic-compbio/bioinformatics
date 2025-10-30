
## build instructions:
# 
# THISVERSION=`date "+%Y%m%d"` 
# docker build --tag csicunam/bioinformatics_iamz:$THISVERSION --tag csicunam/bioinformatics_iamz:latest --force-rm --compress .
#
# docker login 
# docker push csicunam/bioinformatics_iamz:${THISVERSION}
# docker push csicunam/bioinformatics_iamz:latest
#
# mkdir $HOME/vep_data
# docker run -t -i -v $HOME/vep_data:/data -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY csicunam/bioinformatics_iamz:latest

## Clone of ensembl-vep image, locale set to en_US.UTF-8 ,
## see https://github.com/Ensembl/ensembl-vep/blob/release/115/docker/Dockerfile
FROM ensemblorg/ensembl-vep:latest 

LABEL maintainer="bcontreras@eead.csic.es"

## Prevent to open interactive dialogs during the installation process
ENV DEBIAN_FRONTEND=noninteractive

## Install required packages
USER root
RUN apt-get update && apt install --no-install-recommends -y \
	make\
	wget\
	git\
	bwa\
	samtools\
        bedtools\
	picard-tools\
	fastqc\
	freebayes\
	bcftools\
	tabix\
	beagle\
	kallisto\
        trimmomatic

## Clone repo with course materials, including test_data
WORKDIR /home/vep
RUN git clone https://github.com/eead-csic-compbio/bioinformatics.git
RUN cd bioinformatics/test_data/transcripts && ./HOWTO.txt
RUN ln -s bioinformatics/test_data .

## Copy course datafiles, data is a repo folder
RUN mkdir -p /home/vep/variant_data/
COPY data /home/vep/variant_data/

## download and install get_homologues
RUN git clone https://github.com/eead-csic-compbio/get_homologues.git
RUN cd get_homologues && make install_docker 
ENV GETHOMS /home/vep/get_homologues

# add version name to image
ARG version
LABEL version=$version
RUN echo $version

# prepare user env
USER vep
WORKDIR /home/vep/

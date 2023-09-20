
## build instructions:
# 
# THISVERSION=`date "+%Y%m%d"` 
# docker build --tag csicunam/bioinformatics_iamz:$THISVERSION --tag csicunam/bioinformatics_iamz:latest --force-rm --compress .
#
# docker login 
# docker push csicunam/bioinformatics_iamz:${THISVERSION}
# docker push csicunam/bioinformatics_iamz:latest


## Create a clone of ensembl-vep image, locale set to en_US.UTF-8 ,
## see https://github.com/Ensembl/ensembl-vep/blob/release/110/docker/Dockerfile,
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
	igv\
	freebayes\
	bcftools\
	tabix\
	beagle


# Clone repo with course materials, including test_data
WORKDIR /home/vep
RUN git clone https://github.com/eead-csic-compbio/bioinformatics.git

## Copy course datafiles
RUN mkdir -p /home/vep/variant_data/
COPY data /home/vep/variant_data/

# add version name to image
ARG version
LABEL version=$version
RUN echo $version

# prepare user env
USER vep
WORKDIR /home/vep/


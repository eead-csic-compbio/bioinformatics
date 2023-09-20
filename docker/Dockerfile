## Create a clone of VEP image, locale set to en_US.UTF-8
## see https://github.com/Ensembl/ensembl-vep/blob/release/110/docker/Dockerfile
FROM ensemblorg/ensembl-vep:latest 

LABEL maintainer="bcontreras@eead.csic.es"

## Prevent to open interactive dialogs during the installation process
ENV DEBIAN_FRONTEND=noninteractive

## Install required packages
USER root
RUN apt-get update && apt install --no-install-recommends -y \
	curl\
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
	make

## Copy course data



# add version name to image
ARG version
LABEL version=$version
RUN echo $version

# prepare user env
USER vep
WORKDIR /data



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

# Install Miniconda on x86 or ARM platforms
RUN arch=$(uname -m) && \
    if [ "$arch" = "x86_64" ]; then \
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"; \
    elif [ "$arch" = "aarch64" ]; then \
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"; \
    else \
    echo "Unsupported architecture: $arch"; \
    exit 1; \
    fi && \
    wget $MINICONDA_URL -O miniconda.sh && \
    mkdir -p /root/.conda && \
    bash miniconda.sh -b -p /root/miniconda3 && \
    rm -f miniconda.sh

ENV PATH="/root/miniconda3/bin::${PATH}"
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

## Clone repo with course materials, including test_data
WORKDIR /home/vep
RUN git clone https://github.com/eead-csic-compbio/bioinformatics.git
RUN cd bioinformatics/test_data/transcripts && ./HOWTO.txt
RUN ln -s bioinformatics/test_data .

## Copy course datafiles, data is a repo folder
RUN mkdir -p /home/vep/variant_data/
COPY data /home/vep/variant_data/

# stuff installed from bioconda
COPY environment.yml environment.yml
RUN conda env create -f environment.yml
RUN echo "conda activate bioinfo" >> ~/.bashrc
ENV PATH /opt/conda/envs/bioinfo/bin:$PATH
ENV CONDA_DEFAULT_ENV bioinfo

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

ENV PATH="${GETHOMS}:${GETHOMS}/bin:${PATH}"
SHELL ["conda", "run", "--no-capture-output", "-n", "bioinfo", "/bin/bash", "-c"]

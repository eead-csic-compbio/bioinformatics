# Open Linux (WSL)


# Create the folder "transcriptome_class"

mkdir transcriptome_class

# Extend permissions

chmod 777 transcriptome_class
-----------------------------------------------------------------

# Docker environment

Web: https://github.com/eead-csic-compbio/bioinformatics/tree/main/docker

# Run (in Ubuntu):

docker pull csicunam/bioinformatics_iamz:latest

# bind local "transcriptome_class" folder to container folder "/data"

docker run -it -v $HOME/transcriptome_class:/data -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY csicunam/bioinformatics_iamz:latest


## From Windows, the "transcriptome_class" folder is located --> \\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class

-----------------------------------------------------------------
# Move into /data --> cd /data


# Create symbolic link:

ln -s /home/vep/test_data/transcripts/brachy/ brachy

# Create report folders:

mkdir report report_filtered

-----------------------------------------------------------------

# Run FastQC in raw reads

fastqc brachy/01_RNAseq_raw/*.gz -o report

# Check reports in Windows:
\\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\report

-----------------------------------------------------------------

# Run Trimmomatic according to the reports

# Create output folder in /data (Linux WSL)

mkdir 02_RNAseq_filtered

## Trimming occurs in the order which the steps are specified on the command line.
## It is recommended in most cases that adapter clipping, if required, is done as early as possible.

## Example: Filtering of ABR3_W_03

TrimmomaticSE -threads 4 -phred33 \
brachy/01_RNAseq_raw/ABR3_W_03.fastq.gz \
02_RNAseq_filtered/ABR3_W_03_trimmed.fastq.gz \
ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## Example: Filtering of ABR3_D_03

TrimmomaticSE -threads 4 -phred33 \
brachy/01_RNAseq_raw/ABR3_D_03.fastq.gz \
02_RNAseq_filtered/ABR3_D_03_trimmed.fastq.gz \
ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

-----------------------------------------------------------------

# Run FastQC in filtered/trimmed reads

fastqc 02_RNAseq_filtered/*.gz -o report_filtered

# Check reports in Windows:
\\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\report_filtered


-----------------------------------------------------------------

# Kallisto

## Index reference (kallisto index):

kallisto index -i Bdistachyon_314_v3.1.transcript.TagSeq500b.idx \
brachy/00_Reference_transcriptome/Bdistachyon_314_v3.1.transcript.TagSeq500b.fa.gz


## Quantification for one sample (kallisto quant)

mkdir 03_quant && mkdir 03_quant/ABR3_W_03

kallisto quant -i Bdistachyon_314_v3.1.transcript.TagSeq500b.idx \
-o 03_quant/ABR3_W_03 \
-b 100 --single -l 100 -s 20 <(zcat 02_RNAseq_filtered/ABR3_W_03_trimmed.fastq.gz)


# Check abundances: more 03_quant/ABR3_W_03/abundance.tsv

-----------------------------------------------------------------
-----------------------------------------------------------------


# For Sleuth and WGCNA:

# Into /data (Linux WSL)

cp -r /home/vep/test_data/transcripts/brachy/05_WGCNA/ .

cp -r /home/vep/test_data/transcripts/brachy/04_Sleuth/ .

# Extent permissions
chmod -R 777 04_Sleuth 05_WGCNA

# Copy Rstudio scripts:

Copy Rstudio script "Practical_Sleuth.rmd" in \\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\04_Sleuth

Copy Rstudio script "Practical_WGCNA_W_dataset" in \\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\05_WGCNA

Copy Rstudio script "Practical_WGCNA_D_dataset" in \\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\05_WGCNA


#########################
######### SLEUTH ########
#########################

# We will work in \wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\04_Sleuth  *Change rsancho by your username (iamz)

# Open Practical_Sleuth.rmd



#########################
######### WGCNA ########
#########################

# We will work in \\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\05_WGCNA  *Change rsancho by your username (iamz)

# Open Practical_WGCNA.rmd
































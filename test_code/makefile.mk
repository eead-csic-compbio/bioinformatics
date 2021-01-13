###########################################################################################
## This makefile is to retrieve 5 random gene_ids from a list of arabidopsis genes, 
## then will assign the corresponding fasta sequence to each gene. 
## Gene list will be downloaded from a github and fasta sequence from EnsemblPlants
###########################################################################################

## Define the direcotries
DATA=~/demo_project/sequences
INPUT=~/demo_project/input
RESLTS=~/demo_project/results
FIGURES=~${RESLT}/figures
ARAB_IDS=${INPUT}/arabidopsis_genes.txt
ARAB_SEQ=${DATA}/Arabidopsis_thaliana.TAIR10.cdna.all.fa


## Start with the targets

# Create the directories
# this target is called 'dir'
# all command lines or targets need to be indented with a tab
dir:
	@echo "Creating sequences and results directories"
	mkdir -p ${DATA}
	mkdir -p ${RESLTS}
	mkdir -p ${INPUT}

# Download the data
# Note: this uses wget, which might not be installed in MacOS systems
# wget can often be substituted with curl with minimal changes
get-data:
	@echo "Getting the gene list from github"
	wget https://raw.githubusercontent.com/eead-csic-compbio/bioinformatics/main/test_data/arabidopsis_genes.txt
	@echo "getting the Arabidopsis_thaliana cdna sequences from ensemblPlants"
	wget ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
	@echo "Unzipping the cdna file"
	gunzip Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
	@echo
	@echo "Moving the data to their corresponging folders"
	mv arabidopsis_genes.txt ${INPUT}
	mv Arabidopsis_thaliana.TAIR10.cdna.all.fa ${DATA}

# Get random gene ids
get-ids:
	@echo
	@echo "Generating 5 random genes from the gene list"
	shuf -n 5 ${ARAB_IDS} -o ${RESLTS}/random_genes.txt
	@echo "${RESLTS}/random_genes.txt"

# Assign fasta
get-fasta:
	@echo
	@echo "Assign fasta sequence to each gene"
	cat ${RESLTS}/random_genes.txt | xargs -n 1 samtools faidx ${ARAB_SEQ} > ${RESLTS}/random-sequences.fna
	@echo "${RESLTS}/random-sequences.fna"

# Run all targets at one go
all: dir get-data get-ids get-fasta
	@echo "tasks successfully done"


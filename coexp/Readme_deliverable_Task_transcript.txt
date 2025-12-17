EXERCISE:

CONSTRUCT A GENE CO-EXPRESSION NETWORK (AND DEFINE THE COEXPRESSION MODULES) USING THE RNA-SEQ DATA FROM THE WATERED (W) SAMPLES AND ANSWER THE FOLLOWING QUESTIONS:

(1) How many samples and isoforms are included in the dataset "TPM_counts_Drought_W_dataset.csv"?

(2) How many samples are discarded after outlier analysis?

(3) What power value have you set as appropriate for calculating adjacency?

(4) How many co-expression modules are established before and how many after the module merging process?

(5) What is the hub isoform (or hub gene) of the cyan module?

(6) According to the module-trait association heat map, which module has the highest positive correlation with the "blwgrd (below ground biomass)" trait?


TIPS:

- All the instructions, step by step, are included in the file "Readme_deliverable_Task_transcript" available on Github.
- To answer the questions, run each chunk one by one and interpret the result obtained.



###############################################################################################
###############################################################################################


IMPORTANT:

To complete the deliverable exercise, the inputs you need are:

- TPM_counts_Drought_W_dataset.csv
- TRAITS_W.csv

And the R script you will use is the one you can download from:

https://github.com/eead-csic-compbio/bioinformatics/blob/main/coexp/Practical_WGCNA_W_dataset.Rmd

----------------------------------------------------------------------------------------

NOTE: If you did not download the input files during the class on 12/16/2025, you can download them by following these steps:

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

## From Windows, the "transcriptome_class" folder is located in:

\\wsl.localhost\Ubuntu-22.04\home\USERNAME\transcriptome_class

(Replace USERNAME with your username. The Ubuntu version shown (22.04) may also change depending on the version you installed).

-----------------------------------------------------------------
# Move into /data

cd /data

# Create symbolic link:

ln -s /home/vep/test_data/transcripts/brachy/ brachy

-----------------------------------------------------------------

# For WGCNA:

# Into /data (Linux WSL)

cp -r /home/vep/test_data/transcripts/brachy/05_WGCNA/ .

# Extent permissions

chmod -R 777 04_Sleuth 05_WGCNA

# Copy Rstudio scripts:


Copy Rstudio script "Practical_WGCNA_W_dataset" in:

\\wsl.localhost\Ubuntu-22.04\home\USERNAME\transcriptome_class\05_WGCNA


#########################
######### WGCNA ########
#########################

# We will work in \\wsl.localhost\Ubuntu-22.04\home\USERNAME\transcriptome_class\05_WGCNA 

# Inside the 05_WGCNA folder you should have the following files:

- TPM_counts_Drought_W_dataset.csv
- TRAITS_W.csv
- Practical_WGCNA_W_dataset.Rmd


###############################################################################################

You must submit the solution (and files) to the exercise following these instructions:

- Create a folder called ‘transcripts/’ in the same GitHub repository as session 1 and save all the relevant files created during the execution of the task there.

- Include a text file named "Answers_transcript_YOURNAME.txt" that displays the questions (listed again below) and the answers:

(1) How many samples and isoforms are included in the dataset "TPM_counts_Drought_W_dataset.csv"?

(2) How many samples are discarded after outlier analysis?

(3) What power value have you set as appropriate for calculating adjacency?

(4) How many co-expression modules are established before and how many after the module merging process?

(5) What is the hub isoform (or hub gene) of the cyan module?

(6) According to the module-trait association heat map, which module has the highest positive correlation with the "blwgrd (below ground biomass)" trait?
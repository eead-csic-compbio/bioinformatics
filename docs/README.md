# Bioinformatics resources and best-practices for plant breeders

This material is maintained by Najla Ksouri, Chesco Montardit, Ernesto Igartua, Bruno Contreras Moreira 
and the Ensembl outreach team

##  Summary

Here we review some bioinformatics resources and databases which can be useful in plant breeding and genomics. 
We will use both standalone and Web-based tools and will also review reproducible analysis practices and software benchmarks.
Test data used in sessions 1-4 can be obtained from <https://github.com/eead-csic-compbio/bioinformatics>.

## Docker image

A [Docker image](https://hub.docker.com/r/csicunam/bioinformatics_iamz) 
is available with most of the software used in the sessions, excluding R,
which we expect to be installed elsewhere.
After installing Docker, it can be run as follows:

    docker pull csicunam/bioinformatics_iamz

    # persistent folder for results files
    mkdir $HOME/vep_data 
    chmod a+w $HOME/vep_data

    docker run -t -i -v $HOME/vep_data:/data -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY csicunam/bioinformatics_iamz:latest



## Contents

|session|title|required time|URL|
|-------|-----|-------------|---|
|1|Annotation of coding sequences|4h|[session 1](./session1.html)|
|2|Analysis of non-coding sequences|4h|[session 2](./session2.html)|
|3|Reproducible analysis practices|2h|[session 3](./session3.html)|
|4|Benchmarks|2h|[session 4](./session4.html)|
|5|Mapping, variant calling & effect prediction||session 5|
|6|Genotyping|3h|[session 6](./session6.html)|
|7|Genome-Wide Association Analysis|2h|[session 7](./session7.html)| 



We post regularly about these and related bioinformatics topics at the [#!/perl/bioinfo](https://bioinfoperl.blogspot.com) blog, mostly in Spanish.

Check also this course to learn how to [script in Linux](https://github.com/eead-csic-compbio/scripting_linux_shell).

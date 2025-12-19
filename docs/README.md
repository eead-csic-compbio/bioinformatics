# Bioinformatics resources and best-practices for plant breeders

This material is maintained by Najla Ksouri, Chesco Montardit, Rubén Sancho, Ernesto Igartua, Ricardo Ramírez González, Bruno Contreras Moreira, MGnify and the Ensembl outreach team

##  Summary

Here we review some bioinformatics resources and databases which can be useful in plant breeding and genomics. 
We will use both standalone and Web-based tools and will also review reproducible analysis practices and software benchmarks.
Test data used in sessions 1-5 can be obtained from <https://github.com/eead-csic-compbio/bioinformatics>.

## Docker image

A [Docker image](https://hub.docker.com/r/csicunam/bioinformatics_iamz) 
is available with most of the software used in the sessions, excluding R,
which we expect to be installed elsewhere.
After installing Docker, it can be run as follows, note that you might require **sudo**:

    docker pull csicunam/bioinformatics_iamz

    # Set persistent folder for results files to avoid data loss when docker is turned off
    # For instance, for the VEP class you could create one named 'vep_class'
    mkdir $HOME/vep_class 
    chmod a+w $HOME/vep_class
    
    # launch docker binding the persistent folder to internal folder (/data)
    sudo docker run -t -i -v $HOME/vep_class:/data -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY csicunam/bioinformatics_iamz:latest



## Contents

|session|title|required time|URL|
|-------|-----|-------------|---|
|1|Annotation of coding sequences|4h|[session 1](./session1.html)|
|2|Analysis of non-coding sequences|4h|[session 2](./session2.html)|
|3|Analysis of transcriptomes|4h|[session 3](./session3.html)|
|4|Benchmarks|2h|[session 4](./session4.html)|
|5|Mapping, variant calling & effect prediction|6h|[session 5](./session5.html) , [session 5a](./session5a.html)|
|6|Genotyping|3h|[session 6](./session6.html)|
|7|Genome-Wide Association Analysis|2h|[session 7](./session7.html)|
|8|metagenomics and amplicon intro|4h|[session 8](https://plants-breeding.mgnify.org)
|R|Reproducible analysis practices|2h|[session R](./sessionR.html)|


## Exercises and report

Each session contains exercises (Exe1, Exe2, ...) which you can solve and document in a report.
When we teach this material, we ask students to create a GitHub repository, with a dedicated folder per session,
explaining how each exercise was solved adding code, comments, even figures, and literature references if needed.
Moreover, any AI resources used in your work, such as ChatGPT, should be properly cited and the relevant queries
**included in the report**. 
[Markdown](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) 
is learned quickly and is a good format for reports, as it is supported by GitHub.

The idea is that students **log their work as they go**, as opposed to uploading a set of solutions on the last day,
so that the thought process and progress is visible, as well as challenges. 
The resulting repo can be evaluated by teachers but also serves as a portfolio of skills and knowledge for potential employers.

[session R](https://eead-csic-compbio.github.io/bioinformatics/sessionR.html#225_Use_version_control_systems) 
provides examples on setting up a GitHub repository and using Git for version control, and also on
the slightly advanced
[Rmarkdown](https://eead-csic-compbio.github.io/bioinformatics/sessionR.html#224_Turn_your_scripts_into_reproducible_reports) format.

## More resources

We post regularly about these and related bioinformatics topics at the [#!/perl/bioinfo](https://bioinfoperl.blogspot.com) blog, mostly in Spanish.

Check also this course to learn how to [script in Linux](https://github.com/eead-csic-compbio/scripting_linux_shell).

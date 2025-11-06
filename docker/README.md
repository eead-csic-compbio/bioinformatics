
# Docker container for the bioinformatics module at CIHEAM Zaragoza

This [Dockerfile](./Dockerfile) is used to build a standalone container for the bioinformatics module at CIHEAM Zaragoza,
available at <https://hub.docker.com/r/csicunam/bioinformatics_iamz>, which can be installed and run as follows:

    docker pull csicunam/bioinformatics_iamz:latest

    # bind local vep folder to container folder data/
    docker run -it -v $HOME/vep_data:/data -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY csicunam/bioinformatics_iamz:latest
 
    # or any other local data folder to container folder data/
    docker run -it -v /path/to/data:/data -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY csicunam/bioinformatics_iamz:latest

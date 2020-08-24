---
title: Setup
---

# Softwares

## Option 1: using a Docker image

The preferred option to install all softwares and packages is to use a tailor-made Docker image. This image will also contain the necessary datasets ([they are described and available here](https://zenodo.org/record/3755641) if you want to look at them independently of the Docker image).   
See [this nice introduction to Docker here](https://aws.amazon.com/docker/).   

This image is based on the [rocker verse Docker image](https://hub.docker.com/r/rocker/verse) image with extra R libraries. The latest image can be [found at the Science Park Study Group DockerHub](https://hub.docker.com/repository/docker/scienceparkstudygroup/master-gls:microbiome-latest).


> ## Before you start
>
> Before the training, please make sure you have done the following: 
>
> 1. First, install [Docker desktop](https://www.docker.com/products/docker-desktop) for your operating system.  
> 2. If needed, install Shell Bash: [follow these instructions](http://swcarpentry.github.io/shell-novice/setup.html).
> 3. Open a new Shell Bash window and navigate to a folder that will be your workspace. For instance, you could create a folder named `r-tutorial/` on your Desktop and move inside with the Shell using `cd ~/Desktop/r-tutorial/`. 
> 4. In a Shell Bash window, type the following command: `docker run --rm --name rstudio_instance -v $PWD:/home/rstudio/ -e PASSWORD=mypwd -p 8787:8787 scienceparkstudygroup/master-gls:microbiome-latest`. This will download a Docker image for the course, create and run a container where RStudio will be running.  
> 4. Navigate to [http://localhost:8787](http://localhost:8787) in your web browser. You should have an RStudio session running. Type `rstudio` as the user name and `mypwd` as your password. 
> 5. To quit, close the web browser window where RStudio is running and exit the Shell too. 
{: .prereq}

__Docker command-line explanations:__  
- The `--rm` removes the container when it has been run. No need to store it into your computer after use.      
- The `--name` gives a name to the running container for easy retrieval.  
- The `-p 8787:8787` follow the format `-p host_port:container_port`. Therefore the port 8787 inside the container will be exposed to the outside port on the host machine. That way, the running instance of RStudio can be access through the <IP address>:port format.

> ## Important note
>
> You can save files to your disk when working inside the Docker-powered R session. You need to save them as you would normally. The files (e.g. `my_plot.png`) will be where you were working (the directory from which you launched the Docker container). If you are in `~/Desktop/r-tutorial/` then your saved files should be there. 
{: .callout}




## Option 2: manual installation
This is the second way to install softwares and packages. It _should_ work but there is no guarantee that it _will_ work since R and packages versions on your machine might be different from the software and package versions used in this lesson. Thus, the preferred way is still to use the Docker image (option 1).  



> ## Before you start: install the softwares
>
> Before the training, please make sure you have done the following: 
>
> 1. Download and install **up-to-date versions** of:
>    - R: [https://cloud.r-project.org](https://cloud.r-project.org)
>    - RStudio: [http://www.rstudio.com/download](http://www.rstudio.com/download) 
>    - Git: [https://git-scm.com/downloads](https://git-scm.com/downloads) 
>
> 2. Within R/RStudio, install the required R packages. To do so, open R and use the `install.packages()` function. 
>
> ~~~
> # these packages come from the "Comprehensive R Archive Network (CRAN)"
> # The most recent version of the package will be installed 
> install.packages("tidyverse")
> install.packages("vegan")
> install.packages("ade4")
> install.packages("patchwork")
> install.packages("agricolae")
> install.packages("FSA")  
> install.packages("rcompanion")  
> install.packages("BiocManager")

>
> # this package comes from Bioconductor 	 
> BiocManager::install('phyloseq')
> ~~~
> {: .language-r}
{: .prereq}


<br>

Once you have the softwares installed, you'll need the files as well. Follow the instructions below.

> ## Before you start: download input files
> All datasets can be found [at the Zenodo archive](https://zenodo.org/record/3755641) and consist of:
> 1. An OTU table: `data_loue_16S_nonnorm.txt`  
> 2. A sample metadata table: `data_loue_16S_nonnorm_grp.txt`  
> 3. A taxonomy table: `data_loue_16S_nonnorm_taxo.txt`  
{: .prereq}

<br>
<br>

{% include links.md %}

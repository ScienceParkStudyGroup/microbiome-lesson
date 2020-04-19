---
title: "Data exploration and properties"
teaching: 60
exercises: 0
questions:
- Which packages do I need for the analyses?  
- How can I import my data in R?  
- What the data sets look like?  
- What are the properties for the occurrence table?  
objectives:
- Install and load the required packages in R  
- Import the data sets  
- Create a phyloseq object  
- Explore the data sets  
- Assess sparsity of the occurrence table  
- Assess the sequencing depth  
keypoints:
- ""
---

## Table of Contents
<<<<<<< HEAD:_episodes/02-data-exploration-and-properties.md
1. [Install and load the required packages  ](#install-and-load-the-required-packages)
2. [Import the data in R  ](#import-the-data-in-r)
3. [Create a phyloseq object  ](#create-a-phyloseq-object)
4. [Global exploration of the data sets  ](#global-exploration-of-the-data-sets)
5. [OTU table (based on the raw data: *data_otu* table or *OTU* from *data_phylo* object)](#otu-table-based-on-the-raw-data-data_otu-table-or-otu-from-data_phylo-object)
6. [Sample metadata table ](#sample-metadata-table)
7. [OTU data properties: sparsity ](#otu-data-properties-sparsity)

=======
1. [Install and load the required packages](#install-and-load-the-required-packages)
2. [Import the data in R](#import-the-data-in-r)
3. [Create a phyloseq object](#create-a-phyloseq-object)
4. [Global exploration of the data sets](#global-exploration-of-the-data-sets)
5. [OTU data properties: sparsity](#otu-data-properties-sparsity)
>>>>>>> 996fb4deaa077dea642eb57ca247b2939e03e3b4:_episodes/02.data-exploration-and-properties.md

## 1. Install and load the required packages  

~~~
# install.packages("Vegan")     # only if you need to install the package
# install.packages("phyloseq")
# install.packages("tidyverse")
# install.packages("patchwork")
# install.packages("agricolae")
# install.packages("FSA")
# install.packages("rcompanion")

library(vegan)
library(phyloseq)
library(tidyverse)
library(patchwork)
library(agricolae)
library(FSA)
library(rcompanion)
~~~
{: .language-r}


## 2. Import the data in R  

~~~
data_otu <- read.table("data_loue_16S_nonnorm.txt", header = TRUE)

data_grp <- read.table("data_loue_16S_nonnorm_grp.txt", header = TRUE)

data_taxo <- read.table("data_loue_16S_nonnorm_taxo.txt", header = TRUE)
~~~
{: .language-r}
  
The OTU table (`data_otu` or `data_loue_16S_nonnorm.txt`) is the occurrence table that contains counts corresponding to the number of times each OTU (Operational Taxonomic Unit, here bacteria) is observed in each sample.  

The sample metadata table (`data_grp` or `data_loue_16S_nonnorm_grp.txt`) represents the metadata information on the different samples, in this case, where (*site*) and when (*month*) the samples were harvested.    

The taxonomy table (`data_taxo` or `data_loue_16S_nonnorm_taxo.txt`) represents the assignment for the different OTU at domain, phylum, class, order, family, genus and species levels. Usually, microbial ecologists use phylum level (or class level for Proteobacteria) to describe the global bacterial communities. Then, they usually use OTU, ASV (Amplicon Sequence Variant) or genus level to check which bacteria are differentially occuring in the different treatments.  

## 3. Create a phyloseq object  
  
Microbial ecologists usually use the vegan and phyloseq packages to analyse the occurrence table. Such as `biom` format, `phyloseq` format support encapsulation of core study data (occurrence table data and sample/observation metadata) in a single file. From the three tables `data_otu`, `data_grp` and `data_taxo` we will create a phyloseq object `data_phylo`.  
  
~~~
OTU = otu_table(as.matrix(data_otu), taxa_are_rows = FALSE) # create the occurrence table object in phyloseq format

SAM = sample_data(data_grp, errorIfNULL = TRUE) # create the sample metadata object in phyloseq format

TAX = tax_table(as.matrix(data_taxo)) # create the observation metadata object (OTU taxonomy) in phyloseq format

data_phylo <- phyloseq(OTU, TAX, SAM) # create the phyloseq object including occurrence table data and sample/observation metadata

data_phylo # print information about the phyloseq object
~~~~
{: .language-r}

## 4. Global exploration of the data sets  
  
### 4.1. OTU table  
Based on raw data: `data_otu` table or `OTU` from `data_phylo` object  
The OTU table is the occurrence table that contains counts.  
  
Check how the data look like showing the first 5 rows and the first 6 columns.  
~~~
data_otu[1:5, 1:6] # as we have here a high number of variables it is better to not use head function for a better readability
~~~
{: .language-r}


The OTU table has samples in rows and variables (OTU) in columns.  

> ## Remark
> If you want to use phyloseq, you can do the same running `otu_table(data_phylo)[1:5, 1:6]` 
{: .callout}
  
  
Check how many samples and variables are in the OTU table?  
  
~~~
nb_samples <- dim(data_otu)[1] # nb of raws, here samples
nb_samples

nb_var <- dim(data_otu)[2] # number of columns, here variables (OTU)
nb_var
~~~
{: .language-r}


> ## Remark
> If you want to use phyloseq, you can use `nsamples()` and `ntaxa()` functions, such as `nsamples(data_phylo)` and `ntaxa(data_phylo)`.  
{: .callout}

### 4.2. Sample metadata table  
Based on `data_grp` table or `SAM` from `data_phylo` object  
This table represents the metadata information on the different samples.  

> ## Exercise
>
> How does the data look like if you show the first 6 rows using head function?   
> 
> > ## Solution
> > `# The sample metadata table has samples in rows and factors in columns.`     
> > `head(data_grp)`
> {: .solution}
{: .challenge}  

> ## Remark
>If you want to use phyloseq, you can do the same running `head(sample_data(data_phylo))`.  
{: .callout}

> ## Exercise
>
> How many samples and factors are in the metadata table?   
> 
> > ## Solution
> > `dim(data_grp)[1] # number of samples.`     
> > `nb_factors <- dim(data_grp)[2] # number of factors`  
> > `nb_factors`
> {: .solution}
{: .challenge}  


  
  
Check how many samples per treatment do we have.  
~~~
summary(data_grp)
~~~
{: .language-r}
  
There are nine samples harvested or replicates (ignoring the month of harvest) for each of the two sites: Cleron (located at the upstream area of the river) and Parcey (located at the downstream area of the river). There are six replicates per month (ignoring the site of harvest). There are three replicates considering both site and month of harvest.  
  
> ## Remark
> If you want to use phyloseq, you can do the same running dim(sample_data(data_phylo))[1], dim(sample_data(data_phylo))[2] and summary(sample_data(data_phylo)$site)  
{: .callout}
 
 
> ## Discussion point
>
> This is a multivariate data set with an underlying experimental design. What would be your first idea to analyze this data?  
{: .discussion}

In the following part of this tutorial, we will discuss some problematic issues of these data that makes it necessary to use a specific strategy to deal with these data.  
  
  

## 5. OTU data properties: sparsity 
Based on the raw data: `data_otu` table  
  
Microbiome data sets are usually sparse.  
  
### 5.1. Number of zeros and percentage of zeros in the OTU table  
  
~~~
sum(data_otu == 0)
sum(data_otu == 0) / (nb_var * nb_samples) * 100
~~~
{: .language-r}


> ## Question
> 1. According to your knowledge, what could be the meaning of a zero value in the OTU/ASV table? 
> 2. Can we consider it as a true zero? 
> 3. Do you think that it is important to filter the data and why?  
> 
> > ## Solution
> > 1. A zero value in the OTU table means that we could not sequence any read for a specific OTU/ASV in a given sample. This means that there is not such OTU/ASV in the given sample or that the sequencing depth was not enough to find it.  
> > 2. So, we usually don't know if it is a true zero or not.  
> > 3. First, some statistical approcahes cannot handle highly sparse data set (for example Principal Component Analysis).  
> > Then, there is often problem when we have to compare rare OTU/ASV. For example, if for one sample you cannot find any 
> > count for a specific OTU but for another sample you can find 1 count for this OTU, can you say that there is more of this 
> > OTU in the second sample than in the first one?  
> > Therefore, in the context of microbiota data analysis, zeros values can be difficult to analyse and interpret and data 
> > should often be filtered.  
> {: .solution} 
{: .challenge}


> ## Remark 
> Here, the percentage of zeros is relatively high. In order to be able to apply specific statistical approaches later, we should think about filtering the OTU data.  
{: .callout}   
   
### 5.2. Counts frequency  
  
Visualize the count frequency in the OTU table using a histogram.  
~~~
hist(as.matrix(data_otu), 
     max(data_otu), 
     right = FALSE, 
     las = 1, 
     xlab = "Occurrence value", 
     ylab = "Frequency", 
     main = "Occurrence frequency")
~~~
{: .language-r}
  
> ## Question 
> How do you interpret this histogram?  
> > ## Solution 
> > This plot represents the occurence frequency. The x axis represents the occurence value (for the all OTU table, so for 
> > each sample and OTU) and y the frequency.  
> > You can see that for x = 0, y = 68633. So, there are 68633 zeros in the OTU table.  
> > You can also see that the histogram is rapidely decreasing when x increases. That means that there is a lot of low counts > > occurences (such as 0, 1 or 2 counts per OTU and per sample) and few abundant OTU.  
> {: .solution}
{: .challenge}

> ## Remark 
> One option for data fitering could be to remove OTU that have a number of counts for all the samples lower than a define value.   
{: .callout}   
  
  
### 5.3. Minimum of counts per OTU for all the samples  
  
> ## Question
> What is the minimum number of counts per OTU for all the samples in this data set?
> > ## Solution
> > `min(colSums(data_otu))`
> {: .solution}
{: .challenge}

   
> ## Remark 
> As illustrated in this tutorial, after the first bioinformatic analysis that generate the OTU table, microbial ecologists usually remove singletons (OTU for which only one read has been identified in the full OTU table, including all the samples).  
> Microbial ecologists consider that reads which are present only once should be treated as artifacts and removed. As a 
> reminder, chimeric amplification products have been already removed in the bioinformatic analysis before generating the OTU 
> table. Chimeras are hybrid products between multiple parent sequences that can be falsely interpreted as novel organisms, thus inflating apparent diversity.  
> thus inflating apparent diversity.  
{: .callout} 
 
  
### 5.4. Non-zero values per OTU  
  
In order to check how the different OTU/ASV are shared bewteen samples, plot the number of non-zero values for each OTU  
~~~
non_zero <- 0*1:nb_var

for (i in 1:nb_var){
  non_zero[i]<-sum(data_otu[,i] != 0)
  }

plot(non_zero, xlab = "OTU", ylab = "Frequency", main = "Number of non zero values", las = 1)
~~~
{: .language-r} 
  
> ## Question 
> How do you interpret this plot?  
> > ## Solution:
> > This plot represents the number of non-zero values. The x axis represents the different OTU (we have here 5248 OTU) and 
> > the y axis represents the frequency (from 1 to 18, as we have 18 samples).  
> > You can see that there is a lot of dots when y < 5, which means that the major part of the non-zeros values are found in 
> > less than 5 samples. We can also say that there are few non-zero counts that are shared in more than 5 samples.  
> > The major part of the OTU are thus specific to only a few number of samples.  
> {: .solution}
{: .challenge}
  
> ## Remark
> Remember that, here, the number of replicates per treatment is 3 per site and per harvesting date, 6 per harvesting date, 
> and 9 per site.  
> We can think about different ways to filter the OTU data. For example, here, we can think about removing rare OTU that are 
> not shared with less than 2 or 3 samples. Try to think about how you can filter your OTU data.  
{: .callout}
  
  
  
  



  
  

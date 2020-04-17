---
title: "Beta-diversity"
teaching: 60
exercises: 
questions:
- "Question one."
objectives:
- "Objective one."
keypoints:
- "Keypoint one."

---

## Data filtering  
  
In the remainder of the data analyses, we want to focus on the differences between samples, bacterial community composition and beta-diversity. 

For this it is necessary to prepare the data in a way that improves the comparability of the samples. As microbiome data sets are usually sparse, it is important to filter the data set. Indeed, data filtering aims to remove low quality or uninformative variables (OTU/ASV) to improve downstream statistical analysis.  

For example, we can filter variables (OTU/ASV) with very low number of counts in only a few samples, which are likely due to sequencing errors or low-level contaminations. In this example, we will keep OTU that have at least 2 counts in at least 11% of the samples. Indeed, here we have 18 samples and 3 replicates per treatment and we want, for example, to have at least 2 counts in at least two samples (2/18 = 0.111).  
  
> ## Remark
> First, data filtering should be applied on the raw data. Then, we are using here the phyloseq object: `data_phylo` 
> because we are using a plyloseq function to filter the raw data.  
{: .callout}

~~~
# filter the OTU data using filter_taxa function included in phyloseq package
data_phylo_filt = filter_taxa(data_phylo, function(x) sum(x > 2) > (0.11 * length(x)), TRUE) 

data_otu_filt = data.frame(otu_table(data_phylo_filt)) # create a separated file
~~~
{: .language-r}   
  
> ## Question 
> 1. How many zeros are present in the filtered OTU table? 
> 2. What is the percentage of zeros in the filtered OTU table? 
> 3. Plot the number of non-zero values for each OTU (in the same way as for the unfiltered data). Interpret these results.  
> 
> > ## Solution  
> > Q1: How many zeros are present in the filtered OTU table? 
> > ~~~
> > sum(data_otu_filt == 0)
> > ~~~
> > {: .language-r}
> > 
> > Q2. What is the percentage of zeros in the filtered OTU table? 
> > ~~~
> > sum(data_otu_filt == 0) / (dim(data_otu_filt)[2] * dim(data_otu_filt)[1]) * 100
> > ~~~
> > {: .language-r}
> > 
> > Q3. Plot the number of non-zero values for each OTU (in the same way as for the unfiltered data). Interpret these results.  
> > ~~~
> > hist(as.matrix(data_otu_filt), max(data_otu_filt), right = FALSE, las = 1, xlab = "Occurrence value", ylab = "Frequency", main = "Occurrence frequency")
> > min(colSums(data_otu_filt))
> > non_zero <- 0 * 1:dim(data_otu_filt)[2]
> > for (i in 1:dim(data_otu_filt)[2]){
> >   non_zero[i]<-sum(data_otu_filt[,i] != 0)
> >   }
> > plot(non_zero, xlab = "OTU", ylab = "Frequency", main = "Number of non zero values", las = 1)
> > ~~~
> > {: .language-r}
> {: .solution}
{: .challenge}  

You can see that we removed 3867 OTU, which were extremely rare and mostly found in only a few samples. We decreased the number and percentage of zeros in our data set by removing these OTU.
  
## Sequencing depth  
  
Sequencing depth and library size represent the number of counts per sample. This number can vary a lot among the different samples, usually up to 10-fold between some samples mainly due to several technical biases. Some of these biases are explained below.  

After extracting the DNA from each sample independently, the biologist measures the DNA concentration for each sample. First, the sensitivity of the technique used to measure DNA concentration is not always very high. Then, according to this concentration, the biologist pools the different samples in one equimolar mix using a pipette which also adds some error. Moreover, the PCR amplification step can also add some difference among samples because PCR amplification can be more/less efficient according to the sequence itself.  
  
> ## Question
> Why do you think that it is important to look for sequencing depth? Can you find two reasons why?  
> > ## Solution
> > 1. First, it is important to know if we have sequenced enough to catch all the diversity present in our environment. 
> > 2. Secondly, if you want to be able to compare the samples to each other, you will need to have a comparable library size between samples.  
> {: .solution}
{: .challenge}  
  
### Rarefaction curve  
Microbial ecologists explore sequencing depth though a rarefaction curve. 
The rarefaction curve shows how many new OTUs are observed when we obtain new reads for a given sample. If the sequencing depth is sufficient, we should observe a plateau, suggesting that even if we had sequenced new reads, they would belong to OTUs already observed. In other words, all the diversity present in a sample is already described and the bacterial community is sufficiently sequenced. This analysis should be executed on the raw data.  
  
> ## Remark
> We check sequencing depth on the raw data: `data_otu`.  
{: .callout} 

~~~
# the parameter step is used here to decrease computation time 
# step = 100 means that the rarefaction curve will be calculated and plotted only 
# for sample size = 1, 101, 201, (...), up to the maximum value.
rarecurve(data_otu, step = 100, cex = 0.75, las = 1) 
~~~
{: .language-r}
  
We can also customize the plot: adding title, legend, etc.  
~~~
par(cex = 1, las = 1)

leg.txt <- c("Cleron", "Parcey")

lty_vector <- c(rep(x = 1, times = 12), rep(x = 2,times = 12))

rarecurve(data_otu, 
	step = 100, 
	lwd = 1.3, 
	xlab = "Number of reads", 
	ylab = "Number of OTU observed", 
	xlim = c(-50, 23000), 
	ylim = c(-5, 4000), 
	label = FALSE, 
	lty = lty_vector) 

legend(15000, 3900, leg.txt, lty = c(2,1), lwd = 1.3, box.lwd = 0.6, cex = 1)
~~~
{: .language-r}
  
> ## Question   
> How do you interpret this plot? Do you think that the sequencing depth was enough for this experiment? Do you think it is possible to compare the samples using this data set?  
> > ## Solution
> > You can see that none of the samples reaches a plateau which means that sequencing depth was not sufficient but we should not be so far from it. Our previous comparison between Richness and Chao1 gave us the same interpretation.    
> > You can also see with this plot that the total number of reads per sample varies between samples so it will be difficult to compare samples.  
> {: .solution}
{: .challenge}
  
### Library size  
We will now plot the total number of counts per sample.  
~~~
sum_seq <- rowSums(data_otu)
plot(sum_seq, ylim = c(0,25000), main = c("Number of counts per sample"), xlab = c("Samples"))
sum_seq
min(sum_seq)
max(sum_seq)
~~~
{: .language-r}
  
> ## Question 
> 1. How do you interpret this plot? 
> 2. Do you observe similar results on the filtered data set? 
> 3. Do you think it is an issue to have variation in the library size?  
> 
> > ## Solution 
> > 1. We can see that there are differences in the library size for the different samples. The library size go from around 9000 reads to a bit more of 20000 reads (more than 2 fold change).  
> > 2. If you look at the filtered data, you can see similar results: filtering has a small impact on library size because you have only removed a few OTUs
>  that were rare.  
> > 3. If you want to compare samples to each other, it will be an issue to have different library sizes. For example, consider the comparison of:  
> >   - *Sample#1* has 1000 reads in total: 500 reads for OTU1, 300 reads for OTU2 and 200 reads for OTU3.
> >   - *Sample#2* has 100 reads in total:  50 reads for  OTU1, 30 reads for  OTU2 and 20 reads for OTU3.      
> > If you compare the two samples without correcting for library size, you will say that sample#2 has 10 times less of OTU1, OTU2 and OTU3 than sample#1. In this case, OTU abundance differences would only be due to differences in library size. Indeed, if you calculate the percentage, you will have 50% of OTU 1, 30% of OTU 2 and 20% of OTU 3 for both samples. So, it is important to normalize your data per sample.
> {: .solution}
{: .challenge}  
  

> ## Remark
> We observe here a difference of around 2.3-fold change between the lowest and the highest library size. However, usually differences among samples are bigger (up to 10-fold). This is probably due to the fact that we have only a subset of the data here (not all treatments were kept). 
{: .callout} 
  
  
  
## Normalization methods
  
In the next analyses, we will compare counts between samples. To do so, we will first have to normalise the filtered occurrence table by sample in order to obtain the same library size for every sample.  

Different methods exist to normalise microbiome data: proportions and rarefying were commonly used for a long time but other methods were also developed, such as DESeq2 or edgeR‐TMM, which are commonly used in RNA-seq data analyses. While there is still no consensus on which normalisation method should be used for microbiome data, microbial ecologists prefer the proportion or the rarefying methods. Indeed, as explained by McKnight and collaborators [in their paper](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13115), DESeq2 or edgeR‐TMM are recommended based on studies that focused on differential abundance testing and variance standardization, rather than community-level comparisons (*i.e.* beta-diversity). Moreover, standardizing the within-sample variance across samples may suppress differences in species evenness, potentially distorting community level patterns, and the log transformations can exaggerate the importance of differences among rare OTU, while suppressing the importance of differences among abondant OTU.  

For proportion normalisation, each read count is divided by the total sum of all reads of the corresponding sample. As a results, read sums for all the different samples are then equal to 1 and OTU occurrences are between 0 and 1.     

For rarefying normalisation, the library size is arbritary defined for every samples as the smallest library size observed for all the samples (here, `min(sum_seq)`). Then, simple random samplings without replacement are performed for every sample using the raw data.  
  
> ## Discussion
> What do you think are the pros and cons of these two methods?  
{: .discussion}

> ## Remark
> When differences between library sizes are important (e.g. 10-fold difference), it is recommended to use rarefying. As we usually observe 10-fold differences in the library sizes in our dataset, we will normalize using the rarefying normalisation method.    
The normalization should be done on the filtered data, the phyloseq object `data_phylo_filt`.  
{: .callout} 
  
  
Rarefy the data  
~~~
set.seed(1782) # set seed for analysis reproducibility

# rarefy the raw data using the phyloseq package
OTU_filt_rar = rarefy_even_depth(otu_table(data_phylo_filt), rngseed = TRUE, replace = FALSE) 

# create a separated file
data_otu_filt_rar = data.frame(otu_table(OTU_filt_rar)) 

# create a phyloseq object
data_phylo_filt_rar <- phyloseq(OTU_filt_rar, TAX, SAM) 
~~~
{: .language-r}
  
> ## Question
> What is the number of counts per sample for the rarefied data set?  
> 
> > ## Solution
> > ~~~
> > rowSums(data_otu_filt_rar)
> > ~~~
> > {: .language-r}
> > 7750
> {: .solution}
{: .challenge}


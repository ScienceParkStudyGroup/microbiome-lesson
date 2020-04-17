---
title: "Data filtering and normalisation"
teaching: 30
exercises: 0
questions:
- "" 
objectives:
- ""
keypoints:
- ""
---
  
## Data filtering  
  
In the remainder of the data analyses we want to focus on the differences between the samples, the bacterial community composition and the beta-diversity. For this it is necessary to prepare the data in a way that improves the comparability of the samples.  
As microbiome data sets are usually sparse, it is important to filter the data set. Indeed, data filtering aims to remove low quality or uninformative variables (OTU/ASV) to improve downstream statistical analysis.  
For example, we can filter variables (OTU/ASV) with very low number of counts in only a few samples, which are likely due to sequencing errors or low-level contaminations. In this example, we will keep OTU that have at least 2 counts in at least 11% of the samples. Indeed, we have here 18 samples and 3 replicates per treatment and we want, for example, to have at least 2 counts in at least two samples (2/18=0.111).  
  
> ## Remark
>  First, data filtering should be applied on the raw data. Then, we are using here the phyloseq object: *data_phylo* because we are using a plyloseq function to filter 
the raw data.  
{: .callout}
  
~~~
{: .language-r}{r}
data_phylo_filt = filter_taxa(data_phylo, function(x) sum(x > 2) > (0.11*length(x)), TRUE) # filter the OTU data using filter_taxa function included in phyloseq package
data_otu_filt = data.frame(otu_table(data_phylo_filt)) # create a separated file
~~~
{: .language-r}   
  
**QUESTIONS:** How many zeros are present in the filtered OTU table? What is the percentage of zeros in the filtered OTU table? Plot the number of non zero values for each OTU (in the same way as for the nonfiltered data). Interpret these results.  
**SOLUTIONS:**  
~~~
{: .language-r}{r}
sum(data_otu_filt==0)
sum(data_otu_filt==0)/(dim(data_otu_filt)[2]*dim(data_otu_filt)[1])*100
hist(as.matrix(data_otu_filt), max(data_otu_filt), right=FALSE, las=1, xlab = "Occurrence value", ylab = "Frequency", main = "Occurrence frequency")
min(colSums(data_otu_filt))
non_zero<-0*1:dim(data_otu_filt)[2]
for (i in 1:dim(data_otu_filt)[2]){
  non_zero[i]<-sum(data_otu_filt[,i] != 0)
  }
plot(non_zero, xlab = "OTU", ylab = "Frequency", main="Number of non zero values", las=1)
# You can see that we removed 3867 OTU, which were extremly rare and mostly found in only few samples. We deacresed the number and percentage of zeros in our data set by removing these OTU.
~~~
{: .language-r}
  
  
## Sequencing depth  
  
Sequencing depth and library size represent the number of counts per sample. This number can vary a lot among the different samples, usually up to 10 fold between some samples and it is especially due to several technical biais. Some of these biais are explained below.  
After extracting the DNA from each sample independently, the biologist measures the DNA concentration for each sample. First, the sensitivity of the technique used to measure DNA concentration is not always very hight. Then, according to this concentration, the biologist pools the different samples in one equimolar mix using pipette, which add also some error. Moreover, the PCR amplification step can also add some difference among samples because PCR amplification can be more/less efficient according to the sequence itself.  
  
**QUESTION:** Why do you think that it is important to look for sequencing depth?  
**SOLUTION:** First, it is important to know if we sequenced enough to catch all the diversity present in our environment. Then, if you want to be able to compare the samples to each others, we need to have a comparable library size between samples.  
  
  
### Rarefaction curve  
Microbial ecologists explore sequencing depth though a rarefaction curve. The rarefaction curve shows how many new OTU are observed when we obtain new reads for a given sample. If the sequencing depth is enough, we should observe a plateau, meaning that even if we sequence new reads they will belong to OTUs already observed. In other word, all the diversity present in a sample is already described and we have sequenced the community deeply enough. This analysis should be execute on the raw data.  
  
> ## Remark 
> We check sequencing depth on the raw data: *data_otu*.  
{: .callout}

~~~
rarecurve(data_otu, step = 100, cex=0.75, las=1) # the parameter step is used here to decrease computation time - step = 100 means that the rarefaction curve will be calculated and plotted only for sample size = 1, 101, 201, (...), and maximun value.
~~~
{: .language-r}
  
We can also customize the plot: adding title, legend, etc.  
~~~
par(cex=1, las=1)
leg.txt <- c("Cleron", "Parcey")
lty_vector <- c(2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1)
rarecurve(data_otu, step = 100, lwd=1.3, xlab = "Number of reads", ylab = "Number of OTU observed", xlim = c(-50, 23000), ylim = c(-5, 4000), label = FALSE, lty = lty_vector) 
legend(15000, 3900, leg.txt, lty=c(2,1), lwd=1.3, box.lwd=0.6, cex=1)
~~~
{: .language-r}
  
**QUESTIONS:** How do you interpret this plot? Do you think that the sequencing depth was enough for this experiment? Do you think it is possible to compare the samples using this data set?  
**SOLUTIONS:** You can see that none of the samples reach a plateau, so the sequencing depth was not enough, but we should not be so far from it. Our previous comparison between Richness and Chao1 give us the same interpretation.  
You can also see with this plot that the total number of reads per sample vary between samples, so it will be difficult to compare samples.  
  
  
### Library size  
We will now plot the total number of counts per sample.  
~~~
sum_seq <- rowSums(data_otu)
plot(sum_seq, ylim=c(0,25000), main=c("Number of counts per sample"), xlab=c("Samples"))
sum_seq
min(sum_seq)
max(sum_seq)
~~~
{: .language-r}
  
**QUESTIONS:** How do you interpret this plot? Do you observe similar results on the filtered data set? Do you think it is an issue to have variation in the library size?  
**SOLUTIONS:** We can see that there is differences in the library size for the different samples. The library size go from around 9000 reads to a bit more of 20000 reads (more than 2 fold change).  
If you lokk at the filtered data, you can see similar result, the filtering  does not change library size so much because you just removed few OTU that were rare.  
Of course if you want to compare samples to each other, it will be an issue to have different library sizes. For example, if you have a sample 1 with 1000 reads (500 reads for OTU1, 300 reads for OTU 2 and 200 reads for OTU 3) and a sample 2 with only 100 reads (50 reads for OTU1, 30 reads for OTU 2 and 20 reads for OTU 3). If you are comparing the two samples without correcting for the library size, you will say that sample 2 as 10 times less OTU 1, 2 and 3 than the sample 1, while it is just due to the difference in library size. Indeed, if you calculate the percentage, you will have 50% of OTU 1, 30% of OTU 2 and 20% of OTU 3 for both samples. So, it is important to normalize your data per sample.  
  
> ## Remark 
> We observe here a diffrence of around 2.3 fold change between the lowest and the highest library size. However, usually differences among samples are bigger (up to 10 fold). This is probably due to the fact that we have only a subset of the data here (not all the treatments were kept).  
 {: .callout} 
  
  
## Normalization per sample  
  
In the further analyses, we will compare counts between samples. To do so, we will have first to normalise the filtered occurrence table by sample in order to obtain the same library size for every samples.  
Different methods exist to normalise microbiome data: proportions and rarefying were commonly used for long time but other methods were also developed, such as DESeq2 or edgeR‐TMM, which are commonly used in RNA-seq data analyses. While there is still no consensus on which normalisation method should be used for microbiome data, microbial ecologists prefer the proportion or the rarefying methods. Indeed, as explained McKnight and colaborators (DOI: 10.1111/2041-210X.13115) DESeq2 or edgeR‐TMM are recommended based on studies that focused on differential abundance testing and variance standardization, rather than community-level comparisons (*i.e.* beta-diversity). Moreover, standardizing the within-sample variance across samples may suppress differences in species evenness, potentially distorting community level patterns, and the log transformations can exaggerate the importance of differences among rare OTU, while suppressing the importance of differences among abondant OTU.  
For proportion, each read count is divided by the total sum of all reads of the corresponding sample. As a results, read sums for all the different samples are then equal to 1 and OTU occurrences are between 0 and 1.  
For rarefying, the library size is arbritary defined for every samples as the smallest library size observed for all the samples (here, *min(sum_seq)*). Then, simple random samples without replacement are performed for every sample using the raw data.  
  
**QUESTION:** What do you think are the pros and cons of these two methods?  
  
> ## Remark 
> When differences between library sizes is high (such as 10 fold change), it is recommended to use rarefying. As we usually observe 10 fold change in the library sizes, we will normalize here using the rarefying method.    
> The normalization should be done on the filtered data, the phyloseq object *data_phylo_filt*.  
 {: .callout}
  
Rarefy the data. 
~~~
set.seed(1782) # set seed for analysis reproducibility
OTU_filt_rar = rarefy_even_depth(otu_table(data_phylo_filt), rngseed=T, replace=FALSE) # rarefy the raw data using Phyloseq package
data_otu_filt_rar = data.frame(otu_table(OTU_filt_rar)) # create a separated file
data_phylo_filt_rar <- phyloseq(OTU_filt_rar, TAX, SAM) # create a phyloseq object
~~~
{: .language-r}
  
**QUESTION:** What is the number of counts per sample for the rarfied data set?  
**SOLUTIONS:** 7750

~~~
rowSums(data_otu_filt_rar)
~~~
{: .language-r}  
---
title: "Data preprocessing: from fastq to ASV"
teaching: 60
exercises: 30
questions:
- How can I go from fasta files to ASV data?  
objectives:
- Check the reads' quality  
- Trim and filter the fasta files  
- Denoise the data  
- Merge the forward and reverse reads  
- Create an ASV table  
- Remove chimeras  
- Assign taxonomic affiliation for each ASV  
keypoints:
- ""
---

## Table of Contents  
- [1. Install and load the required packages  ](#1-install-and-load-the-required-packages)
- [2. Getting ready  ](#2-getting-ready)
- [3. Inspect read quality profiles  ](#3-inspect-read-quality-profiles)
- [4. Filter and trim  ](#4-filter-and-trim)
- [5. Denoising  ](#5-denoising)
- [6. Merge paired reads  ](#6-merge-paired-reads)
- [7. Construct sequence table  ](#7-construct-sequence-table)
- [8. Remove chimeras  ](#8-remove-chimeras)
- [9. Track reads through the pipeline  ](#9-track-reads-through-the-pipeline)
- [10. Assign taxonomy  ](#10-assign-taxonomy)
  
  
## 1. Install and load the required packages  
  
~~~
library(dada2)
packageVersion("dada2")
~~~
{: .language-r}
  
## 2. Getting ready  
  
In the Shell unzip the fasta file.  
untar("MiSeq_data.tar.gz")  
  
~~~
# Set the path
path <- "~/MiSeq_data"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
~~~
{: .language-r}
  
## 3. Inspect read quality profiles  
  
~~~
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
~~~
{: .language-r}
  
## 4. Filter and trim  
  
~~~
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
~~~
{: .language-r}
  
## 5. Denoising 
  
~~~
# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
~~~
{: .language-r}
  
## 6. Merge paired reads  
  
~~~
# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])
~~~
{: .language-r}
  
## 7. Construct sequence table  
  
~~~
# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
~~~
{: .language-r}
  
## 8. Remove chimeras  
  
~~~
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
~~~
{: .language-r}
  
## 9. Track reads through the pipeline  
  
~~~
# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
~~~
{: .language-r}
  
## 10. Assign taxonomy  
  
~~~
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
~~~
{: .language-r}
  
  

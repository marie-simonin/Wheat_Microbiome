---
title: "Wheat Microbiome paper - Bioinformatic processing"
author: "Marie Simonin"
date: "7/31/2019"
output: html_document
---

  1. Processing of 16S rRNA gene dataset from raw data
  2. Processing of 18S rRNA gene dataset from raw data


# 1. Processing of 16S rRNA gene dataset from raw data
## In Qiime2
## Import fastq files in Qiime2
```{r}
source activate qiime2-2019.1

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 'Data to import' \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path 16S-demux-paired-end.qza
```

## Visualize and verify sequence quality (total number of sequences = 14 799 693). Open qzv file on Qiime2View online
```{r}
qiime demux summarize \
  --i-data 16S-demux-paired-end.qza \
  --o-visualization 16S-demux-paired-end.qzv
```


## Sequence Denoising with DADA2 (after DADA2, number of SV = 24515 and total number of sequences = 5 793 329)
```{r}
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16S-demux-paired-end.qza \
  --p-trim-left-f 18 \
  --p-trim-left-r 18 \
  --p-trunc-len-f 239 \
  --p-trunc-len-r 208 \
  --o-table 16S-table-2.qza \
  --p-n-threads 0 \
  --o-representative-sequences 16S-rep-seqs-2.qza \
  --o-denoising-stats 16S-denoising-stats-2.qza
```


## Visualize the outputs of DADA2 in Qiime2View online
```{r}
qiime metadata tabulate \
  --m-input-file 16S-denoising-stats-2.qza \
  --o-visualization 16S-denoising-stats-2.qzv

qiime feature-table summarize \
  --i-table 16S-table-2.qza \
  --o-visualization 16S-table-2.qzv \
  --m-sample-metadata-file metadata_16S.txt

qiime feature-table tabulate-seqs \
  --i-data 16S-rep-seqs-2.qza \
  --o-visualization 16S-rep-seqs-2.qzv
```

## Export SV table (biom file) and representative sequences (fasta file)
```{r}
qiime tools export \
  --input-path 16S-rep-seqs-2.qza \
  --output-path Dada2-output

qiime tools export \
  --input-path 16S-table-2.qza \
  --output-path Dada2-output
```


## In Qiime 1
## convert biom file in txt file in Qiime 1
```{r}
source activate qiime1
biom convert -i 16S-SV-table.biom -o 16S-SV-table.txt --to-tsv --table-type="OTU table"
```

## Assign taxonomy with Qiime1 using SILVA 132 database
```{r}
assign_taxonomy.py -i 16S_rep-seq.fasta.txt -r silva_132_99_16S.fna -t 16S_consensus_taxonomy_7_levels.txt
```

## Remove Chloroplast, Mitochondria, Eukaryotic SVs and SVs with low abundance (present in only 1 sample and/or with less than 10 total observation count) from SV table and rep-seq fasta file

## Remove Contaminant SV using negative controls and Decontam package in R

## Perform rarefaction on filtered SV table (non-microbial, low abundance sequences and contaminants removed) and export the results as a txt file
```{r}
single_rarefaction.py -i 16S-SV-table-filtered.biom -o 18S-SV-table-filtered-rarefied1018.biom -d 1018

biom summarize-table -i 16S-SV-table-filtered-rarefied1018.biom

biom convert -i 16S-SV-table-filtered-rarefied1018.biom -o 16S-SV-table-filtered-rarefied1018.txt --to-tsv
```


## In Qiime2
## Import filtered rep-seq fasta file (non-microbial, low abundance sequences and contaminants removed) in Qiime2 to make phylogenetic tree
```{r}
source activate qiime2-2019.1
qiime tools import \
  --input-path 16S-rep-seqs-filtered.fasta \
  --output-path 16S-rep-seqs-filtered.qza \
  --type 'FeatureData[Sequence]'
```

## Make phylogenic tree in Qiime2
```{r}
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 16S-rep-seqs-filtered.qza \
  --o-alignment aligned-16S-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-16S-tree-allSVs.qza \
  --o-rooted-tree rooted-16S-tree-allSVs.qza
```

## Export trees
```{r}
qiime tools export \
  --input-path unrooted-16S-tree-allSVs.qza \
  --output-path exported-tree-16S

qiime tools export \
  --input-path rooted-16S-tree-allSVs.qza \
  --output-path exported-tree-16S
```





######################################################################################

# 2. Processing of 18S rRNA gene dataset from raw data
## In Qiime2
## Import fastq files in Qiime2
```{r}
source activate qiime2-2019.1

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 'Data to import' \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path 18S-demux-paired-end.qza
```

## Visualize and verify sequence quality (total number of sequences = 14 799 693). Open qzv file on Qiime2View online
```{r}
qiime demux summarize \
  --i-data 18S-demux-paired-end.qza \
  --o-visualization 18S-demux-paired-end.qzv
```


## Sequence Denoising with DADA2 (after DADA2, number of SV = 24515 and total number of sequences = 5 793 329)
```{r}
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 18S-demux-paired-end.qza \
  --p-trim-left-f 18 \
  --p-trim-left-r 18 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 190 \
  --o-table 18S-table-2.qza \
  --p-n-threads 0 \
  --o-representative-sequences 18S-rep-seqs-2.qza \
  --o-denoising-stats 18S-denoising-stats-2.qza
```


## Visualize the outputs of DADA2 in Qiime2View online
```{r}
qiime metadata tabulate \
  --m-input-file 18S-denoising-stats-2.qza \
  --o-visualization 18S-denoising-stats-2.qzv

qiime feature-table summarize \
  --i-table 18S-table-2.qza \
  --o-visualization 18S-table-2.qzv \
  --m-sample-metadata-file metadata_18S.txt

qiime feature-table tabulate-seqs \
  --i-data 18S-rep-seqs-2.qza \
  --o-visualization 18S-rep-seqs-2.qzv
```

## Export SV table (biom file) and representative sequences (fasta file)
```{r}
qiime tools export \
  --input-path 18S-rep-seqs-2.qza \
  --output-path Dada2-output

qiime tools export \
  --input-path 18S-table-2.qza \
  --output-path Dada2-output
```


## In Qiime 1
## convert biom file in txt file in Qiime 1
```{r}
source activate qiime1
biom convert -i 18S-SV-table.biom -o 18S-SV-table.txt --to-tsv --table-type="OTU table"
```

## Assign taxonomy with Qiime1 using Protist Ribosomal Reference database (PR2, Guillou et al. 2013) 
```{r}
assign_taxonomy.py -i 18S_rep-seq.fasta.txt -r pr2_version_4.11.1_mothur.fasta -t pr2_version_4.11.1_mothur_tax.txt
```

##Filter Chloroplast, Mitochondria, Chloroplastida and Animalia SVs and SVs with low abundance (present in only 1 sample and/or with less than 10 total observation count) from SV table and rep-seq fasta file

## Remove Contaminant SV using negative controls and Decontam package in R

## Perform rarefaction on filtered SV table (non-microbial, low abundance sequences and contaminants removed) and export the results as a txt file
```{r}
single_rarefaction.py -i 18S-SV-table-filtered.biom -o 18S-SV-table-filtered-rarefied1127.biom -d 1127

biom summarize-table -i 18S-SV-table-filtered-rarefied1127.biom

biom convert -i 18S-SV-table-filtered-rarefied1127.biom -o 18S-SV-table-filtered-rarefied1127.txt --to-tsv
```


## In Qiime2
## Import filtered rep-seq fasta file (non-microbial and low abundance sequences removed) in Qiime2 to make phylogenetic tree
```{r}
source activate qiime2-2019.1
qiime tools import \
  --input-path 18S-rep-seqs-filtered.fasta \
  --output-path 18S-rep-seqs-filtered.qza \
  --type 'FeatureData[Sequence]'
```

## Make phylogenic tree in Qiime2
```{r}
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 18S-rep-seqs-filtered.qza \
  --o-alignment aligned-18S-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-18S-tree-allSVs.qza \
  --o-rooted-tree rooted-18S-tree-allSVs.qza
```

## Export trees
```{r}
qiime tools export \
  --input-path unrooted-18S-tree-allSVs.qza \
  --output-path exported-tree-18S

qiime tools export \
  --input-path rooted-18S-tree-allSVs.qza \
  --output-path exported-tree-18S
```

# Alyssa Pybus
# Script for pre-processing data ################

rm(list=ls())   # clears global environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # tells R that the folder where this code resides is our working directory

pacman::p_load(tidyverse)


# TRANSCRIPT DATA FILTERING AND SORTING ####################################

# read in experiment data and metadata
metadata.us = readxl::read_excel("Data/22131-02 metadata.xlsx") %>%
  mutate(Injury = factor(Injury,levels=c("Sham","1xCHI","pre-3xCHI","3xCHI","pre-5xCHI","5xCHI")))
count.in = rio::import("Data/22131-02-gene-count-all.csv")
gene.df.uf = count.in[,2:3]  # unfiltered gene names (non-unique) and EnsemblIDs (unique)
count.uf = count.in[,4:42]  # unfiltered gene counts
colnames(count.uf) = metadata.us$Sample
rownames(count.uf) = gene.df.uf$EnsemblID

# filter criteria as recommended on WGCNA FAQ (genes where 90% of samples have count >10)
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html para 4
keepInd=which(rowSums(count.uf>10) > {ncol(count.uf)*.9})
count.us = count.uf[keepInd,]   # unsorted columns (samples), filtered genes (rows)
gene.df = gene.df.uf[keepInd,]  # filtered genes (rows)

# sort by AbsoluteTime
metadata = dplyr::arrange(metadata.us,AbsoluteTime)  # sort sample metadata by ascending time point since first injury
sort_ind = match(metadata$Sample,metadata.us$Sample)
count <- count.us[,sort_ind]  # sort count matrix columns by ascending time point since first injury (match metadata tibble)

save(count,metadata,gene.df,file="R Data/genecounts.rda")  # save to access later


# PROTEIN DATA FILTERING AND SORTING ####################################

# load protein data
df_in <- readxl::read_excel("Data/DoD TBI Acute All.xlsx") %>%
  dplyr::arrange(Region,AbsoluteTime) %>%
  dplyr::filter(Sex == "F") %>%  # Female samples only (RNAseq was only on female samples)
  dplyr::filter(Region == "FC")

# define the analytes to treat as numeric variables by setting first and last in range
first_analyte="GFAP"
last_analyte="TNF-a"
analyte_range={which(colnames(df_in)==first_analyte)}:{which(colnames(df_in)==last_analyte)}

#create array of average background for each analyte
df_b <- df_in %>%
  filter(Sample == "background") %>%
  dplyr::select(all_of(analyte_range)) %>%
  colMeans()

#create data frame of sample data
df_samples <- df_in %>%
  dplyr::filter(Sample != "background") %>%
  dplyr::filter(Sample %in% metadata$Sample) %>%
  mutate(Sample = as.numeric(Sample))

# Subtract background, set negative values to zero
df_samples[,analyte_range] <- sweep(df_samples[,analyte_range],2,df_b)
df_samples[df_samples<0] <- 0 # analytes below background set to zero

# create unicode labels for greek characters
greek = colnames(df_samples)
greek[which(greek=="AB40")] = "A\u03b240"
greek[which(greek=="AB42")] = "A\u03b242"
greek[which(greek=="AB42/40")] = "A\u03b242/40"
greek[which(greek=="tTau")] = "Total Tau"
greek[which(greek=="pTau")] = "Phospho-Tau T181"
greek[which(greek=="IFN-g")] = "IFN-\u03b3"
greek[which(greek=="IL-1a")] = "IL-1\u03b1"
greek[which(greek=="IL-1b")] = "IL-1\u03b2"
greek[which(greek=="MIP-1a")] = "MIP-1\u03b1"
greek[which(greek=="MIP-1b")] = "MIP-1\u03b2"
greek[which(greek=="TNF-a")] = "TNF-\u03b1"

# create color bars for use in figures
InjPalette = rev(RColorBrewer::brewer.pal(6,"Spectral"))
TimePalette = RColorBrewer::brewer.pal(4,"Purples")

InjColors = case_when(
  df_samples$Injury=="Sham" ~ InjPalette[1],
  df_samples$Injury== "1xCHI" ~ InjPalette[2],
  df_samples$Injury== "pre-3xCHI" ~ InjPalette[3],
  df_samples$Injury== "3xCHI" ~ InjPalette[4],
  df_samples$Injury== "pre-5xCHI" ~ InjPalette[5],
  df_samples$Injury== "5xCHI" ~ InjPalette[6])

TimeColors = case_when(
  df_samples$TimePoint== 0 ~ TimePalette[1],
  df_samples$TimePoint== 0.5 ~ TimePalette[2],
  df_samples$TimePoint== 4 ~ TimePalette[3],
  df_samples$TimePoint== 24 ~ TimePalette[4],)


save(df_samples,analyte_range,greek,InjPalette,TimePalette,InjColors,TimeColors,file="R Data/proteins.rda") #  save to access later


# LOAD AND FORMAT GSVA GENE SETS ####################################

# MSigDB c2
geneSetsIn=readxl::read_excel("Reference Files/c2.cp.v7.0.symbols.xlsx",col_names=FALSE) %>%
  dplyr::select(-2) %>%
  t() 
colnames(geneSetsIn) = geneSetsIn[1,]
gs_MSigDB_C2 = geneSetsIn[-1,] %>%
  as.data.frame() %>%
  as.list() %>%
  lapply(function(x) x[!is.na(x)])

# Astrocyte Annotations
gs_astro=readxl::read_excel("Reference Files/GeneSets_Astrocytes_FXN1.xlsx") %>%
  as.list() %>%
  lapply(function(x) x[!is.na(x)])

# Neuron Annotations
gs_neuron=readxl::read_excel("Reference Files/GeneSets_Neurons_FXN1.xlsx") %>%
  as.list() %>%
  lapply(function(x) x[!is.na(x)])

save(gs_MSigDB_C2,gs_astro,gs_neuron,file="R Data/gs.rda")


# CALCULATION OF DESEQ2 NORMALIZED / VARIANCE STABILIZED GENE DATA FOR WGCNA ####################################

rm(list=ls())   # clears "workspace"
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # tells R that the folder where this code resides is our working directory

load("R Data/genecounts.rda")

# DESeq2 Normalization 
des_mat = count %>%
  round(digits=0) %>%
  as.matrix()
metadata$AbsoluteTime=factor(metadata$AbsoluteTime)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = des_mat,
                              colData = metadata,
                              design = ~ AbsoluteTime)
dds <- DESeq2::DESeq(dds)

# Variance Stabilizing Transform
vsd <- DESeq2::getVarianceStabilizedData(dds)
datExpr = t(vsd)

cleanDat = vsd
rownames(cleanDat) = str_c(gene.df$geneName,"|",rownames(vsd))

save(datExpr,cleanDat,metadata,gene.df,file="R Data/geneVSD.rda")

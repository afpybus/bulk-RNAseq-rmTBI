# File:    WGCNA.R
# Project:  DoD rmTBI Acute Response

# PRELIMINARIES ################################
rm(list=ls())   # clears "workspace"
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # tells R that the folder where this code resides is our working directory

if (!require("pacman")) install.packages("pacman")  # install pacman if you don't have it already
if (!require("BiocManager")) install.packages("BiocManager") # install BiocManager if you don't have it already
if (!require("WoodLabFunctions")) pacman::p_load_gh("afpybus/WoodLabFunctions")
#BiocManager::install("DESeq2")
#BiocManager::install("GSVA")
#BiocManager::install("ropls")

# Load contributed packages with pacman
pacman::p_load(pacman,rio,tidyverse,readxl,matrixStats,ggpubr,heatmap3,Rtsne,Hmisc,biomaRt,
               gplots,GSVA,ropls,RColorBrewer,DESeq2,gridExtra,WGCNA,Rtsne) # install/load these packages
# pacman: for package handling
# rio: for smart import/export functions
# tidyverse: for handling tibble data frames and data wrangling
# readxl: for excel importing
# matrixStats: for matrix operations
# ggpubr: for publication-ready figures
# heatmap3: our favorite heatmap function
# Rtsne: t-distributed stochastic neighbor embedding
# gplots: more graphical things for figures
# GSVA: gene set variation analysis, Bioconductor
# ropls: for PCA, PLSR, PLSDA, etc., Bioconductor
# RColorBrewer: for color palettes
# DESeq2: for normalizing RNAseq read counts, Bioconductor


load("R Data/geneVSD.rda")
load("R Data/proteins.rda")

# # Run WGCNA #################
# # parameters adopted from Johnson et al, Nature Neuroscience (2022) (Dammer, Lah, Levey, Seyfried)
# if (!require("WGCNA")) install.packages("WGCNA")
# net = blockwiseModules(datExpr, power = 4, networkType = "signed", 
#                        deepSplit = 4,
#                        minModuleSize = 20, 
#                        TOMDenom = "mean", 
#                        corType = "bicor",
#                        mergeCutHeight = 0.07, 
#                        reassignThreshold = 0.05,
#                        numericLabels = T, verbose = 3,maxBlockSize = 15000)
# saveRDS(net,"R Data/net.RDS")
# # running WGCNA after loading later packages will mask cor() and cause issues





#### All Time Points WGCNA ################

# png("WGCNA Output/hm.png",units="in",height=10,width=7,res=600)
# heatmap_wl(cleanDat,clust_c = F,clust_r=F)
# dev.off()

# # Sample Cluster and Threshold Power ####
# # ~~~ Show clustered samples to detect possible outliers ~~~~~~~~~
# sampleTree = hclust(dist(datExpr), method = "average")
# png(paste0(output_folder,"Sample Cluster.png"),height=5,width=9,units="in",res=600)
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# dev.off()
# 
# # r1 = outlierRemoval(datExpr,0.001)
# # rownames(datExpr)[r1$removeIndex]
# 
# oplsOut = opls(datExpr,predI=2) %>%
#   rotate_opls(degrees=2)
# scores_plot(oplsOut$T1,oplsOut$T2,color=metadata$Injury,analysis = "PCA")
# # Sample 247,244 are possible outliers, use mahalanobis detection to decide whether to exclude


# # ~~~ Select threshold power ~~~~~~~~~
# sft = pickSoftThreshold(datExpr, verbose = 5)
# png("WGCNA Output/ScaleFreeThresholdPower.png",height=5,width=5,units="in",res=1000)
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=sft$fitIndices$Power,cex=0.9,col="red");
# abline(h=0.80,col="red")
# dev.off()
# # A threshold power of 4 was chosen since it was the smallest threshold that resulted in a scale-free local max R^2 value greater than 0.80
# 
# rio::export(sft$fitIndices,"WGCNA Output/SoftThresholdPowers.csv")


# Load in WGCNA Net ~~~~~~~~~~~~~~~ ####

net = readRDS("R Data/net.RDS")
table(net$colors)
og_ME_nums =net$colors
mergedColors = labels2colors(net$colors)
net[["colors"]] = mergedColors


# Export Lists of Genes (NOT for Gene Ontology, those are filtered for kME>0.60) ##############

df=data.frame(gene=gene.df[,2],color=net$colors)
wb=xlsx::createWorkbook()
sh=xlsx::createSheet(wb, "Sheet1")
for (i in 1:length(unique(net$colors))){
  df_i = filter(df,color==unique(net$colors)[i]) %>% dplyr::select(gene)
  colnames(df_i)=unique(net$colors)[i]
  xlsx::addDataFrame(df_i,startColumn=i,sheet=sh,row.names = FALSE)
}
xlsx::saveWorkbook(wb, paste0("WGCNA Output/geneList.xlsx"))


# Dendro plots ~~~~~~~~ ####

# Plot the dendrogram and the module colors underneath

pdf("WGCNA Output/Dendrogram.pdf",height=3.5,width=7)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    dendroLabels = FALSE,hang=0.03,addGuide = TRUE,guideHang = 0.05)
dev.off()

png("WGCNA Output/Dendrogram.png",height=3.5,width=7,units="in",res=1000)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    dendroLabels = FALSE,hang=0.03,addGuide = TRUE,guideHang = 0.05)
dev.off()

# ~~~~~ Eigengene Network / Module Eigengenes ~~~~~~~~~~~~~~~~ ####

sampleTree = hclust(dist(t(net$MEs)), method = "average")
MEs = net$MEs[,sampleTree$order]

identical(as.character(df_samples$Sample),rownames(MEs))  # make sure this is true
identical(as.character(metadata$Sample),rownames(MEs))  # make sure this is true

MEcolors = matrix(nrow=length(unique(mergedColors)),ncol=3)
for (i in 1:nrow(MEcolors)){
  module = unique(og_ME_nums)[i]
  MEcolors[i,2] = module
  MEcolors[i,1] = unique(mergedColors[which(og_ME_nums==module)])
}
MEcolors[,3]=str_c("ME",MEcolors[,2])
sampleTree$labels = MEcolors[,1][match(sampleTree$labels,MEcolors[,3])]
colnames(MEs) = str_c("ME",MEcolors[,1][match(colnames(MEs),MEcolors[,3])])


# ~~~~~~~~~ ME Expression Heatmaps ~~~~~~~~~~~~~~~~~~ ####

ME_hm = as_tibble(MEs) %>%
  mutate(Sample = as.numeric(rownames(MEs))) %>%
  left_join(metadata) %>%
  arrange(AbsoluteTime)

hm_mat = ME_hm[,1:ncol(MEs)] %>% dplyr::select(-MEgrey) %>% as.matrix()
hclust_out = hclust(dist(t(hm_mat)), method = "average")

hm_matZ = apply(hm_mat,2,"scale")

#Set heatmap3 color bar parameters
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red") # create spectrum of blue to red for heatmaps

pdf("WGCNA Output/Z_colorbar.pdf", width=2,height=5,pointsize = 10)
color.bar(barColors,-1.5,1.5,nticks=7)
dev.off()

# png(paste0(output_folder,"ME_ExpGroup_heatmapZ.png"),height=7,width=10,units="in",res=1000)
pdf("WGCNA Output/ME_ExpGroup_heatmapZ.pdf",height=7,width=8)
heatmap3(hm_matZ,
         Rowv=NA,Colv=as.dendrogram(hclust_out),cexCol = 1.5,
         scale = "none",col=barColors, breaks=breakBarColors,margins=c(5,15),
         legendfun=function()showLegend(legend=c(NA),col=c(NA)),labRow=NA,
         labCol =  MEcolors[,3][match(gsub("ME","",colnames(hm_matZ)),MEcolors[,1])],
         RowSideColors = cbind(c(Injury=InjColors),c(Time=TimeColors)),ColSideWidth=1,RowSideLabs=NA,
         ColSideColors = cbind(c(substring(colnames(hm_mat),3))),ColSideLabs = NA)
dev.off()

png("WGCNA Output/ME_ExpGroup_heatmapZ.png",height=7,width=8,units="in",res=1000)
heatmap3(hm_matZ,
         Rowv=NA,Colv=as.dendrogram(hclust_out),cexCol = 1.5,
         scale = "none",col=barColors, breaks=breakBarColors,margins=c(5,15),
         legendfun=function()showLegend(legend=c(NA),col=c(NA)),labRow=NA,
         labCol =  MEcolors[,3][match(gsub("ME","",colnames(hm_matZ)),MEcolors[,1])],
         RowSideColors = cbind(c(Injury=InjColors),c(Time=TimeColors)),ColSideWidth=1,RowSideLabs=NA,
         ColSideColors = cbind(c(substring(colnames(hm_mat),3))),ColSideLabs = NA)
dev.off()

# ~~~~~~~~~ ME Expression Line Graphs ~~~~~~~~~~~~~~~~~~ ####

# dir.create("MEs by Injury/",recursive = TRUE,showWarnings = FALSE)
# for(i in 1:ncol(MEs)){
#   png(paste0("MEs by Injury/",colnames(MEs)[i],".png"),height=4,width=5,units="in",res=1000)
#   # pdf(paste0(output_folder,"/MEs by Injury/",colnames(MEs)[i],".pdf"),height=4,width=5)
#   print(new_meanSEM(x=metadata$Injury,y=MEs[,i],title = colnames(MEs)[i],axis.text.x.size = 12,axis.text.y.size = 12)%>%
#           ggpar(palette=InjPalettePar,legend="none"))
#   dev.off()
# }

ME_tib = as_tibble(MEs) %>%
  mutate(Sample = as.numeric(rownames(MEs))) %>%
  left_join(metadata) %>%
  gather(key="ME",value="Score",1:ncol(MEs)) %>%
  mutate(color_code = substring(ME,3)) 

ME_tib$ME_num = MEcolors[,3][match(ME_tib$color_code,MEcolors[,1])]

ME_sub = ME_tib %>%
  arrange(ME_num) %>%
# filter(ME %in% "MEsalmon")
# filter(ME %in% c("MEblue","MEturquoise"))
# filter(ME %in% c("MEred","MEbrown"))
# filter(ME %in% c("MEred","MEbrown","MEsalmon"))
# filter(ME %in% c("MEgreenyellow","MEpurple"))
# filter(ME %in% c("MEgreen","MEtan"))
# filter(ME %in% c("MEmagenta","MEpink","MEgrey","MEblack","MEyellow"))
# filter(ME %in% c("MEmagenta","MEpink","MEblack"))
# filter(ME %in% c("MEcyan"))
# filter(ME %in% c("MEblue"))
filter(ME %in% c("MEturquoise"))

# png(paste0("WGCNA Output/ME by Inj/",paste0(unique(ME_sub$ME),collapse=" "),".png"),height=4,width=7,units="in",res=1000)
pdf(paste0("WGCNA Output/ME by Inj/",paste0(unique(ME_sub$ME),collapse=" "),".pdf"),height=2,width=3)
ggline(ME_sub,x = "Injury",y="Score",color="ME_num",add=c("mean_se"),xlab = "",
       palette = unique(ME_sub$color_code))+theme(legend.title = element_blank(),text = element_text(size=6))
dev.off()

compare_means(formula = Score ~ Injury, data = ME_sub,paired = FALSE,group.by = "color_code",ref.group = "Sham",method="t.test")


df_perm = ME_tib %>%
  filter(Injury %in% c("Sham","5xCHI"),color_code=="red")

sham_true = filter(df_perm,Injury=="Sham") %>% dplyr::select(Score)
hit_true = filter(df_perm,Injury=="5xCHI") %>% dplyr::select(Score)
true_dif = colMeans(hit_true) - colMeans(sham_true)
scores = df_perm$Score

perm_res = matrix(1:10000)
for(i in 1:10000){
  sham_perm = sample(scores,count(df_perm$Injury=="Sham"))
  hit_perm = scores[!{scores %in% sham_perm}]
  avg_dif = mean(hit_perm) - mean(sham_perm)
  
  perm_res[i] = true_dif <= avg_dif
}
sum(perm_res) / 10000



# ind = which(net$colors=="turquoise")[5]
# error_plot(x=metadata$Injury,y=cleanDat[ind,],title = rownames(cleanDat)[ind])
# 
# ind = which(net$colors=="blue")[5]
# error_plot(x=metadata$Injury,y=cleanDat[ind,],title = rownames(cleanDat)[ind])

# dir.create(paste0(output_folder,"MEs by Time Point/"),recursive = TRUE,showWarnings = FALSE)
# for(i in 1:ncol(MEs)){
#   png(paste0(output_folder,"/MEs by Time Point/",colnames(MEs)[i],".png"),height=4,width=5,units="in",res=1000)
#   # pdf(paste0(output_folder,"/MEs by Time Point/",colnames(MEs)[i],".pdf"),height=4,width=5)
#   print(new_meanSEM(x=factor(metadata$AbsoluteTime),y=MEs[,i],colors = metadata$Injury,title = colnames(MEs)[i],axis.text.x.size = 12,axis.text.y.size = 12)%>%
#           ggpar(palette=InjPalettePar,legend="none"))
#   dev.off()
# }
# pdf(paste0(output_folder,"/MEs by Time Point/",paste0(unique(ME_sub$ME),collapse=" "),".pdf"),height=4,width=7)
# ggline(ME_sub,x = "AbsoluteTime",y="Score",color="ME",add="mean_se",palette = sort(unique(ME_sub$color_code)))
# dev.off()


# ~~~~~ Trait Heatmaps ~~~~~~~~~~~~~~~~ ####

datTraits = dplyr::select(df_samples,`L hemi 3mm`,GFAP:`TNF-a`) %>% 
  as.matrix()
rownames(datTraits) = df_samples$Sample

MEs_nogrey = dplyr::select(MEs,-MEgrey)

moduleTraitCor = WGCNA::cor(MEs_nogrey, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs_nogrey));
# FDR_p = apply(moduleTraitPvalue,2,function(x)p.adjust(x,method="fdr"))
# which(FDR_p<0.05)

# c = floor((which(abs(moduleTraitCor)>0.5)-1) / nrow(moduleTraitCor)) + 1
c = floor((which(moduleTraitPvalue<0.01)-1) / nrow(moduleTraitCor)) + 1

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

png("WGCNA Output/traitHeatmap.png",height=7,width=12,units="in",res=1000)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs_nogrey),
               ySymbols = names(MEs_nogrey),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

justCor = moduleTraitCor[,unique(c)]
justCorP = moduleTraitPvalue[,unique(c)]

# Will display correlations and their p-values
textMatrix2 = paste(signif(justCor, 2), "\n(",
                   signif(justCorP, 1), ")", sep = "");
dim(textMatrix2) = dim(justCor)

png("WGCNA Output/traitHeatmap_justCor.png",height=7,width=10,units="in",res=1000)
labeledHeatmap(Matrix = justCor,
               xLabels = colnames(justCor),
               yLabels = names(MEs_nogrey),
               ySymbols = names(MEs_nogrey),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix2,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#Set heatmap3 color bar parameters
breakBarColors=c(-200,seq(-0.75, 0.75, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red") # create spectrum of blue to red for heatmaps

pdf("WGCNA Output/R_colorbar.pdf", width=2,height=5,pointsize = 10)
color.bar(barColors,-0.75,0.75,nticks=7)
dev.off()

png("WGCNA Output/traitHeatmap2_justCor.png",height=7,width=10,units="in",res=1000)
heatmap3(t(justCor), Colv=as.dendrogram(hclust_out),
         scale = "none",col=barColors, breaks=breakBarColors,margins=c(5,15),
         legendfun=function()showLegend(legend=c(NA),col=c(NA)),cexRow = 1.2,cexCol = 1.3,
         ColSideWidth=1,RowSideLabs=NA,
         labCol =  MEcolors[,3][match(substring(rownames(justCor),3),MEcolors[,1])],
         ColSideColors = cbind(c(substring(rownames(justCor),3))),ColSideLabs = NA)
dev.off()

pdf("WGCNA Output/traitHeatmap2_justCor.pdf",height=7,width=8)
heatmap3(t(justCor), Colv=as.dendrogram(hclust_out),
         scale = "none",col=barColors, breaks=breakBarColors,margins=c(5,15),
         legendfun=function()showLegend(legend=c(NA),col=c(NA)),cexRow = 1.2,cexCol = 1.3,
         ColSideWidth=1,RowSideLabs=NA,
         labCol =  MEcolors[,3][match(substring(rownames(justCor),3),MEcolors[,1])],
         ColSideColors = cbind(c(substring(rownames(justCor),3))),ColSideLabs = NA)
dev.off()






# Hub Genes ~~~~~~~~~~~~~~~~~~~ ####

datExpr_genenames = datExpr
colnames(datExpr_genenames) = gene.df[,2][match(colnames(datExpr),gene.df[,1])]
kME = signedKME(datExpr = datExpr_genenames,datME = MEs)

ME9_AD=import("Data/MagentaAD_MAGMA.txt",header=FALSE)

kME_ME9_AD = kME[which(toupper(rownames(kME)) %in% unlist(ME9_AD)),]
kME_ME9_AD$genes = rownames(kME_ME9_AD)
kME_ME9_AD = kME_ME9_AD[,c(16,14)] %>%
  arrange(kMEmagenta)

ME1_AD=import("Data/TurquoiseAD_MAGMA.txt",header=FALSE)

kME_ME1_AD = kME[which(toupper(rownames(kME)) %in% unlist(ME1_AD)),]
kME_ME1_AD$genes = rownames(kME_ME1_AD)
kME_ME1_AD = kME_ME1_AD[,c(16,6)] %>%
  arrange(kMEturquoise)

hubgenes = list()


for(col_i in 1:{ncol(MEs)-1}){
  ME=substring(colnames(kME)[col_i],first = 4)
  kME_thisME = kME[which(net$colors==ME),]
  
  sorted = sort(kME_thisME[,col_i],decreasing = TRUE,index.return=TRUE)
  top_kME = kME_thisME[sorted$ix,][1:25,]
  this_ME = data.frame(gene=rownames(top_kME),kME=top_kME[,col_i])
  hubgenes$this_ME = this_ME
  names(hubgenes)[which(names(hubgenes)=="this_ME")] = ME
}

hub_ul = unlist(hubgenes)
just_genes=hub_ul[which(str_detect(names(hub_ul),pattern = "gene"))] %>% as.data.frame()
export(just_genes,"WGCNA Output/hubgenes.csv",row.names=TRUE)



# kME Filter Gene List for GO

datExpr_genenames = datExpr
colnames(datExpr_genenames) = gene.df[,2][match(colnames(datExpr),gene.df[,1])]
kME = signedKME(datExpr = datExpr_genenames,datME = MEs)

for(col_i in 1:{ncol(MEs)-1}){
  ME=substring(colnames(kME)[col_i],first = 4)
  kME_thisME = kME[which(net$colors==ME),]
  kME_filt=kME_thisME[which(kME_thisME[,which(colnames(kME)==paste0("kME",ME))]>0.60),]
  go_genes=rownames(kME_filt) %>% as.data.frame()
  export(go_genes,paste0("WGCNA Output/GO Gene Lists/",ME,".txt"),col.names=FALSE)
}






# Cell Type FET ~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# downloaded CellTypeFET repository from edammer: https://github.com/edammer/CellTypeFET.git


# Cross-species FET modified to optionally adjust for symbol lookup inefficiency/loss (Wrapper R Script)
# by Eric Dammer & Divya Nandakumar
#=================================================================================================================#
# As published in: 
# Seyfried et al Cell Syst, 2017	https://www.sciencedirect.com/science/article/pii/S2405471216303702
# Johnson ECB et al, Nat Med, 2020	https://www.nature.com/articles/s41591-020-0815-6
# Johnson ECB et al, Nat Neurosci, 2022	https://www.nature.com/articles/s41593-021-00999-y
# ...and others
#=================================================================================================================#
# Provides and charts statistics for significance of hypergeometric overlap of gene symbol or 'UniqueID' lists.
# Can be used for cell type marker enrichment in modules, or any gene list overlap.
#
# UniqueIDs are formatted as "Symbol|ID...", where symbol is the official gene symbol (case is species specific).
#
# Sample marker list data provided are for 5 cell types of mammalian brain
# derived from thresholded filtering of supplemental data in acutely
# isolated purified cells from mouse brain measured as mRNA, in:
#  Ye Zhang et al, J Neurosci, 2014
#  https://www.jneurosci.org/content/34/36/11929.short
#  (Barres RNA app now at https://www.brainrnaseq.org/ )
#
# or measured as protein, in:
#  Kirti Sharma et al, Nat Neurosci, 2015 https://www.nature.com/articles/nn.4160
#
# Sample gene clusters are provided for single cell cluster markers in mRNA
# from fly brain, as published in:
#  Kristofer Davie et al, Cell, 2018
#  https://www.sciencedirect.com/science/article/pii/S0092867418307207
#
# Sample WGCNA modules representing proteome organization of the motor cortex
# of the brains of ALS and healthy individuals is provided for testing use only.
# Please reach out if you are interested in repurposing this unpublished data.
#=================================================================================================================#
# edammer@emory.edu Eric Dammer, Bioinformatic Scientist, NT Seyfried Systems Biology Group
# Emory University School of Medicine (2023)
#=================================================================================================================#

options(stringsAsFactors=FALSE)

# CellType FET PREPARATION OF IN-MEMORY VARIABLES (optional) ####

## Load into memory 3 variables for your network
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/CellTypeFET-main"))
# load("motorCtx-ALS_data_forORA.Rdata")  # sample Rdata file provided to demonstrate the script using motor cortex brain total homogenate in an unpublished ALS cohort
# DO NOT RELEASE OR USE WITHOUT PERMISSION


numericMeta = metadata

# ## Standardize the variable names of for modulesInMemory option below
# net<-netALS  		# list output of WGCNA::blockwiseModules() function building your coexpression network. The one necessary item in the list is the colors variable, i.e. net$colors or net[["colors"]]
# cleanDat<-cleanDatALS	# data.frame with rows representing gene products measured and all making it into the WGCNA network as goodGenes.
# # Rownames are expected to have the official gene symbol followed by a pipe and any other information after, i.e. "NEFL|..."
# numericMeta<-numericMetaALS   #sample traits/metadata

## Standardize the rownames of cleanDat here in the form "Symbol|...additional info", if needed.
#betterRownames<-read.csv(file="FlyNet_cleanDat6467_betterNames.csv",row.names=1,header=TRUE)
#rownames(cleanDat)<-betterRownames$UniqueID


# CellType FET CONFIGURATION PARAMETERS FULL LIST ####
##             WITH SAMPLE VALUES GIVEN FOR SAMPLE DATA PROVIDED                   ##


refDataFiles <- c("MyGene-Mouse-SharmaZhangUnion.csv")	# Originally mined brain cell type mRNA and protein lists were based on experiments in MOUSE.

# Use Modules in memory OR a .csv file with your input gene lists.
modulesInMemory=TRUE                              		# Load modules as categories? (If TRUE, categoriesFile not used, but you need cleanDat, net[["colors"]] and numericMeta variables)
# categoriesFile="Fly_Seurat_87cluster_BrainCellTypes.csv"	# File Name of Categories (Lists of Fly genes), only loaded if modulesInMemory=FALSE
# NOTE this file format has a column for official gene symbols of each module or cluster, with the cluster name/ID as column names in row 1
# categorySpeciesCode="mmusculus"				# What species are the gene symbols in categoriesFile?

# Other Options
allowDuplicates=TRUE				# Allow duplicate symbols across different lists for overlap?
# (should be true if you have general cell type lists and e.g. disease-associated phenotype cell type lists)
resortListsDecreasingSize=FALSE			# resort categories/modules and reference data lists? (decreasing size order)
barOption=FALSE					# draw bar charts for each list overlap instead of a heatmap.
adjustFETforLookupEfficiency=FALSE		# adjust p FET input for cross-species lookup inefficiency/loss of list member counts?
verticalCompression=3				# Plot(s) are squeezed into 1 row out of this many in each PDF page, compressing the heatmap tracks vertically (or the bar chart heights) for each reference list)
reproduceHistoricCalc=FALSE			# should be FALSE unless trying to reproduce exact calculations of prior publications listed.
#####################################################################################


## Generate Sample Outputs

# Load Seyfried/Emory pipeline FET as function geneListFET() having all the parameters described above, many with defaults used.
source("./geneListFET.R")


## output enrichment significance of provided test network data in memory as -log10(FDR) heatmap; enrichments checked are in provided 5 brain cell type marker gene lists
geneListFET(FileBaseName="1_ES",
            heatmapTitle=" Cell Type Marker Enrichment",
            modulesInMemory=TRUE,categorySpeciesCode="mmusculus",  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles,speciesCode="mmusculus",refDataDescription="NeuronVsNot")  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?


## Same as above, but output a bar chart instead of a heatmap.
geneListFET(FileBaseName="2_ES_bars", barOption=TRUE,
            heatmapTitle="Cell Type Marker Enrichment",
            modulesInMemory=TRUE,categorySpeciesCode="mmusculus",  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles,speciesCode="mmusculus",refDataDescription="NeuronVsNot")  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?
#NOTE: 2 cell type enrichments per row (with verticalCompression # of rows per PDF page); sample output #2 has 5 bar charts for the ALS motor cortex network's modules' gene enrichment into each of the human of 5 cell type marker lists.
#      and the next pages has the same network's modules' gene enrichment in the mouse symbol lists of the same 5 cell type markers.


# ## output enrichment significance of provided "categories" (or fly brain gene clusters) input from a .csv file instead of using coexpression modules in memory; -log10(FDR) heatmap; enrichments are checked against provided 5 mammalian brain cell type marker gene lists
# geneListFET(FileBaseName="3.FlyBrainSCgene87SeuratClusters_FET_to_5mammalianBrainCellTypes",
#             heatmapTitle="Fly Brain 87 Cluster Overlaps with 5 Mammalian Brain Cell Type Marker Reference Lists",
#             modulesInMemory=FALSE,categoriesFile="Fly_Seurat_87cluster_BrainCellTypes.csv",categorySpeciesCode="dmelanogaster",  #use clusters of gene symbols in a provided file, categoriesFile; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
#             refDataFiles=refDataFiles,speciesCode="mmusculus",refDataDescription="5brainCellTypes")  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?


CTE_in = rio::import("2_Enrichment Bar Charts_no1x.Overlap.in.MyGene-Mouse-SharmaZhangUnion.csv-DuplicatesALLOWED-hitListStats.csv")
FDR_CTE = CTE_in[9:13,2:ncol(MEs)] %>% sapply(as.numeric)
rownames(FDR_CTE) = CTE_in[9:13,1]
colnames(FDR_CTE) = CTE_in[1,2:ncol(MEs)]
log_FDR_CTE = -log(FDR_CTE,base=10)

for(i in 1:nrow(log_FDR_CTE)){
  CT = rownames(log_FDR_CTE)[i]
  CT_ES = sort(log_FDR_CTE[i,],decreasing = TRUE)
  pdf(paste0("../WGCNA Output/Cell Type Enrichment/",CT,".pdf"),height=3,width=5.5)
  barplot(CT_ES,col=names(CT_ES), ylab="-log(FDR p-value)",
          names.arg = MEcolors[,3][match(names(CT_ES),MEcolors[,1])],las=2,main=CT)+
    abline(h=1,lty=2, lwd=1.5, col="red")+
    calibrate::textxy(X=ncol(MEs)-7,Y=1,paste0("FDR p-value = 0.01"),cex=1,col="red")
  dev.off()
  
  png(paste0("../WGCNA Output/Cell Type Enrichment/",CT,".png"),height=3,width=5.5,res=600,units="in")
  barplot(CT_ES,col=names(CT_ES), ylab="-log(FDR p-value)",
          names.arg = MEcolors[,3][match(names(CT_ES),MEcolors[,1])],las=2,main=CT)+
    abline(h=1,lty=2, lwd=1.5, col="red")+
    calibrate::textxy(X=ncol(MEs)-7,Y=1,paste0("FDR p-value = 0.01"),cex=1,col="red")
  dev.off()
}
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # tells R that the folder where this code resides is our working directory


################ MAGMA ################

setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/MAGMA.SPA-main"))

#MAGMA-SPA (Seyfried Pipeline Adaptation for MAGMA)
#---------------------------------

# Required parameters, variables, and data must be set as shown above in .GlobalEnv before calling function; currently no defaults are automatic.
##################################
MAGMAinputDir= "MAGMAinput"

MAGMAinputs= c(	"AD_GWAS_ENSEMBLE_averageMinusLogP(Plt0.05).csv",
                "ALS_GWAS_ENSEMBLE_avgMinusLogP(Plt0.05).csv")     #These files must be in MAGMAinputDir

maxP=0.05                 #no genes with a MAGMA summarized p value greater than this will be considered even if in the MAGMA-derived input files.
FDR=0.10                  #FDR or q value (0 < FDR < 1); recommend 0.10, i.e. 10%
barcolors= c("darkslateblue","mediumorchid")  #specify one unique color for each of above MAGMAinputs
#common colors: "darkslateblue","mediumorchid", "seagreen3","hotpink","goldenrod","darkorange","darkmagenta", ...
relatednessOrderBar=TRUE  #Plot mean scaled enrichment bar plot in column order (relatedness) of MEs?  If FALSE, they will be plotted in size rank order M1, M2, ...

# Data created during the Seyfried Analysis Pipeline
##################################


# NETcolors= netTBI$colors     #module color assignments, vector of length equal to number of rows in cleanDat; should have all colors for modules from 1:minimumSizeRank as printed by WGCNA::labels2colors(1:nModules)
# cleanDat #rownames must start with HUMAN gene symbols, separated by any other rowname information using ';' or '|' character
rownames(cleanDat) = str_to_upper(rownames(cleanDat))

# Other variables
#################################
outFilePrefix="5"         #Filename prefix; step in the pipeline -- for file sorting by name.
outFileSuffix="AD_ALS_GWAS_ensembleAvg"
parallelThreads=8         #Each permutation analysis is run on a separate thread simultaneously, up to this many threads.
calculateMEs=TRUE        #Recalculate MEs and their relatedness order, even if the data already exists.
plotOnly=FALSE            #If plotOnly is TRUE, the variables created by MAGMA.SPA function holding plot data should already exist (xlabels, allBarData).
##################################

# Run the permutation analysis and generate all outputs
source("MAGMA.SPA.R")
MAGMAoutList <- MAGMA.SPA()
# Outputs XLSX, PDF, and list of barplot y values (allBarData), barplot labels (xlabels), and all permutation statistics and gene symbol hits (all_output)


# Rerun function, just to plot previously calculated statistics
allBarData<-MAGMAoutList$allBarData
xlabels<-MAGMAoutList$xlabels
all_output<-MAGMAoutList$all_output
plotOnly=TRUE
MAGMAoutList <- MAGMA.SPA()


FDR=0.025

magmaOut = read_xlsx("5.MAGMA-SPA-AD_ALS_GWAS_ensembleAvg.xlsx",sheet = "ALS_GWAS_ENSEMBLE_avgMinusLogP(")
scaledES = magmaOut[2,-1] #%>% as.numeric() %>% sort(decreasing = TRUE)
sortInd = scaledES %>% as.numeric() %>% sort(decreasing=TRUE,index.return=TRUE)


pdf(paste0("../WGCNA Output/MAGMA GWAS Enrichment/ALS Enrichment Scores.pdf"),height=3,width=5.5)
barplot(sortInd$x,col = colnames(scaledES)[sortInd$ix],las=2,
        main="MAGMA Enrichment: ALS",ylab="E.S.",
        names.arg = MEcolors[,3][match(colnames(scaledES)[sortInd$ix],MEcolors[,1])]) +
  abline(h=qnorm(FDR,lower=F), lty=2, lwd=1.5, col="red")+
  calibrate::textxy(X=ncol(MEs)-4,Y=qnorm(FDR,lower=F),paste0("FDR = ",FDR*100,"%"),cex=1,col="red")
dev.off()

png(paste0("../WGCNA Output/MAGMA GWAS Enrichment/ALS Enrichment Scores.png"),height=3,width=5.5,res=600,units="in")
barplot(sortInd$x,col = colnames(scaledES)[sortInd$ix],las=2,
        main="MAGMA Enrichment: ALS",ylab="E.S.",
        names.arg = MEcolors[,3][match(colnames(scaledES)[sortInd$ix],MEcolors[,1])]) +
  abline(h=qnorm(FDR,lower=F), lty=2, lwd=1.5, col="red")+
  calibrate::textxy(X=ncol(MEs)-4,Y=qnorm(FDR,lower=F),paste0("FDR = ",FDR*100,"%"),cex=1,col="red")
dev.off()



magmaOut = read_xlsx("5.MAGMA-SPA-AD_ALS_GWAS_ensembleAvg.xlsx",sheet = "AD_GWAS_ENSEMBLE_averageMinusLo")
scaledES = magmaOut[2,-1] #%>% as.numeric() %>% sort(decreasing = TRUE)
sortInd = scaledES %>% as.numeric() %>% sort(decreasing=TRUE,index.return=TRUE)

pdf(paste0("../WGCNA Output/MAGMA GWAS Enrichment/AD Enrichment Scores.pdf"),height=3,width=5.5)
barplot(sortInd$x,col = colnames(scaledES)[sortInd$ix],las=2,
        main="MAGMA Enrichment: AD",ylab="E.S.",
        names.arg = MEcolors[,3][match(colnames(scaledES)[sortInd$ix],MEcolors[,1])]) +
  abline(h=qnorm(FDR,lower=F), lty=2, lwd=1.5, col="red")+
  calibrate::textxy(X=ncol(MEs)-4,Y=qnorm(FDR,lower=F),paste0("FDR = ",FDR*100,"%"),cex=1,col="red")
dev.off()

png(paste0("../WGCNA Output/MAGMA GWAS Enrichment/AD Enrichment Scores.png"),height=3,width=5.5,res=600,units="in")
barplot(sortInd$x,col = colnames(scaledES)[sortInd$ix],las=2,
        main="MAGMA Enrichment: AD",ylab="E.S.",
        names.arg = MEcolors[,3][match(colnames(scaledES)[sortInd$ix],MEcolors[,1])]) +
  abline(h=qnorm(FDR,lower=F), lty=2, lwd=1.5, col="red")+
  calibrate::textxy(X=ncol(MEs)-4,Y=qnorm(FDR,lower=F),paste0("FDR = ",FDR*100,"%"),cex=1,col="red")
dev.off()


# What's going on with ME1 ###############

load("R Data/geneVSD.rda")
load("R Data/proteins.rda")

turq_genes = gene.df[which(net$colors=="turquoise"),2]
nge = datExpr_genenames %>% t()
nge = nge[which(rownames(nge) %in% turq_genes),]

dataZ = apply(t(nge),2,"scale")
hc=hclust(dist(t(dataZ), method = "euclidean"), method = "average")
Colv = as.dendrogram(hc)
rownames(dataZ) = rownames(t(nge))
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers, used for z-scored data
barColors = gplots::colorpanel(length(breakBarColors)-1, "blue", "white", "red2")

heatmap3(dataZ, 
         col=barColors, breaks=breakBarColors,legendfun=function()heatmap3::showLegend(legend=c(NA),col=c(NA),cex=2.5),
         Rowv=NA,Colv=Colv,
         # labRow = labRow, labCol = labCol,
         scale="none", cexCol=1,cexRow=1, margins=c(5,2),
         main = "", xlab = "", ylab = "",RowSideColors = cbind(TimeColors,InjColors))

heatmap_wl_sidebar(df_samples[,analyte_range],RowSideColors = cbind(TimeColors,InjColors),clust_c=TRUE)


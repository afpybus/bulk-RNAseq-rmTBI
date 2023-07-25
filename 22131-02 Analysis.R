# File:    22131-02 Analysis.R
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
               gplots,GSVA,ropls,RColorBrewer,DESeq2,gridExtra,WGCNA) # install/load these packages
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
library(tidyverse)
library(heatmap3)


# Load in RDA objects produced by Wrangle.R

load("R Data/Subset 0 0.5 4 24 48 48.5 52 72 96 96.5 100 120.rda")

output_folder = "Figures/"
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)




########### GSVA MDSigDB C2 #####################
# 
# breakBarColors=c(-200,seq(-2, 2, 0.01),200) #Outside numbers clip outliers, used for z-scored data
# barColors = gplots::colorpanel(length(breakBarColors)-1, "blue", "white", "red2")
# 
# # GSVA Clustered Heatmap
# dir.create(paste0(output_folder,"GSVA MSigDB/"),recursive = TRUE,showWarnings = FALSE)
# png(paste0(output_folder,"GSVA MSigDB/GSVA MSigDB Clustered Heatmap.png"), width=7,height=8,pointsize = 14, units="in",res=600)
# heatmap3(gseZ.C2, ColSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,  
#          col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
#          Rowv=as.dendrogram(hrGSVA.C2), 
#          #Rowv=NA,
#          Colv=NA,  scale="none",
#          margins = c(3,10),
#          labRow = NA, labCol=NA,
#          cexCol=1.3,
#          main = "GSVA MSigDB") 
# dev.off()

# 
# if(pdf){pdf(paste0(output_folder,"GSVA MSigDB/GSVA MSigDB Clustered Heatmap_Clust2.pdf"), width=7,height=8,pointsize = 14, useDingbats = FALSE)}else(
#   png(paste0(output_folder,"GSVA MSigDB/GSVA MSigDB Clustered Heatmap_Clust2.png"), width=7,height=8,pointsize = 14, units="in",res=600))
# heatmap3(gseZ.C2[which(cutree(hrGSVA,k=3)==2),], ColSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,  
#          col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
#          #Rowv=as.dendrogram(hrGSVA), 
#          #Rowv=NA,
#          Colv=NA,  scale="none",
#          margins = c(3,10),
#          labRow = NA, labCol=NA,
#          cexCol=1.3)
# title("GSVA Cluster", adj = 0.5, line = 2)
# dev.off()
# 
# clust2 = which(cutree(hrGSVA,k=3)==2) %>%
#   names() %>%
#   as.matrix()
# 
# export(clust2,paste0(output_folder,"GSVA MSigDB/GSVA MSigDB clust2.csv"))

sets_of_interest = c("KEGG_ALZHEIMERS_DISEASE",
                     "KEGG_PARKINSONS_DISEASE",
                     "REACTOME_AMYLOID_FIBER_FORMATION",
                     "REACTOME_INTERLEUKIN_1_SIGNALING",
                     "REACTOME_CELLULAR_RESPONSES_TO_STRESS",
                     "REACTOME_MAPK6_MAPK4_SIGNALING")

for (j in c("AbsoluteTime","Injury","TimePoint")){
  dir.create(paste0(output_folder,"GSVA MSigDB/",j),recursive=TRUE,showWarnings = FALSE)
  x = dplyr::select(metadata,all_of(j)) %>% as.matrix()
  if(j=="Injury"){x=metadata$Injury}
  for (i in 1:length(sets_of_interest)){
    set_ind = which(rownames(gse.C2)==sets_of_interest[i])
    p = error_plot(x=x,y=gseZ.C2[set_ind,],color=metadata$Injury,title = sets_of_interest[i],
                   ylab = paste0("Expression [a.u.]")) %>%
      ggpar(palette=unique(InjColorBar))
    saveRDS(p,paste0(output_folder,"GSVA MSigDB/",j,"/Mean SEM ",sets_of_interest[i],".rds"))
    pdf(paste0(output_folder,"GSVA MSigDB/",j,"/Mean SEM ",sets_of_interest[i],".pdf"),height=6,width=7.8,useDingbats = FALSE)
    print(p)
    dev.off()
    png(paste0(output_folder,"GSVA MSigDB/",j,"/Mean SEM ",sets_of_interest[i],".png"),height=6,width=7.8,units="in",res=600)
    print(p)
    dev.off()
  }
}

for(j in c("AbsoluteTime","Injury","TimePoint")){
  grid.plots = list()
  for(i in 1:length(sets_of_interest)){
    grid.plots$new = readRDS(paste0(output_folder,"GSVA MSigDB/",j,"/Mean SEM ",sets_of_interest[i],".rds"))
    names(grid.plots)[which(names(grid.plots)=="new")] = sets_of_interest[i]
  }
  nrow=floor(sqrt(length(sets_of_interest)))
  ncol=ceiling(length(sets_of_interest)/nrow)
  
  
  dir.create(paste0(output_folder,"Summary Panels/MSigDB GSVA"),recursive = TRUE,showWarnings = FALSE)
  png(paste0(output_folder,"GSVA MSigDB/",j,"/MSigDB GSVA Sets by ",j," Panel.png"),res=600,units="in",height=6*nrow+0.2,width=7.8*ncol)
  grid.arrange(grobs=grid.plots,nrow=nrow,top=text_grob(paste0("Gene Set Expression by ",j,"\n"),size = 24,face = "bold"))
  dev.off()
  png(paste0(output_folder,"Summary Panels/MSigDB GSVA/MSigDB GSVA Sets by ",j," Panel.png"),res=600,units="in",height=6*nrow+0.2,width=7.8*ncol)
  grid.arrange(grobs=grid.plots,nrow=nrow,top=text_grob(paste0("Gene Set Expression by ",j,"\n"),size = 24,face = "bold"))
  dev.off()
}

metadata$Injury2 = as.character(metadata$Injury)
metadata$Injury2[metadata$Injury=="pre-3xCHI"] = "2xCHI"
metadata$Injury2[metadata$Injury=="pre-5xCHI"] = "4xCHI"
metadata$Injury2=factor(metadata$Injury2,levels=c("Sham","1xCHI","2xCHI","3xCHI","4xCHI","5xCHI"))

alz_ind=which(str_detect(rownames(gse.C2),"KEGG_ALZHEIMERS_DISEASE"))
alz = gseZ.C2[alz_ind,]
pdf("Figures/GSVA MSigDB/Pathology Figure/Kegg Alz.pdf",height=4.5,width=5)
error_plot(x=metadata$Injury2,y=alz,color = metadata$Injury,ylab="Expression [a.u.]",title = "Kegg Alzheimer's Disease") +
  stat_compare_means(method="t.test",ref.group = "Sham",label="p.signif") +
  scale_color_manual(values=unique(InjColorBar)) +
  theme(legend.position = "none")
dev.off()

par_ind=which(str_detect(rownames(gse.C2),"KEGG_PARKINSONS_DISEASE"))
par = gseZ.C2[par_ind,]
pdf("Figures/GSVA MSigDB/Pathology Figure/Kegg PD.pdf",height=4.5,width=5)
error_plot(x=metadata$Injury2,y=par,color = metadata$Injury,ylab="Expression [a.u.]",title = "Kegg Parkinson's Disease") +
  stat_compare_means(method="t.test",ref.group = "Sham",label="p.signif") +
  scale_color_manual(values=unique(InjColorBar)) +
  theme(legend.position = "none")
dev.off()

stats_df = data.frame(alz=alz,Inj=metadata$Injury2)
compare_means(alz~Inj,data=stats_df,ref.group = "Sham",method="t.test",p.adjust.method = "bonferroni")


########### GSVA Astrocyte Sets #####################

breakBarColors=c(-200,seq(-2, 2, 0.01),200) #Outside numbers clip outliers, used for z-scored data
barColors = gplots::colorpanel(length(breakBarColors)-1, "blue", "white", "red2")

# GSVA Clustered Heatmap
dir.create(paste0(output_folder,"GSVA Astro/"),recursive = TRUE,showWarnings = FALSE)
pdf(paste0(output_folder,"GSVA Astro/GSVA Astro Clustered Heatmap.pdf"), width=7,height=8,pointsize = 14, useDingbats = FALSE)
# png(paste0(output_folder,"GSVA Astro/GSVA Astro Clustered Heatmap.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(gseZ.astro, ColSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA.astro), 
         #Rowv=NA,
         Colv=NA,  scale="none",
         margins = c(3,10),
         #labRow = NA, 
         labCol=NA,
         cexCol=1.3)
title("GSVA Astro Cluster", adj = 0.5, line = 2)
dev.off()


# GSVA Clustered Heatmap
pdf(paste0(output_folder,"GSVA Astro/GSVA Astro Clustered Heatmap2.pdf"), width=8,height=7,pointsize = 14, useDingbats = FALSE)
# png(paste0(output_folder,"GSVA Astro/GSVA Astro Clustered Heatmap.png"), width=8,height=7,pointsize = 14, units="in",res=600)
heatmap3(t(gseZ.astro), RowSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Colv=as.dendrogram(hrGSVA.astro),
         Rowv=NA,
         scale="none",
         margins = c(10,10),
         labRow = NA,
         # labCol=NA,
         cexCol=1.3)
title("GSVA Astro Cluster", adj = 0.5, line = 2)
dev.off()

# GSVA Clustered Heatmap
dir.create(paste0(output_folder,"Summary Panels/Astro GSVA"),recursive = TRUE,showWarnings = FALSE)
png(paste0(output_folder,"Summary Panels/Astro GSVA/Astro GSVA Clustered Heatmap.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(gseZ.astro, ColSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA.astro), 
         #Rowv=NA,
         Colv=NA,  scale="none",
         margins = c(3,10),
         #labRow = NA, 
         labCol=NA,
         cexCol=1.3)
title("GSVA Astro Cluster", adj = 0.5, line = 2)
dev.off()

metadata=metadata.sub
for (j in c("AbsoluteTime","Injury","TimePoint")){  # c("AbsoluteTime","Injury","TimePoint")
  dir.create(paste0(output_folder,"GSVA Astro/",j),recursive=TRUE,showWarnings = FALSE)
  x = dplyr::select(metadata,all_of(j)) %>% as.matrix()
  if(j=="Injury"){x=metadata$Injury}
  for (i in 1:dim(gse.astro)[1]){
    p = error_plot(x=x,y=gseZ.astro[i,],color=metadata$Injury,title = rownames(gseZ.astro)[i],
                   ylab= paste0(rownames(gseZ.astro)[i]," [a.u.]")) %>%
      ggpar(palette=unique(InjColorBar))
    saveRDS(p,paste0(output_folder,"GSVA Astro/",j,"/Mean SEM ",rownames(gseZ.astro)[i],".rds"))
    pdf(paste0(output_folder,"GSVA Astro/",j,"/Mean SEM ",rownames(gseZ.astro)[i],".pdf"),height=6,width=7.8,useDingbats = FALSE)
    print(p)
    dev.off()
    png(paste0(output_folder,"GSVA Astro/",j,"/Mean SEM ",rownames(gseZ.astro)[i],".png"),height=6,width=7.8,units="in",res=600)
    print(p)
    dev.off()
    
    for(gene in c("Gfap","S100b")){
      dir.create(paste0(output_folder,"GSVA Astro/",j,"/",gene),showWarnings = FALSE)
      gene.y.log = norm.genes.sub[which(gene.df$geneName==gene),] %>% as.matrix() %>% log() %>% scale()
      p = error_plot_gradient(x=x,y=gseZ.astro[i,],color=gene.y.log,title = rownames(gseZ.astro)[i],
                              ylab = paste0(rownames(gseZ.astro)[i]," [a.u.]"),color.str=paste0("z-score of\nlog(",gene,")"))
      saveRDS(p,paste0(output_folder,"GSVA Astro/",j,"/",gene,"/Mean SEM ",rownames(gseZ.astro)[i],".rds"))
      pdf(paste0(output_folder,"GSVA Astro/",j,"/",gene,"/Mean SEM ",rownames(gseZ.astro)[i],".pdf"),height=6,width=7.8,useDingbats = FALSE)
      print(p)
      dev.off()
      png(paste0(output_folder,"GSVA Astro/",j,"/",gene,"/Mean SEM ",rownames(gseZ.astro)[i],".png"),height=6,width=7.8,units="in",res=600)
      print(p)
      dev.off()
    }
  }
}

cutree.astro.out = cutree(hrGSVA.astro,k=2)
sets=list()
for(i in 1:length(unique(cutree.astro.out))){
  sets$new = names(which(cutree.astro.out==i))
  names(sets)[which(names(sets)=="new")] = paste0("set",i)
}

for(set_ind in 1:length(sets)){
  sets_of_interest=sets[[set_ind]]
  for(j in c("AbsoluteTime","Injury","TimePoint")){
    
    # Normal panels colored by Injury group
    grid.plots = list()
    for(i in 1:length(sets_of_interest)){
      grid.plots$new = readRDS(paste0(output_folder,"GSVA Astro/",j,"/Mean SEM ",sets_of_interest[i],".rds"))
      names(grid.plots)[which(names(grid.plots)=="new")] = sets_of_interest[i]
    }
    
    nrow=floor(sqrt(length(sets_of_interest)))
    ncol=ceiling(length(sets_of_interest)/nrow)
    png(paste0(output_folder,"GSVA Astro/",j,"/Cluster Set ",set_ind," Panel.png"),res=600,units="in",height=6*nrow+0.2,width=7.8*ncol)
    grid.arrange(grobs=grid.plots,nrow=nrow,top=text_grob(paste0("Astrocyte Gene Set Expression in Cluster ",set_ind," by ",j,"\n"),size = 24,face = "bold"))
    dev.off()
    dir.create(paste0(output_folder,"Summary Panels/Astro GSVA"),recursive = TRUE,showWarnings = FALSE)
    png(paste0(output_folder,"Summary Panels/Astro GSVA/Astro GSVA Cluster ",set_ind," by ",j," Panel.png"),res=600,units="in",height=6*nrow+0.2,width=7.8*ncol)
    grid.arrange(grobs=grid.plots,nrow=nrow,top=text_grob(paste0("Astrocyte Gene Set Expression in Cluster ",set_ind," by ",j,"\n"),size = 24,face = "bold"))
    dev.off()
    
    # Panels colored by gene expression
    for(gene in c("Gfap","S100b")){
      dir.create(paste0(output_folder,"Summary Panels/Astro GSVA/",gene),showWarnings = FALSE)
      grid.plots = list()
      for(i in 1:length(sets_of_interest)){
        grid.plots$new = readRDS(paste0(output_folder,"GSVA Astro/",j,"/",gene,"/Mean SEM ",sets_of_interest[i],".rds"))
        names(grid.plots)[which(names(grid.plots)=="new")] = sets_of_interest[i]
      }
      
      nrow=floor(sqrt(length(sets_of_interest)))
      ncol=ceiling(length(sets_of_interest)/nrow)
      png(paste0(output_folder,"GSVA Astro/",j,"/",gene,"/Cluster Set ",set_ind," Panel.png"),res=600,units="in",height=6*nrow+0.2,width=7.8*ncol)
      grid.arrange(grobs=grid.plots,nrow=nrow,top=text_grob(paste0("Astrocyte Gene Set Expression in Cluster ",set_ind," by ",j," and ",gene,"\n"),size = 24,face = "bold"))
      dev.off()
      png(paste0(output_folder,"Summary Panels/Astro GSVA/",gene,"/Astro GSVA Cluster ",set_ind," by ",j," Panel.png"),res=600,units="in",height=6*nrow+0.2,width=7.8*ncol)
      grid.arrange(grobs=grid.plots,nrow=nrow,top=text_grob(paste0("Astrocyte Gene Set Expression in Cluster ",set_ind," by ",j," and ",gene,"\n"),size = 24,face = "bold"))
      dev.off()
    }
  }
}

# Astro GSVA MLR ~~~~~~~~~~~
mlr.data = gse.astro  # columns are samples in same order as rows in metadata, rows are genes/gene sets/proteins

inj_num = case_when(
  metadata$Injury == "Sham" ~0,
  metadata$Injury == "1xCHI" ~1,
  metadata$Injury == "pre-3xCHI" ~2,
  metadata$Injury == "3xCHI" ~3,
  metadata$Injury == "4xCHI" ~4,
  metadata$Injury == "5xCHI" ~5,
)
inj_binary = case_when(
  metadata$Injury == "Sham" ~"Sham",
  metadata$Injury != "Sham" ~"Injured"
)
inj_tp = factor(metadata$TimePoint)
mlr.out = lm(t(mlr.data)~inj_num + inj_binary+ inj_tp)
rownames(mlr.out$coefficients) = c("Intercept","Number of Injuries", "Injured vs Sham","30min Time Point","4hr Time Point","24 hr Time Point")
lm_summary = summary(mlr.out)
MLR_pmat = matrix(nrow=dim(mlr.data)[1],ncol=mlr.out$rank)
rownames(MLR_pmat) = rownames(mlr.data)
colnames(MLR_pmat) = rownames(mlr.out$coefficients)
for (i in 1:nrow(MLR_pmat)){
  MLR_pmat[i,] = lm_summary[[i]]$coefficients[,4]}
sig_analytes = which(t(MLR_pmat)<0.05,arr.ind=TRUE) %>%
  data.frame() %>%
  mutate(color="red") %>%
  mutate(lwd=1)

# FDR p-value heatmap
breakBarColors=c(0,10^seq(-3, 0, 0.05)) 
barColors = colorpanel(length(breakBarColors)-1, "blue","white")
hc.mlr = hclust(dist(MLR_pmat, method = "euclidean"), method = "average")

dir.create(paste0(output_folder,"GSVA Astro/"),recursive = TRUE,showWarnings = FALSE)
pdf(paste0(output_folder,"GSVA Astro/GSVA Astro MLR heatmap.pdf"),useDingbats=F,height=4,width=8)
heatmap3(t(MLR_pmat),
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
         Rowv=NA, highlightCell = sig_analytes, Colv=as.dendrogram(hc.mlr),
         scale="none",margins=c(12,2),
         cexCol=0.5,cexRow=0.5,
         main = "Astro MLR p-values")
dev.off()

dir.create(paste0(output_folder,"Summary Panels/Astro GSVA"),recursive = TRUE,showWarnings = FALSE)
png(paste0(output_folder,"Summary Panels/Astro GSVA/GSVA Astro MLR heatmap.png"),res=600,units="in",height=4,width=8)
heatmap3(t(MLR_pmat),
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
         Rowv=NA, highlightCell = sig_analytes, Colv=as.dendrogram(hc.mlr),
         scale="none",margins=c(12,2),
         cexCol=0.5,cexRow=0.5,
         main = "Astro MLR p-values")
dev.off()

# p-value log color bar
png(paste0(output_folder,"GSVA Astro/Clusterbar_p.png"), width=2,height=5,units="in",res=600,pointsize = 14)
color.bar(barColors,-3,0,nticks=4)
dev.off()

df_astro = gse.astro %>%
  t() %>%
  as_tibble() %>%
  cbind(metadata) %>%
  mutate(Injury2 = case_when(
    Injury == "pre-3xCHI" ~ "3xCHI",
    Injury == "pre-5xCHI" ~ "5xCHI",
    TRUE ~ Injury
  ))

colnames(df_astro)[which(colnames(df_astro)=="Stress Response")] = "StressResponse"
colnames(df_astro)[which(colnames(df_astro)=="Carbohydrate Metabolism")] = "CarbohydrateMetabolism"
colnames(df_astro)[which(colnames(df_astro)=="General Metabolism")] = "GeneralMetabolism"

compare_means(data=df_astro,formula=`Neural Development`~`TimePoint`,method="t.test",
              ref.group="0",p.adjust.method="bonferroni")

scatter_reg(x=proteins.sub$GFAP,y=df_astro$A1)
scatter_reg(x=proteins.sub$GFAP,y=df_astro$StressResponse)






####################################
# 
# df_inj = df_samples %>%
#   filter(AbsoluteTime %in% 120)
# 
# scatter_reg(x=df_inj$`L hemi 3mm`,y=df_inj$pTau,x.str = "CBF [a.u.]",y.str="Phospho-tau T181 [a.u.]",text_position = "top right")
# 
# df_sham = df_samples %>%
#   filter(AbsoluteTime == 0)
# 
# scatter_reg(x=df_sham$GFAP,y=df_sham$`L hemi 3mm`,x.str = "GFAP",y.str="CBF L 3mm")
# 
# 
# # # GSVA LOESS Time Curves ############
# 
# sham_means = gse.astro[,which(metadata$Injury=="Sham")] %>%
#   apply(1,mean)
# 
# loess_df = gse.astro %>%
#   sweep(1,sham_means,"/")
# 
# ggplot(data=metadata,aes(x=AbsoluteTime,y=gse.astro[1,]))+
#   geom_smooth(method = "loess")
# 
# ggline(metadata,x="AbsoluteTime",y=gse.astro[1,])
# gse.astro[1,]



########### GSVA Neuron Sets #####################

breakBarColors=c(-200,seq(-2, 2, 0.01),200) #Outside numbers clip outliers, used for z-scored data
barColors = gplots::colorpanel(length(breakBarColors)-1, "blue", "white", "red2")

# p-value log color bar
pdf(paste0(output_folder,"GSVA Astro/zbar_2.pdf"), width=2,height=5,pointsize = 14)
color.bar(barColors,-2,2,nticks=5)
dev.off()

# GSVA Clustered Heatmap
dir.create(paste0(output_folder,"GSVA Neuron/"),recursive = TRUE,showWarnings = FALSE)
pdf(paste0(output_folder,"GSVA Neuron/GSVA Neuron Clustered Heatmap.pdf"), width=7,height=8,pointsize = 14, useDingbats = FALSE)
# png(paste0(output_folder,"GSVA Neuron/GSVA Neuron Clustered Heatmap.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(gseZ.neuron, ColSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA.neuron), 
         #Rowv=NA,
         Colv=NA,  scale="none",
         margins = c(3,10),
         #labRow = NA, 
         labCol=NA,
         cexCol=1.3)
title("GSVA Neuron Cluster", adj = 0.5, line = 2)
dev.off()

# heatmap3(neuro_comb, ColSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,
#          col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
#          #Rowv=as.dendrogram(hrGSVA.neuron),
#          #Rowv=NA,
#          Colv=NA,  scale="none",
#          margins = c(3,10),
#          #labRow = NA,
#          labCol=NA,
#          cexCol=1.3)

pdf(paste0(output_folder,"GSVA Neuron/GSVA Neuron Clustered Heatmap2.pdf"), width=8,height=7,pointsize = 14, useDingbats = FALSE)
# png(paste0(output_folder,"GSVA Neuron/GSVA Neuron Clustered Heatmap2.png"), width=8,height=7,pointsize = 14, units="in",res=600)
heatmap3(t(gseZ.neuron), RowSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Colv=as.dendrogram(hrGSVA.neuron), 
         Rowv=NA,
         scale="none",
         margins = c(10,10),
         labRow = NA,
         # labCol=NA,
         cexCol=1.3)
title("GSVA Neuron Cluster", adj = 0.5, line = 2)
dev.off()

pdf(paste0(output_folder,"GSVA Neuron/GSVA Neuron Clustered Heatmap.pdf"), width=7,height=8,pointsize = 14, useDingbats = FALSE)
# png(paste0(output_folder,"GSVA Neuron/GSVA Neuron Clustered Heatmap.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(gseZ.neuron, ColSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA.neuron), 
         #Rowv=NA,
         Colv=NA,  scale="none",
         margins = c(3,10),
         #labRow = NA, 
         labCol=NA,
         cexCol=1.3)
title("GSVA Neuron Cluster", adj = 0.5, line = 2)
dev.off()

dir.create(paste0(output_folder,"Summary Panels/Neuron GSVA"),recursive = TRUE,showWarnings = FALSE)
png(paste0(output_folder,"Summary Panels/Neuron GSVA/Neuron GSVA Clustered Heatmap.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(gseZ.neuron, ColSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA.neuron), 
         #Rowv=NA,
         Colv=NA,  scale="none",
         margins = c(3,10),
         #labRow = NA, 
         labCol=NA,
         cexCol=1.3)
title("GSVA Neuron Cluster", adj = 0.5, line = 2)
dev.off()

for (j in c("AbsoluteTime","Injury","TimePoint")){
  dir.create(paste0(output_folder,"GSVA Neuron/",j),recursive=TRUE,showWarnings = FALSE)
  x = dplyr::select(metadata,all_of(j)) %>% as.matrix()
  if(j=="Injury"){x=metadata$Injury}
  for (i in 1:dim(gse.neuron)[1]){
    p = error_plot(x=x,y=gseZ.neuron[i,],color=metadata$Injury,title = rownames(gse.neuron)[i],
                   ylab = paste0(rownames(gse.neuron)[i]," [a.u.]")) %>%
      ggpar(palette=unique(InjColorBar))
    saveRDS(p,paste0(output_folder,"GSVA Neuron/",j,"/Mean SEM ",rownames(gse.neuron)[i],".rds"))
    pdf(paste0(output_folder,"GSVA Neuron/",j,"/Mean SEM ",rownames(gse.neuron)[i],".pdf"),height=6,width=7.8,useDingbats = FALSE)
    print(p)
    dev.off()
    png(paste0(output_folder,"GSVA Neuron/",j,"/Mean SEM ",rownames(gse.neuron)[i],".png"),height=6,width=7.8,units="in",res=600)
    print(p)
    dev.off()
  }
}

df_gseZ.neuron = gseZ.neuron %>%
  t() %>%
  as_tibble() %>%
  mutate(Sample = as.numeric(colnames(gseZ.neuron))) %>%
  left_join(metadata) %>%
  mutate(Injury2 = case_when(
    Injury %in% c("Sham","1xCHI") ~ "1xCHI",
    Injury %in% c("pre-3xCHI","3xCHI") ~ "3xCHI",
    Injury %in% c("pre-5xCHI","5xCHI") ~ "5xCHI"
  ))

for(i in 1:length(rownames(gse.neuron))){
  set = rownames(gse.neuron)[i]
  pdf(paste0(output_folder,"GSVA Neuron/Mean SEM by Inj Number ",set,".pdf"),height=3.5,width=8,useDingbats = FALSE)
  # png(paste0(output_folder,"GSVA Neuron/Mean SEM by Inj Number ",set,".png"),height=3.5,width=8,units="in",res=1000)
  print(ggerrorplot(df_gseZ.neuron,x="TimePoint",y=set,add=c("mean_se","jitter"),color="Injury",facet.by = "Injury2",
                    title = paste0(set," by Injury Number"),ylab = paste0(set," [a.u.]"),xlab="Time Point [hr]") %>%
          ggpar(palette=unique(InjColorBar),legend="none") + theme(plot.title=element_text(hjust=0.5,face = "bold")))
  dev.off()
}

cutree.neuron.out = cutree(hrGSVA.neuron,k=3)
sets=list()
for(i in 1:length(unique(cutree.neuron.out))){
  sets$new = names(which(cutree.neuron.out==i))
  names(sets)[which(names(sets)=="new")] = paste0("set",i)
}

for(set_ind in 1:length(sets)){
  sets_of_interest=sets[[set_ind]]
  for(j in c("AbsoluteTime","Injury","TimePoint")){
    grid.plots = list()
    for(i in 1:length(sets_of_interest)){
      grid.plots$new = readRDS(paste0(output_folder,"GSVA Neuron/",j,"/Mean SEM ",sets_of_interest[i],".rds"))
      names(grid.plots)[which(names(grid.plots)=="new")] = sets_of_interest[i]
    }
    nrow=floor(sqrt(length(sets_of_interest)))
    ncol=ceiling(length(sets_of_interest)/nrow)
    
    png(paste0(output_folder,"GSVA Neuron/",j,"/Cluster Set ",set_ind," Panel.png"),res=600,units="in",height=6*nrow+0.2,width=7.8*ncol)
    grid.arrange(grobs=grid.plots,nrow=nrow,top=text_grob(paste0("Neuron Gene Set Expression in Cluster ",set_ind," by ",j,"\n"),size = 24,face = "bold"))
    dev.off()
    dir.create(paste0(output_folder,"Summary Panels/Neuron GSVA"),recursive = TRUE,showWarnings = FALSE)
    png(paste0(output_folder,"Summary Panels/Neuron GSVA/Neuron GSVA Cluster ",set_ind," by ",j," Panel.png"),res=600,units="in",height=6*nrow+0.2,width=7.8*ncol)
    grid.arrange(grobs=grid.plots,nrow=nrow,top=text_grob(paste0("Neuron Gene Set Expression in Cluster ",set_ind," by ",j,"\n"),size = 24,face = "bold"))
    dev.off()
  }
}


# Neuron GSVA MLR ~~~~~~~~~~~
mlr.data = gse.neuron  # columns are samples in same order as rows in metadata, rows are genes/gene sets/proteins

inj_num = case_when(
  metadata$Injury == "Sham" ~0,
  metadata$Injury == "1xCHI" ~1,
  metadata$Injury == "pre-3xCHI" ~2,
  metadata$Injury == "3xCHI" ~3,
  metadata$Injury == "4xCHI" ~4,
  metadata$Injury == "5xCHI" ~5,
)
inj_binary = case_when(
  metadata$Injury == "Sham" ~"Sham",
  metadata$Injury != "Sham" ~"Injured"
)
inj_tp = factor(metadata$TimePoint)
mlr.out = lm(t(mlr.data)~inj_num + inj_binary+ inj_tp)
rownames(mlr.out$coefficients) = c("Intercept","Number of Injuries", "Injured vs Sham","30min Time Point","4hr Time Point","24 hr Time Point")
MLR_pmat = matrix(nrow=dim(mlr.data)[1],ncol=mlr.out$rank)
rownames(MLR_pmat) = rownames(mlr.data)
colnames(MLR_pmat) = rownames(mlr.out$coefficients)
for (i in 1:nrow(MLR_pmat)){
  MLR_pmat[i,] = lm_summary[[i]]$coefficients[,4]}
sig_analytes = which(t(MLR_pmat)<0.05,arr.ind=TRUE) %>%
  data.frame() %>%
  mutate(color="red") %>%
  mutate(lwd=1)

# FDR p-value heatmap
breakBarColors=c(0,10^seq(-3, 0, 0.05)) 
barColors = colorpanel(length(breakBarColors)-1, "blue","white")
hc.mlr = hclust(dist(MLR_pmat, method = "euclidean"), method = "average")

pdf(paste0(output_folder,"GSVA Neuron/GSVA Neuron MLR heatmap.pdf"),useDingbats=F,height=4,width=8)
heatmap3(t(MLR_pmat),
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
         Rowv=NA, highlightCell = sig_analytes, Colv=as.dendrogram(hc.mlr),
         scale="none",margins=c(12,2),
         cexCol=0.5,cexRow=0.5,
         main = "Neuron MLR p-values")
dev.off()

dir.create(paste0(output_folder,"Summary Panels/Neuron GSVA"),recursive = TRUE,showWarnings = FALSE)
png(paste0(output_folder,"Summary Panels/Neuron GSVA/GSVA Neuron MLR heatmap.png"),res=600,units="in",height=4,width=8)
heatmap3(t(MLR_pmat),
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
         Rowv=NA, highlightCell = sig_analytes, Colv=as.dendrogram(hc.mlr),
         scale="none",margins=c(12,2),
         cexCol=0.5,cexRow=0.5,
         main = "Neuron MLR p-values")
dev.off()

# p-value log color bar
png(paste0(output_folder,"GSVA Neuron/Clusterbar_p.png"), width=2,height=5,units="in",res=600,pointsize = 14)
color.bar(barColors,-3,0,nticks=4)
dev.off()







# # transgenes #########
# 
# trans_id=14960:14963
# transgenes=norm.genesZ[trans_ind,] %>%
#   t() %>%
#   as_tibble() %>%
#   mutate(Sample=as.numeric(colnames(norm.genesZ))) %>%
#   left_join(metadata,by="Sample")
# colnames(transgenes)[1:length(trans_id)] = gene.df$geneName[trans_id]
# 
# ggerrorplot(transgenes,x="AbsoluteTime",y="AC149090.1",add="jitter",color="Injury",palette = unique(InjColorBar)) +
#   stat_compare_means(ref.group = "0",method="t.test",hide.ns = TRUE,label = "p.signif")
# ggerrorplot(transgenes,x="AbsoluteTime",y="CAAA01118383.1",add="jitter",color="Injury",palette = unique(InjColorBar))+
#   stat_compare_means(ref.group = "0",method="t.test",hide.ns = TRUE,label = "p.signif")
# ggerrorplot(transgenes,x="AbsoluteTime",y="CAAA01147332.1",add="jitter",color="Injury",palette = unique(InjColorBar))+
#   stat_compare_means(ref.group = "0",method="t.test",hide.ns = TRUE,label = "p.signif")
# ggerrorplot(transgenes,x="AbsoluteTime",y="MAPT",add="jitter",color="Injury",palette = unique(InjColorBar))+
#   stat_compare_means(ref.group = "0",method="t.test",hide.ns = TRUE,label = "p.signif")


# CREATING SUBSET RDS DATA #################################
# 
# source("Sample Subset.R")
# load("R Data/genecounts.rda") # Transcript Data
# load("R Data/proteins.rda") # Protein Data
# load("R Data/gs.rda") # Gene Sets
# 
# # Computes the gene and metadata subset, GSVA, and clustering code for your filter then saves to RDS in Subsets folder
# subset_DoD(ATOI=c(0,0.5,4,24,48,48.5,52,72,96,96.5,100,120)) # All samples
# subset_DoD(ATOI=c(0,24,48,72,96,120)) # 24hr Time Points
# subset_DoD(ATOI=c(0,4,48,52,96,100)) # 4hr Time Points
# subset_DoD(ATOI=c(0,0.5,48,48.5,96,96.5)) # 30min Time Points


#### All Time Points #################



# # Check Housekeeping Genes ###########
# 
# hk = import("Reference Files/Mouse cell selective reference genes.csv")
# ng = norm.genes.sub
# rownames(ng) = gene.df[match(rownames(ng),gene.df[,1]),2]
# hk_ng = ng[which(rownames(ng) %in% hk$`Gene Symbol`),] %>%
#   t() %>%
#   scale()
# 
# heatmap_wl_sidebar(hk_ng,clust_r = FALSE,clust_c = TRUE,RowSideColors = cbind(InjColorBar,TimeColorBar))

# # tSNE of Genes ######################
# tsneOut = Rtsne(t(norm.genes.subZ),dims=2,perplexity=floor({ncol(norm.genes.subZ)-1}/3))
# color_by = c("AbsoluteTime","Injury","TimePoint","Round")
# dir.create(paste0(output_folder,"tSNE_Genes/"),recursive = TRUE,showWarnings = FALSE)
# dir.create(paste0(output_folder,"Summary Panels/"),recursive = TRUE,showWarnings = FALSE)
# grid.plots=list()
# for(colorcode in color_by){
#   p = scores_plot(tsneOut$Y[,1],tsneOut$Y[,2],color=as.matrix(dplyr::select(metadata.sub,all_of(colorcode))),analysis="tSNE",color.str = colorcode)
#   grid.plots$new = p
#   names(grid.plots)[which(names(grid.plots)=="new")] = colorcode
#   saveRDS(p,paste0(output_folder,"tSNE_Genes/ScorePlot_",colorcode,".rds"))
#   pdf(paste0(output_folder,"tSNE_Genes/ScorePlot_",colorcode,".pdf"),height=6,width=7.8,useDingbats = FALSE)
#   print(p)
#   dev.off()
#   png(paste0(output_folder,"tSNE_Genes/ScorePlot_",colorcode,".png"),height=6,width=7.8,units="in",res=600)
#   print(p)
#   dev.off()
# }
# png(paste0(output_folder,"tSNE_Genes/tSNE Genes Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("tSNE of All Genes\n"),size = 24,face = "bold"))
# dev.off()
# png(paste0(output_folder,"Summary Panels/tSNE Genes Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("tSNE of All Genes\n"),size = 24,face = "bold"))
# dev.off()

# # tSNE of Gene Sets MDSigDB C2 #######################
# tsneOut = Rtsne(t(gseZ.C2),dims=2,perplexity=floor({ncol(norm.genes.subZ)-1}/3))
# color_by = c("AbsoluteTime","Injury","TimePoint","Round")
# dir.create(paste0(output_folder,"GSVA MSigDB/tSNE GSVA MSigDB/"),recursive = TRUE,showWarnings = FALSE)
# dir.create(paste0(output_folder,"Summary Panels/"),recursive = TRUE,showWarnings = FALSE)
# grid.plots=list()
# for(colorcode in color_by){
#   p = scores_plot(tsneOut$Y[,1],tsneOut$Y[,2],color=as.matrix(dplyr::select(metadata.sub,all_of(colorcode))),analysis="tSNE",color.str = colorcode)
#   saveRDS(p,paste0(output_folder,"GSVA MSigDB/tSNE GSVA MSigDB/ScorePlot_",colorcode,".rds"))
#   grid.plots$new = p
#   names(grid.plots)[which(names(grid.plots)=="new")] = colorcode
#   pdf(paste0(output_folder,"GSVA MSigDB/tSNE GSVA MSigDB/ScorePlot_",colorcode,".pdf"),height=6,width=7.8,useDingbats = FALSE)
#   print(p)
#   dev.off()
#   png(paste0(output_folder,"GSVA MSigDB/tSNE GSVA MSigDB/ScorePlot_",colorcode,".png"),height=6,width=7.8,units="in",res=600)
#   print(p)
#   dev.off()
# }
# png(paste0(output_folder,"GSVA MSigDB/tSNE GSVA MSigDB/tSNE GSVA MSigDB Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("tSNE of MSigDB C2 Gene Sets\n"),size = 24,face = "bold"))
# dev.off()
# png(paste0(output_folder,"Summary Panels/tSNE GSVA MSigDB Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("tSNE of MSigDB C2 Gene Sets\n"),size = 24,face = "bold"))
# dev.off()

# # PCA of Genes ######################
# 
# oplsOut <- opls(x = t(norm.genes.subZ), predI=2)
# oplsVar = oplsOut@modelDF$R2X * 100
# rotateOut <- rotate_opls(oplsOut,degrees = 0)
# color_by = c("AbsoluteTime","Injury","TimePoint","Round")
# dir.create(paste0(output_folder,"PCA_Genes/"),recursive = TRUE,showWarnings = FALSE)
# dir.create(paste0(output_folder,"Summary Panels/"),recursive = TRUE,showWarnings = FALSE)
# grid.plots=list()
# for(colorcode in color_by){
#   p = scores_plot(rotateOut$T1,rotateOut$T2,color=as.matrix(dplyr::select(metadata.sub,all_of(colorcode))),analysis="PCA",color.str = colorcode)
#   grid.plots$new = p
#   names(grid.plots)[which(names(grid.plots)=="new")] = colorcode
#   saveRDS(p,paste0(output_folder,"PCA_Genes/ScorePlot_",colorcode,".rds"))
#   pdf(paste0(output_folder,"PCA_Genes/ScorePlot_",colorcode,".pdf"),height=6,width=7.8,useDingbats = FALSE)
#   print(p)
#   dev.off()
#   png(paste0(output_folder,"PCA_Genes/ScorePlot_",colorcode,".png"),height=6,width=7.8,units="in",res=600)
#   print(p)
#   dev.off()
# }
# png(paste0(output_folder,"PCA_Genes/PCA Genes Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("PCA of All Genes\n"),size = 24,face = "bold"))
# dev.off()
# png(paste0(output_folder,"Summary Panels/PCA Genes Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("PCA of All Genes\n"),size = 24,face = "bold"))
# dev.off()

# # PCA of Proteins ######################
# 
# oplsOut <- opls(x = proteins.sub[,analyte_range], predI=2)
# oplsVar = oplsOut@modelDF$R2X * 100
# rotateOut <- rotate_opls(oplsOut,degrees = 0)
# color_by = c("AbsoluteTime","Injury","TimePoint","Round")
# dir.create(paste0(output_folder,"PCA_Proteins/"),recursive = TRUE,showWarnings = FALSE)
# grid.plots=list()
# for(colorcode in color_by){
#   p = scores_plot(rotateOut$T1,rotateOut$T2,color=as.matrix(dplyr::select(metadata.sub,all_of(colorcode))),analysis="PCA",color.str = colorcode)
#   grid.plots$new = p
#   names(grid.plots)[which(names(grid.plots)=="new")] = colorcode
#   saveRDS(p,paste0(output_folder,"PCA_Proteins/ScorePlot_",colorcode,".rds"))
#   pdf(paste0(output_folder,"PCA_Proteins/ScorePlot_",colorcode,".pdf"),height=6,width=7.8,useDingbats = FALSE)
#   print(p)
#   dev.off()
#   png(paste0(output_folder,"PCA_Proteins/ScorePlot_",colorcode,".png"),height=6,width=7.8,units="in",res=600)
#   print(p)
#   dev.off()
# }
# png(paste0(output_folder,"PCA_Proteins/PCA Proteins Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("PCA of All Proteins\n"),size = 24,face = "bold"))
# dev.off()
# png(paste0(output_folder,"Summary Panels/PCA_Proteins Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("PCA of All Proteins\n"),size = 24,face = "bold"))
# dev.off()

# # PCA of Gene Sets MDSigDB ######################
# 
# oplsOut <- opls(x = t(gseZ.C2), predI=2)
# oplsVar = oplsOut@modelDF$R2X * 100
# rotateOut <- rotate_opls(oplsOut,degrees = 0)
# color_by = c("AbsoluteTime","Injury","TimePoint","Round")
# dir.create(paste0(output_folder,"GSVA MSigDB/PCA GSVA MSigDB/"),recursive = TRUE,showWarnings = FALSE)
# grid.plots=list()
# for(colorcode in color_by){
#   p = scores_plot(rotateOut$T1,rotateOut$T2,color=as.matrix(dplyr::select(metadata.sub,all_of(colorcode))),analysis="PCA",color.str = colorcode)
#   saveRDS(p,paste0(output_folder,"GSVA MSigDB/PCA GSVA MSigDB/ScorePlot_",colorcode,".rds"))
#   grid.plots$new = p
#   names(grid.plots)[which(names(grid.plots)=="new")] = colorcode
#   pdf(paste0(output_folder,"GSVA MSigDB/PCA GSVA MSigDB/ScorePlot_",colorcode,".pdf"),height=6,width=7.8,useDingbats = FALSE)
#   print(p)
#   dev.off()
#   png(paste0(output_folder,"GSVA MSigDB/PCA GSVA MSigDB/ScorePlot_",colorcode,".png"),height=6,width=7.8,units="in",res=600)
#   print(p)
#   dev.off()
# }
# png(paste0(output_folder,"GSVA MSigDB/PCA GSVA MSigDB/PCA GSVA MSigDB Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("PCA of MSigDB C2 Gene Sets\n"),size = 24,face = "bold"))
# dev.off()
# png(paste0(output_folder,"Summary Panels/PCA GSVA MSigDB Scores Plot Panel.png"),res=600,units="in",height=8,width=10)
# grid.arrange(grobs=grid.plots,ncol=2,top=text_grob(paste0("PCA of MSigDB C2 Gene Sets\n"),size = 24,face = "bold"))
# dev.off()

# # # GFAP vs S100B Analysis ##########
# gfap.gene = norm.genes[which(gene.df$geneName=="Gfap"),] %>% as.matrix() %>% log() %>% scale()
# s100b.gene = norm.genes[which(gene.df$geneName=="S100b"),] %>% as.matrix() %>% log() %>% scale()
# scatter_reg(x=gfap.gene,y=s100b.gene,x.str = "Gfap",y.str="S100b")
# inj=metadata$Injury
# dataPlot = data.frame(gfap.gene,s100b.gene,inj)
# ggscatter(data=dataPlot,x="gfap.gene",y="s100b.gene",color="inj") %>%
#   ggpar(palette=unique(InjColorBar))
# inj_binary = case_when(
#   metadata$AbsoluteTime<=96 ~ "Less than 5xCHI",
#   metadata$AbsoluteTime>96 ~ "5xCHI"
# )
# scatter_reg_group(x=gfap.gene,y=s100b.gene,x.str = "Gfap",y.str="S100b",grouping = inj_binary,text_position = "bottom right")
# 
# 
# gfap.protein = proteins$GFAP
# s100b.protein = proteins$S100B
# inj = metadata$Injury
# i=4
# new_meanSEM_gradient(x=x,y=gseZ.astro[i,],colors=gfap.protein,title = names(gs_astro)[i],
#                      y_str = paste0(names(gs_astro)[i]," [a.u.]"),legend.name=paste0("z-score"))
# 
# scatter_reg(x=gfap.gene,y=gfap.protein,x.str = "Gfap Gene Expression",y.str="GFAP Protein",text_position = "bottom right")
# dataPlot = data.frame(gfap.gene,gfap.protein,inj)
# ggscatter(data=dataPlot,x="gfap.gene",y="gfap.protein",color="inj") %>%
#   ggpar(palette=unique(InjColorBar))
# scatter_reg_group(x=gfap.gene,y=gfap.protein,x.str = "Gfap Gene",y.str="GFAP Protein",grouping = inj_binary,text_position = "bottom right")
# 
# scatter_reg(gfap.protein,s100b.gene)

# # Sort by mean difference between sham and 120hr ##########
# differenceVec=rowMeans(gseZ[,which(metadata.sub$AbsoluteTime==0)])-rowMeans(gseZ[,which(metadata.sub$AbsoluteTime==120)])
# diffSort=sort(differenceVec, index.return=TRUE)
# 
# # GSVA Mean Diff Sorted Heatmap
# if(pdf){pdf(paste0(output_folder,"GSVA Mean Diff Heatmap.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)} else(
#   png(paste0(output_folder,"GSVA Mean Diff Heatmap.png"), width=7,height=8,pointsize = 14,units="in",res=600))
# heatmap3(gseZ[diffSort$ix,], ColSideColors=cbind(TimePoint=c(TimeColorBar),Injury=c(InjColorBar)),  ColSideWidth=1, ColSideLabs=NA,    
#          col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
#          Rowv=NA,
#          Colv=NA,  scale="none",
#          margins = c(3,10),
#          labRow = NA, labCol=NA,
#          cexRow=0.85, cexCol=1.3)
# title("GSVA Diff Sort", adj = 0.5, line = 2)
# dev.off()









#### All 24hr Time Points #################
# 
# ATOI=c(0,24,48,72,96,120)
# settings=str_c(ATOI,collapse=" ")
# 
# #loads the already computed subset from RDS file saved to Data folder
# subset=readRDS(paste0("Subsets/Subset ",settings,".RDS"))
# norm.genes = subset$norm.genes.sub
# norm.genesZ = subset$norm.genes.subZ
# metadata = subset$metadata.sub
# gse = subset$gse
# gseZ = subset$gseZ
# hrGSVA = subset$hrGSVA
# 
# colorBars = makeColorBars(metadata)
# InjColorBar = colorBars$InjColorBar
# TimeColorBar = colorBars$TimeColorBar
# 
# output_folder = paste0(tag,"Figure Output/",settings,"/")
# dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)
# 
# 
# tsneOut = Rtsne(t(norm.genesZ),dims=2,perplexity=floor({ncol(norm.genesZ)-1}/3))
# 
# if(pdf){pdf(paste0(output_folder,"tSNE_Genes_ScorePlot.pdf"),height=6,width=7.8,useDingbats = FALSE)}else(
#   png(paste0(output_folder,"tSNE_Genes_ScorePlot.png"),height=6,width=7.8,units="in",res=600))
# scores.plot(tsneOut$Y[,1],tsneOut$Y[,2],color=metadata$AbsoluteTime,analysis="tSNE")
# dev.off()
# 
# 
# tsneOut = Rtsne(t(gseZ),dims=2,perplexity=floor({ncol(norm.genesZ)-1}/3))
# 
# if(pdf){pdf(paste0(output_folder,"tSNE_GSVA_ScorePlot.pdf"),height=6,width=7.8,useDingbats = FALSE)}else(
#   png(paste0(output_folder,"tSNE_GSVA_ScorePlot.png"),height=6,width=7.8,units="in",res=600))
# scores.plot(tsneOut$Y[,1],tsneOut$Y[,2],color=metadata$AbsoluteTime,analysis="tSNE")
# dev.off()
# 
# 
# # Sort by mean difference between sham and 120hr
# differenceVec=rowMeans(gseZ[,which(metadata.sub$AbsoluteTime==0)])-rowMeans(gseZ[,which(metadata.sub$AbsoluteTime==120)])
# diffSort=sort(differenceVec, index.return=TRUE)
# 
# # GSVA Mean Diff Sorted Heatmap
# if(pdf){pdf(paste0(output_folder,"GSVA Mean Diff Heatmap.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)} else(
#   png(paste0(output_folder,"GSVA Mean Diff Heatmap.png"), width=7,height=8,pointsize = 14,units="in",res=600))
# heatmap3(gseZ[diffSort$ix,], ColSideColors =cbind(InjColorBar),  ColSideWidth=1, ColSideLabs=NA,    
#          col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
#          Rowv=NA,
#          Colv=NA,  scale="none",
#          margins = c(3,10),
#          labRow = NA, labCol=NA,
#          cexRow=0.85, cexCol=1.3)
# title("GSVA Diff Sort", adj = 0.5, line = 2)
# dev.off()
# 
# # GSVA Clustered Heatmap
# if(pdf){pdf(paste0(output_folder,"GSVA Clustered Heatmap.pdf"), width=7,height=8,pointsize = 14, useDingbats = FALSE)}else(
#   png(paste0(output_folder,"GSVA Clustered Heatmap.png"), width=7,height=8,pointsize = 14, units="in",res=600))
# heatmap3(gseZ, ColSideColors =cbind(InjColorBar), ColSideWidth=1, ColSideLabs=NA,  
#          col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
#          Rowv=as.dendrogram(hrGSVA), 
#          #Rowv=NA,
#          Colv=NA,  scale="none",
#          margins = c(3,10),
#          labRow = NA, labCol=NA,
#          cexCol=1.3)
# title("GSVA Cluster", adj = 0.5, line = 2)
# dev.off()



# CBF Analysis ###############
# 
# tsneOut = Rtsne(t(norm.genesZ),dims=2,perplexity=floor({ncol(norm.genesZ)-1}/3))
# dir.create(paste0(output_folder,"CBF/"),recursive = TRUE,showWarnings = FALSE)
# p = scores.plot.gradient(tsneOut$Y[,1],tsneOut$Y[,2],color=as.matrix(select(proteins,"L hemi 3mm")),analysis="tSNE",color.str = "CBF")
# saveRDS(p,paste0(output_folder,"CBF/tSNE_Genes_ScorePlot_CBF.rds"))
# pdf(paste0(output_folder,"CBF/tSNE_Genes_ScorePlot_CBF.pdf"),height=6,width=7.8,useDingbats = FALSE)
# print(p)
# dev.off()
# png(paste0(output_folder,"CBF/tSNE_Genes_ScorePlot_CBF.png"),height=6,width=7.8,units="in",res=600)
# print(p)
# dev.off()
# 
# oplsOut <- opls(x = t(norm.genesZ), predI=2)
# oplsVar = oplsOut@modelDF$R2X * 100
# rotateOut <- rotate.opls(oplsOut,degrees = 0)
# dir.create(paste0(output_folder,"CBF/"),recursive = TRUE,showWarnings = FALSE)
# p = scores.plot.gradient(rotateOut$T1,rotateOut$T2,color=as.matrix(select(proteins,"L hemi 3mm")),analysis="PCA",color.str = "CBF")
# saveRDS(p,paste0(output_folder,"CBF/PCA_Genes_ScorePlot_CBF.rds"))
# pdf(paste0(output_folder,"CBF/PCA_Genes_ScorePlot_CBF.pdf"),height=6,width=7.8,useDingbats = FALSE)
# print(p)
# dev.off()
# png(paste0(output_folder,"CBF/PCA_Genes_ScorePlot_CBF.png"),height=6,width=7.8,units="in",res=600)
# print(p)
# dev.off()
# 
# 
# dir.create(paste0(output_folder,"CBF/"),recursive = TRUE,showWarnings = FALSE)
# gse.cbf.Z = rbind(gseZ,t(scale(proteins$`L hemi 3mm`)))
# rownames(gse.cbf.Z)[dim(gse.cbf.Z)[1]] = "CBF"
# corMat = rcorr(t(gse.cbf.Z))
# p.cbf = corMat$P[which(rownames(corMat$P)=="CBF"),] 
# p.cbf.fdr = p.cbf %>% p.adjust("fdr")
# R.cbf = corMat$r[which(rownames(corMat$P)=="CBF"),]
# cbf.correlations = cbind(p.cbf,p.cbf.fdr,R.cbf)[-length(p.cbf),]
# cbf.correlations.sort=cbf.correlations[sort(cbf.correlations[,1],index.return=TRUE)$ix,]
# cbf.correlations.sort = cbind(cbf.correlations.sort,rownames(cbf.correlations.sort))
# export(cbf.correlations.sort,paste0(output_folder,"CBF/CBF Correlating Pathways.csv"))
# 
# pathways_of_interest_cbf = c(
#   "REACTOME_INTERLEUKIN_1_PROCESSING",
#   "REACTOME_INTERLEUKIN_18_SIGNALING",
#   "BIOCARTA_CYTOKINE_PATHWAY",
#   "BIOCARTA_NFKB_PATHWAY",
#   "REACTOME_RIP_MEDIATED_NFKB_ACTIVATION_VIA_ZBP1",
#   "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY",
#   "BIOCARTA_PCAF_PATHWAY"
# )
# 
# cbf = proteins$`L hemi 3mm`
# inj=metadata$Injury
# for(i in 1:length(pathways_of_interest_cbf)){
#   pathway = gse[which(rownames(gse)==pathways_of_interest_cbf[i]),]
#   dataPlot = data.frame(cbf,pathway,inj)
#   p1 = new_meanSEM(x=metadata$AbsoluteTime,y=pathway,colors=inj,title = pathways_of_interest_cbf[i]) %>% ggpar(palette = unique(InjColorBar))
#   p2 = new_meanSEM_gradient(x=metadata$AbsoluteTime,y=pathway,colors=cbf,title = pathways_of_interest_cbf[i],legend.name = "CBF")
#   p3 = scatter_reg(x=cbf,y=pathway,
#                    x.str = "CBF [a.u.]",y.str = "GSVA Enrichment Score",title=pathways_of_interest_cbf[i])
#   p4 = ggscatter(data=dataPlot,x="cbf",y="pathway",color="inj",size = 5,xlab = "CBF [a.u.]",ylab = "Enrichment Score",title = pathways_of_interest_cbf[i]) %>%
#     ggpar(palette=unique(InjColorBar))
#   dir.create(paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i]),recursive = TRUE,showWarnings = FALSE)
#   saveRDS(p1,paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/AbsTimeByInj.rds"))
#   saveRDS(p2,paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/AbsTimebyCBF.rds"))
#   saveRDS(p3,paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/RegressionCBF.rds"))
#   saveRDS(p4,paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/ScatterByInj.rds"))
#   png(paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/AbsTimeByInj.png"),width = 8,height=4.5,units = "in",res=600); print(p1); dev.off()
#   png(paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/AbsTimebyCBF.png"),width = 8,height=4.5,units = "in",res=600); print(p2); dev.off()
#   png(paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/RegressionCBF.png"),width = 5,height=4.5,units = "in",res=600); print(p3); dev.off()
#   png(paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/ScatterByInj.png"),width = 5,height=4.5,units = "in",res=600); print(p4); dev.off()
#   pdf(paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/AbsTimeByInj.pdf"),width = 8,height=4.5); print(p1); dev.off()
#   pdf(paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/AbsTimebyCBF.pdf"),width = 8,height=4.5); print(p2); dev.off()
#   pdf(paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/RegressionCBF.pdf"),width = 5,height=4.5); print(p3); dev.off()
#   pdf(paste0(output_folder,"CBF/GSVA Pathways/",pathways_of_interest_cbf[i],"/ScatterByInj.pdf"),width = 5,height=4.5); print(p4); dev.off()
# }
# # gse.cbf.Z = rbind(gseZ,t(scale(proteins$`L hemi 3mm`)))
# # rownames(gse.cbf.Z)[dim(gse.cbf.Z)[1]] = "CBF"
# # hrGSVA.cbf= hclust(dist((gse.cbf.Z),method = "euclidean"), method="ward.D2")
# # cutree.cbf.out = cutree(hrGSVA.cbf,k=15)
# # cutree.cbf.out[which(names(cutree.cbf.out)=="CBF")]
# # cbf.clust = which(cutree.cbf.out == cutree.cbf.out[which(names(cutree.cbf.out)=="CBF")]) %>%
# #   names() %>%
# #   as.matrix()
# # export(cbf.clust,paste0(output_folder,"CBF Co-Varying Pathways.csv"))
# # 







# CIBERSORTx #######################
# 
# # Cell-Type Specific Gene Expression Profiles
# 
# GEPs = import("Reference Files/CIBERSORTx_Job10_output_ExpressionProfiles/CIBERSORTxGEP_Job10_GEPs_Filtered.txt")
# 
# L5IT = select(GEPs,GeneSymbol,"L5 IT")
# L5IT_noNA = L5IT[!is.na(L5IT$`L5 IT`),]
# 
# 
# neuro_GSVA_genes = unlist(gs_neuron)
# genes_OI = neuro_GSVA_genes[which(neuro_GSVA_genes %in% str_to_upper(L5IT_noNA$GeneSymbol))]
# sig_OI = L5IT_noNA %>%
#   filter(str_to_upper(GeneSymbol) %in% genes_OI) %>%
#   dplyr::arrange(-`L5 IT`)
# 
# ggbarplot(sig_OI,x="GeneSymbol",y="L5 IT",add="mean_se")
# 
# ggboxplot(sig_OI,y="L5 IT")

####################################################




#CIBERSTORTx

# CS_out = import("Reference Files/CIBERSORTx_Job9_Results_SBatch_Corrected.csv") %>%
#   dplyr::rename(Sample = Mixture)
# CS_out = CS_out[match(metadata$Sample,CS_out$Sample),] %>%
#   select(all_of(2:18)) %>%
#   as.matrix()
# rownames(CS_out) = metadata$Sample
# 
# cs_t = CS_out %>%
#   apply(2,"scale") %>%
#   t()
# neuro_comb = rbind(gseZ.neuron,cs_t)
# 
# glut = c("L2/3 IT","L5 ET","L5 IT","L5/6 NP", "L6 CT","L6 IT","L6b")
# gaba = c("Lamp5","Meis2","Pvalb","Sncg","Sst","Sst Chodl","Vip")
# nonn = c("Astro","Endo","Micro-PVM","Oligo","OPC","Peri","VLMC")
# 
# 
# bars = colMeans(CS_out) %>% sort(decreasing = TRUE) * 100
# bars_col = case_when(
#   names(bars) %in% glut ~ "palegreen2",
#   names(bars) %in% gaba ~ "orchid",
#   names(bars) %in% nonn ~ "deepskyblue"
# )
# pdf("CIBERSORTx Avg Percent of Each Cell Type.pdf",height=4,width=5)
# barplot(bars,col = bars_col,las=2,ylab = "Average Percent of Sample",legend.text = c("Glutamatergic","Non-Neuronal","GABAergic"))
# dev.off()
# 
# CS_tib = CS_out %>%
#   t() %>%
#   as_tibble() %>%
#   mutate(CellType = colnames(CS_out)) %>%
#   mutate(Class = case_when(
#     CellType %in% glut ~ "Glutamatergic",
#     CellType %in% gaba ~ "GABAergic",
#     CellType %in% nonn ~ "Non-Neuronal"
#   )) %>%
#   gather(key="Sample",value = "Percent",1:39) %>%
#   mutate(CellType = factor(CellType,levels=names(bars)))
# CS_tib$Percent = CS_tib$Percent * 100
# 
# # pdf("CIBERSORTx Avg Percent of Each Cell Type with error bar.pdf",height=4,width=6)
# png("CIBERSORTx Avg Percent of Each Cell Type with error bar.png",height=4,width=6,units="in",res=600)
# ggbarplot(CS_tib,x="CellType",y="Percent",fill="Class",add=c("mean_se"),ylab = "Percent of Sample",xlab="",title="CIBERSORTx - Average Cell Type Composition") +
#   rotate_x_text() +
#   theme(plot.title = element_text(hjust=0.5,face="bold"),legend.position = "right",legend.title = element_blank(),legend.text = element_text(size=10))
# dev.off()




AD_ind = which(toupper(gene.df$geneName) %in% unlist(AD_genes))
AD_gene_ncz = data.frame(norm.genesZ)[AD_ind,]
rownames(AD_gene_ncz) = gene.df$geneName[AD_ind]

InjPalette = rev(RColorBrewer::brewer.pal(6,"Spectral"))
InjPalettePar = c("Sham"=InjPalette[1],"1xCHI"=InjPalette[2],"pre-3xCHI"=InjPalette[3],"3xCHI"=InjPalette[4],"pre-5xCHI"=InjPalette[5],"5xCHI"=InjPalette[6])
TimePalette = RColorBrewer::brewer.pal(4,"Purples")
TimePalettePar = c("0"=TimePalette[1],"0.5"=TimePalette[2],"4"=TimePalette[3],"24"=TimePalette[4])

InjColors = case_when(
  metadata$Injury=="Sham" ~ InjPalette[1],
  metadata$Injury== "1xCHI" ~ InjPalette[2],
  metadata$Injury== "pre-3xCHI" ~ InjPalette[3],
  metadata$Injury== "3xCHI" ~ InjPalette[4],
  metadata$Injury== "pre-5xCHI" ~ InjPalette[5],
  metadata$Injury== "5xCHI" ~ InjPalette[6])

TimeColors = case_when(
  metadata$TimePoint== 0 ~ TimePalette[1],
  metadata$TimePoint== 0.5 ~ TimePalette[2],
  metadata$TimePoint== 4 ~ TimePalette[3],
  metadata$TimePoint== 24 ~ TimePalette[4],)


# With Injury Bar and Time Bar
cairo_pdf(paste0(output_folder,"/AD_Kegg_hm.pdf"),width=20,height=7)
heatmap_wl_sidebar(data=t(AD_gene_ncz),title="",labRow = NA,cex_c=0.3,mar = c(3,.2),
                   RowSideColors = cbind(Injury=InjColors,TimePoint=TimeColors),clust_c = TRUE)  
dev.off()

# With Sex Bar, Injury Bar, and Time Bar
if(pdf){cairo_pdf(paste0(SOI,"/Heatmaps/Heatmap_SexInjTime_",ROI,".pdf"),width=7,height=12)} else
{png(paste0(SOI,"/Heatmaps/Heatmap_SexInjTime_",ROI,".png"),units="in",width=7,height=12,res=600)}
Heatmap_SideCol(hmMat,cbind(Sex=SexColors,Injury=InjColors,TimePoint=TimeColors),title="",cex_r = 0.3,cex_c=0.7,mar = c(2,8))
dev.off()

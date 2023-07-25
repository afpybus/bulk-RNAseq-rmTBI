# Create function to subset the data ######################
subset_DoD <- function(TOI=c(0,0.5,4,24),
                       ATOI=c(0,0.5,4,24,48,48.5,52,72,96,96.5,100,120),
                       IOI=c("Sham","1xCHI","pre-3xCHI","3xCHI","pre-5xCHI","5xCHI")) 
{
  # Selecting groups by creating an index of all samples that match the filter
  subset.index <- which(metadata$TimePoint %in% TOI & metadata$AbsoluteTime %in% ATOI & metadata$Injury %in% IOI)
  
  # Create subsets of genes and metadata for analysis
  count.sub <- count[,subset.index]     # Matrix of subset data
  metadata.sub <- metadata[subset.index,]       # Tibble of subset metadata
  if(identical(metadata$Sample,df_samples$Sample)){proteins.sub = df_samples[subset.index,]}else{stop("Protein and Metadata not indexed together.")}
  
  # DESeq2 Normalization  
  
  des_mat = count.sub %>%
    round(digits=0) %>%
    as.matrix()
  dds <- DESeqDataSetFromMatrix(countData = des_mat,
                                colData = metadata.sub,
                                design = ~ AbsoluteTime)
  dds <- DESeq(dds)
  dds <- estimateSizeFactors(dds)
  norm.genes.sub <- counts(dds, normalized=TRUE)
  
  # Create a tibble of z-scored subset data
  norm.genes.subZ <- norm.genes.sub %>%
    apply(1,"scale") %>%
    t() %>%
    as_tibble()
  colnames(norm.genes.subZ) = colnames(norm.genes.sub)
  rownames(norm.genes.subZ) = rownames(norm.genes.sub)
  norm.genes.subZ[is.na(norm.genes.subZ)] <- 0   # set any NA values to 0
  
  ng.gse = norm.genes.sub %>% as.matrix()
  rownames(ng.gse) = str_to_upper(gene.df$geneName)
  
  ng.gse_mouse = ng.gse
  rownames(ng.gse_mouse) = gene.df$geneName
  
  ##### GSVA #########
  
  # MSigDB C2
  gse.C2=gsva(ng.gse, gs_MSigDB_C2, mx.diff=FALSE)
  gseZ.C2=(apply(gse.C2,1,scale)) %>%
    t()  #z-score the data
  colnames(gseZ.C2)=colnames(gse.C2)
  gseZ.C2[is.nan(gseZ.C2)] = 0
  hrGSVA.C2= hclust(dist((gseZ.C2),method = "euclidean"), method="ward.D2")
  
  # MSigDB M2
  gse.M2=gsva(ng.gse_mouse, gs_MSigDB_M2, mx.diff=FALSE)
  gseZ.M2=(apply(gse.M2,1,scale)) %>%
    t()  #z-score the data
  colnames(gseZ.M2)=colnames(gse.M2)
  gseZ.M2[is.nan(gseZ.M2)] = 0
  hrGSVA.M2= hclust(dist((gseZ.M2),method = "euclidean"), method="ward.D2")
  
  # Astrocyte Annotations
  gse.astro=gsva(ng.gse, gs_astro, mx.diff=FALSE)
  gseZ.astro=(apply(gse.astro,1,scale)) %>%
    t()  #z-score the data
  colnames(gseZ.astro)=colnames(gse.astro)
  gseZ.astro[is.nan(gseZ.astro)] = 0
  hrGSVA.astro= hclust(dist((gseZ.astro),method = "euclidean"), method="ward.D2")
  
  # Neuron Annotations
  gse.neuron=gsva(ng.gse, gs_neuron, mx.diff=FALSE)
  gseZ.neuron=(apply(gse.neuron,1,scale)) %>%
    t()  #z-score the data
  colnames(gseZ.neuron)=colnames(gse.neuron)
  gseZ.neuron[is.nan(gseZ.neuron)] = 0
  hrGSVA.neuron= hclust(dist((gseZ.neuron),method = "euclidean"), method="ward.D2")
  
  
  ##### Set Color Palettes #########
  
  InjPalette = rev(RColorBrewer::brewer.pal(6,"Spectral"))
  InjColorBar <- case_when(
    metadata.sub$Injury == "Sham" ~ InjPalette[1],
    metadata.sub$Injury == "1xCHI" ~ InjPalette[2],
    metadata.sub$Injury == "pre-3xCHI" ~ InjPalette[3],
    metadata.sub$Injury == "3xCHI" ~ InjPalette[4],
    metadata.sub$Injury == "pre-5xCHI" ~ InjPalette[5],
    metadata.sub$Injury == "5xCHI" ~ InjPalette[6],
  )
  
  TimePalette = RColorBrewer::brewer.pal(4,"Purples")
  TimeColorBar <- case_when(
    metadata.sub$TimePoint == 0 ~ TimePalette[1],
    metadata.sub$TimePoint == 0.5 ~ TimePalette[2],
    metadata.sub$TimePoint == 4 ~ TimePalette[3],
    metadata.sub$TimePoint == 24 ~ TimePalette[4]
  )
  
  ####### Save Output ###############
  
  #creates a string to describe filter setting
  settings=paste(str_c(unique(metadata.sub$AbsoluteTime),collapse=" "))
  
  norm.genes = norm.genes.sub
  norm.genesZ = norm.genes.subZ
  count = count.sub
  metadata = metadata.sub
  proteins = proteins.sub
  
  save(norm.genes,norm.genesZ,count,metadata,proteins,gse.C2,gseZ.C2,hrGSVA.C2,gse.M2,gseZ.M2,hrGSVA.M2,
       gse.astro,gse.neuron,gseZ.astro,gseZ.neuron,hrGSVA.astro,hrGSVA.neuron,gene.df,
       InjColorBar,TimeColorBar,analyte_range,file=paste0("R Data/Subset ",settings,".rda"))
}
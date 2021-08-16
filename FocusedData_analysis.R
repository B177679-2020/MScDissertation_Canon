# PRELIMINARY ANALYSIS OF FOCUSED DATASET
# Load in focused dataset from INCISE project, link it to metadata from Batch1, create a targets df
# Preliminary overview of dataset with proportions of polyp type, recurrence...
# PCA analysis.


##Load libraries 
library(limma)
library(scatterplot3d)
library(stringr)
library(ggplot2)
library(edgeR)
library(tidyr)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnrichmentBrowser")
library(EnrichmentBrowser)

BiocManager::install("pathfindR")
library(EnrichmentBrowser)

load_metadata_df <- function(Metadata){
  #Rename columns of the dataframes so that they are named with INCISE code
  rownames(Metadata) <- Metadata[,1]
  #Metadata <- Metadata[,-1]
  return(Metadata)
}

load_data_df <- function(FocusedData, Metadata) {

  
  names(FocusedData)[1] = "Gene_ID"
  rownames(FocusedData) <- FocusedData[,1]
  #FocusedData <- FocusedData[,-1]
  
  ##Processing of sample names. TempOSeq data has names like "CD38_4647". There might be several probes for the
  #same gene with different numeric part. This lines of code aim to separate the gene names from the numeric code
  #and group the repeats adding up their counts. Also ensure the index column is called Gene_ID.
  FocusedData <- separate(data = FocusedData, col = Gene_ID, into = c("Gene_ID", "tag"), sep = "_")
  FocusedData <- FocusedData %>% group_by(Gene_ID) %>% summarise(across(where(is.numeric), sum))
  FocusedData <- as.data.frame(FocusedData)
  rownames(FocusedData) <- FocusedData[,1]
  FocusedData <- FocusedData[,-1]

  
  ##Remove rows that have all 0 (genes not expressed in any sample)
  #FocusedData <- FocusedData[rowSums(FocusedData[])>10,]
  #Quantile10 = rowQuantiles(as.matrix(FocusedData2), probs = 0.1)
  keep.exprs <- filterByExpr(FocusedData, min.count=3)
  FocusedData = FocusedData[keep.exprs,]

  ##Join Data with Metadata
  #Identify metadata samples which are present in the Batch dataset
  intersection_f <- substr(colnames(FocusedData), start = 1, stop = 7) %in% rownames(Metadata)
  #Keep only data that has metadata associated
  FocusedData <- FocusedData[,c(intersection_f)]
  
  return(FocusedData)

  }

create_targets_df <- function(FocusedData, Metadata){
  
  
  ##Create a targets file for the FocusedData with information for each sample
  Sample <- colnames(FocusedData)
  INC_Code <- substr(colnames(FocusedData), start = 1, stop = 7)
  Annotation <- str_sub(colnames(FocusedData), 9, -1)
  FUTURE_POLYP_OR_CRC <- Metadata[INC_Code,]$FUTURE_POLYP_OR_CRC
  FUTURE_LESION <- Metadata[INC_Code,]$FUTURE_LESION
  dysplasia <- Metadata[INC_Code,]$dysplasia
  aden_tubular <- Metadata[INC_Code,]$aden_tubular
  aden_villous <- Metadata[INC_Code,]$aden_villous
  aden_tubulovillous <- Metadata[INC_Code,]$aden_tubulovillous
  aden_serrated <- Metadata[INC_Code,]$aden_serrated
  Sex <- Metadata[INC_Code,]$Sex
  AdvAdenoma <- Metadata[INC_Code,]$AdvancedAdenomas
  
  #Merge all data in a targets dataframe
  Targets_df_focused <- data.frame(Sample, INC_Code, Annotation, FUTURE_POLYP_OR_CRC, FUTURE_LESION, dysplasia, 
                           aden_tubular, aden_villous, aden_tubulovillous, aden_serrated, Sex, AdvAdenoma)
  
  #Add a factor for patients, will be used in DE analysis for blocking
  Targets_df_focused <- transform(Targets_df_focused, patient=as.numeric(factor(INC_Code)))

  #Join diff
  Targets_df_focused$Annotation <- as.character(Targets_df_focused$Annotation)
  Targets_df_focused["Annotation"][Targets_df_focused["Annotation"] == "Dyspl.Epi.1"] <- "Dysplasia"
  Targets_df_focused["Annotation"][Targets_df_focused["Annotation"] == "Dyspl.Epi.2"] <- "Dysplasia"
  Targets_df_focused["Annotation"][Targets_df_focused["Annotation"] == "Interface.1"] <- "Interface"
  Targets_df_focused["Annotation"][Targets_df_focused["Annotation"] == "Interface.2"] <- "Interface"
  
  levels(Targets_df_focused$Sex) <- c(levels(Targets_df_focused$Sex), "Male", "Female")
  Targets_df_focused$Sex[Targets_df_focused$Sex == "1"] <- "Male"
  Targets_df_focused$Sex[Targets_df_focused$Sex == "0"] <- "Female"
  
  levels(Targets_df_focused$AdvAdenoma) <- c(levels(Targets_df_focused$AdvAdenoma), "NoAdv", "Adv")
  Targets_df_focused$AdvAdenoma[Targets_df_focused$AdvAdenoma == "1"] <- "Adv"
  Targets_df_focused$AdvAdenoma[Targets_df_focused$AdvAdenoma == "0"] <- "NoAdv"
  
  levels(Targets_df_focused$FUTURE_POLYP_OR_CRC) <- c(levels(Targets_df_focused$FUTURE_POLYP_OR_CRC), "Rec", "NoRec")
  Targets_df_focused$FUTURE_POLYP_OR_CRC[Targets_df_focused$FUTURE_POLYP_OR_CRC == "1"] <- "Rec"
  Targets_df_focused$FUTURE_POLYP_OR_CRC[Targets_df_focused$FUTURE_POLYP_OR_CRC == "0"] <- "NoRec"

  
  
  return(Targets_df_focused)
  
  }

fdata_graphs <- function(Targets_df_focused, Metadata, intersection_f, FocusedData){
  
  ##Create simple graphs to get a view of the samples from the Focused dataset
  TableAnnot_f <- table(Targets_df_focused$Annotation)
  lbls <- paste(names(TableAnnot), "\n", TableAnnot, sep="")
  pie(TableAnnot, labels = lbls,
      main="Pie Chart of number of samples for each polyp section")

  ##Barplot for polyp sections
  barplot(TableAnnot_f, main = "Number of samples for each polyp section in focused dataset",
          ylab = "Number of samples")
  
  ##Barplot for polyp types
  barplot(c(count(Metadata[intersection_f,]$aden_tubulovillous), count(Metadata[intersection_f,]$aden_tubular),
            count(Metadata[intersection_f,]$aden_serrated), count(Metadata[intersection_f,]$aden_villous)), 
          names.arg = c("tubulovillous", "tubular", "serrated", "villous"), 
          main = "Number of polyp types in focused samples")
  
  ##Barplot for recurrence
  TableRecurrence_f <- table(Count_recurrence_f$FUTURE_POLYP_OR_CRC)
  Count_recurrence_f <- Targets_df_focused[match(unique(Targets_df_focused$patient), Targets_df_focused$patient), ]
  barplot(TableRecurrence_f, main = "Recurrence of polyps in focused samples", names.arg = c("No", "Yes"))
  
  ##This is all to prepare for the PCA analysis
  ##Find way of doing this better
  
  Dyspl1 <- grep("*Dyspl.Epi.1",colnames(FocusedData))
  names(Dyspl1) <- rep("red", length(Dyspl1))
  Dyspl2 <- grep("*Dyspl.Epi.2",colnames(FocusedData))
  names(Dyspl2) <- rep("red", length(Dyspl2))
  
  Interface1 <- grep("*Interface.1",colnames(FocusedData))
  names(Interface1) <- rep("blue", length(Interface1))
  Interface2 <- grep("*Interface.2",colnames(FocusedData))
  names(Interface2) <- rep("blue", length(Interface2))
        
  Normal <- grep("*Normal.1",colnames(FocusedData))
  names(Normal) <- rep("green", length(Normal))
  
  NotAnnotated <- grep("*Not.Annotated" ,colnames(FocusedData))
  names(NotAnnotated) <- rep("yellow", length(NotAnnotated))
  
  
  Types <- c(Dyspl1, Dyspl2, Interface1, Interface2, Normal, NotAnnotated)
  Types <- sort(Types)
  Colours <- names(Types)
  
  ## Perform PCA
  
  pca_f <- prcomp(t(FocusedData), scale=T)
  # Plot the PCA results
  #png("PCA_FocusedDataset.png")
  s3d_f<-scatterplot3d(pca_f$x[,1:3], pch=19, color=Colours)
  s3d.coords_f <- s3d$xyz.convert(pca$x[,1:3])
  text(s3d.coords$x, s3d.coords$y,pos = 3,offset = 0.5, col=Colours)
  legend("topleft", pch = 20, col=unique(Colours), legend = c("Dysplasia", "Interface", "Normal", "Non Annotated"), bty='n', cex=1)
  #dev.off()
  
  Dysp_dend <- FocusedData[Dyspl1]
  Norm_dend <- FocusedData[Normal]
  dend_data_f <- cbind(Dysp_dend, Norm_dend)
  hc_f<-hclust(as.dist(1-cor(dend_data_f, method="pearson")), method="average")
  plot(hc_f, cex = 0.7, xlab = "Pearson distance between samples")
  hcd_f <- as.dendrogram(hc_f)
  png("dendrogram_focused.png")
  plot(hcd_f, horiz = TRUE)
  dev.off()
  specific_leaf <- hcd_f[[1]][[1]][[1]]
}

## PERFORM DE ANALYSIS WITH EDGER

simpleDE_EdgeR <- function(FocusedData, Targets_df_focused) {
  ##Simple DE analysis normal vs recurrence
  #Create a DGElist object from a table of counts
  DElist <- DGEList(counts=as.matrix(FocusedData), group=Targets_df_focused$Annotation)
  #Estimate normalization factors
  DElist = calcNormFactors(DElist, method = "upperquartile")
  
  # CAN'T GET THIS GRAPH TO WORK
  #plotMDS(DElist)
  #plotMDS(DElist, labels = Targets_df$Sample,
  #        col = c("darkgreen","blue")[factor(Targets_df$FUTURE_POLYP_OR_CRC)])
  DElist <- estimateCommonDisp(DElist)
  DElist <- estimateTagwiseDisp(DElist)
  
  png("MeanVar.png")
  plotMeanVar(DElist, show.tagwise.vars = TRUE, NBline = TRUE)
  dev.off()
  
  png("plotBCV.png")
  plotBCV(DElist)
  dev.off()
  
  #Simple DE analysis (normal vs. dysplasia)
  de = exactTest(DElist, pair = c("Dysplasia","Normal.1"))
  tt = topTags(de, n = nrow(DElist))
  head(tt$table, 200)
  rn1 = rownames(tt$table)
  deg1 = rn1[tt$table$FDR < .05]
  
  #Plot smear plot of results
  png("smear_normaldyspl.png")
  plotMD(de, de.tags = deg1, main = "Dysplasia vs. normal zones")
  dev.off()
  
  plotMD(de, de.tags = deg1, main = "Dysplasia vs. normal zones")
  
}

#Build a more complex DE analysis
#Multilevel analysis. Section 3.5 EdgeR User's Guide

multilvl_design_EdgeR <- function(Targets_df) {
  
    Patient_f <- factor(Targets_df_focused$patient)
    Annotation_f <- factor(Targets_df_focused$Annotation)
    Recurrence_f <- factor(Targets_df_focused$FUTURE_POLYP_OR_CRC)
    Sex_f <- factor(Targets_df_focused$Sex)
    
    #Define contrasts of interest
    Group <- factor(paste(
                          Targets_df_focused$FUTURE_POLYP_OR_CRC, 
                          Targets_df_focused$Annotation,
                          sep="."))
    
    cbind(Targets_df_focused, Group)
    
    design_f <- model.matrix(~ Patient)
    
    NoRec.Dysplasia <- FUTURE_POLYP_OR_CRC== "0" & Annotation=="Dysplasia"
    NoRec.Interface <- FUTURE_POLYP_OR_CRC== "0" & Annotation=="Interface"
    NoRec.Normal.1 <- FUTURE_POLYP_OR_CRC== "0" & Annotation=="Normal.1"
    NoRec.Not.Annotated <- FUTURE_POLYP_OR_CRC== "0" & Annotation=="Not.Annotated"
    Rec.Dysplasia <- FUTURE_POLYP_OR_CRC== "1" & Annotation=="Dysplasia"
    Rec.Interface <- FUTURE_POLYP_OR_CRC== "1" & Annotation=="Interface"
    Rec.Normal.1 <- FUTURE_POLYP_OR_CRC== "1" & Annotation=="Normal.1"
    Rec.Not.Annotated <- FUTURE_POLYP_OR_CRC== "1" & Annotation=="Not.Annotated"
    
    design_f <- cbind(design_f, NoRec.Dysplasia, NoRec.Normal.1, 
                    Rec.Dysplasia, Rec.Normal.1)
    
    return(design_f)
    
 
  
  }


complexDE_EdgeR <- function(FocusedData, design_f){
    #Account for type of tissue and recurrence
    ComplexDE_f <- DGEList(counts=FocusedData)
    
    #Estimate normalizing factors
    ComplexDE_f <- calcNormFactors(ComplexDE_f, method="upperquartile")
    
    #Estimate dispersion values using CR-adjusted likelihood
    DElist_f = estimateGLMTrendedDisp(ComplexDE_f, design_f)
    # DElist2 = estimateGLMTagwiseDisp(DElist2, design)
    fitcomplex_f <- glmQLFit(DElist_f, design_f)
    
    #Establish the contrasts to study in DE analysis
    contrasts_f <- makeContrasts(NoRec.DysVsNormal = NoRec.Dysplasia - NoRec.Normal.1,
                               Rec.DysvsNormal = Rec.Dysplasia - Rec.Normal.1,
                               Dys.RecvsNoRec = Rec.Dysplasia - NoRec.Dysplasia,
                               Normal.RecvsNoRec = Rec.Normal.1 - NoRec.Normal.1,
                               
                               levels=design_f)

    #CONTRAST 1: NON RECURRENT POLYPS, DYSPLASIA VS NORMAL ZONE
    cat("Now comparing NON RECURRENT PATIENTS, DYSPLASIA VS NORMAL ZONE")
    Contrast1_f = glmQLFTest(fitcomplex_f, contrast = contrasts_f[, "NoRec.DysVsNormal"])
    tt_Contrast1_f = topTags(Contrast1_f, n = nrow(DElist_f))
    head(tt_Contrast1_f$table, 20)
    rn_Contrast1_f = rownames(tt_Contrast1_f$table)
    # deg_Contrast1_f = rn_Contrast1_f[tt_Contrast1_f$table$FDR < 0.05]
    png("FocusedDS_NoRec_DysVSNormal.png")
    plotMD(Contrast1_f, de.tags = deg_Contrast1_f, main = "Non recurrent patients, dysplasia vs normal area of polyp")
    dev.off()
    
    write.csv(tt_Contrast1_f$table, file = "toptags_FOCUSED_NonRec_AdvVSNAdv.csv")
    
    # CONTRAST 2: DYSPLASIA AREAS, RECURRENT VS NON RECURRENT PATIENTS
    cat("now comparing dysplasia areas of rekjcurrent vs non recurrent patients")
    Contrast2_f = glmQLFTest(fitcomplex_f, contrast = contrasts_f[, "Dys.RecvsNoRec"])
    tt_Contrast2_f = topTags(Contrast2_f, n = nrow(DElist_f))
    head(tt_Contrast2_f$table, 20)
    rn_Contrast2_f = rownames(tt_Contrast2_f$table)
    deg_Contrast2_f = rn_Contrast2_f[tt_Contrast2_f$table$FDR < 0.05]
    png("FocusedDS_Dys_RecvsNR.png")
    plotMD(Contrast2_f, de.tags = deg_Contrast2_f, main = "Dysplasia polyp sections, recurrent vs non recurrent patients")
    dev.off()
    
    write.csv(tt_Contrast2_f$table, file = "toptags_FOCUSED_Dysp_RecvsNoRec.csv")
    
    # CONTRAST 3: NORMAL AREAS, RECURRENT VS NON RECURRENT PATIENTS
    cat("now comparing dysplasia areas of recurrent vs non recurrent patients")
    Contrast3_f = glmQLFTest(fitcomplex_f, contrast = contrasts_f[, "Normal.RecvsNoRec"])
    tt_Contrast3_f = topTags(Contrast3_f, n = nrow(DElist_f))
    head(tt_Contrast3_f$table, 20)
    rn_Contrast3_f = rownames(tt_Contrast3_f$table)
    deg_Contrast3_f = rn_Contrast3_f[tt_Contrast3_f$table$FDR < 0.05]
    png("FocusedDS_Norm_RecvsNR.png")
    plotMD(Contrast3_f, de.tags = deg_Contrast3_f, main = "Normal polyp section, recurrent vs non recurrent patients")
    dev.off()
    
    write.csv(tt_Contrast3_f$table, file = "toptags_FOCUSED_Normal_RecvsNoRec.csv")
    
    #CONTRAST 4: RECURRENT POLYPS, DYSPLASIA VS NORMAL ZONE
    cat("Now comparing NON RECURRENT PATIENTS, DYSPLASIA VS NORMAL ZONE")
    Contrast4_f = glmQLFTest(fitcomplex_f, contrast = contrasts_f[, "Rec.DysvsNormal"])
    tt_Contrast4_f = topTags(Contrast4_f, n = nrow(DElist_f))
    head(tt_Contrast4_f$table, 20)
    rn_Contrast4_f = rownames(tt_Contrast4_f$table)
    deg_Contrast4_f = rn_Contrast4_f[tt_Contrast4_f$table$FDR < 0.05]
    png("FocusedDS_Rec_DysVSNormal.png")
    plotMD(Contrast4_f, de.tags = deg_Contrast4_f, main = "Recurrent patients, dysplasia vs normal area of polyp")
    dev.off()
    write.csv(tt_Contrast4_f$table, file = "toptags_FOCUSED_Rec_AdvVSNAdv.csv")
    
    
    Contrast1_f_hm <- tt_Contrast1_f$table[1:30,]
    Indexes1_f <- rownames(Contrast1_f_hm[1:30,])
    plotthis1_f <- cbind(Contrast1_f_hm[Indexes1_f,]$logFC, Contrast4_f$table[Indexes1_f,]$logFC)
    annot_col <- as.matrix(c("Non_rec", "Rec"))
    #Save heatmap as png for report
    png("heatmap_BMP4.png")
    pheatmap(plotthis1_f, labels_row = Indexes1_f, annotation_col = annot_col)
    dev.off()
    
    
    library(VennDiagram)
    
    venn.diagram(
      x = list(deg_Contrast1_f, deg_Contrast2_batch1),
      category.names = c("Focused dataset, dysplasia vs normal" , "Batch 1 dataset, advanced vs non advanced"),
      filename = 'venn_diagramm._focusedvsb1.png',
      output=TRUE )
    
    overlap <- calculate.overlap(x = list(deg_Contrast1_f, deg_Contrast2_batch1))
    
 
    
 
    }

csv_to_rnk <- function(myfile){
    x<-read.table("toptags_edgeR_dysplrecurrence.csv",sep = ",", header=T)
    attach(x)
    x$fcSign=sign(logFC)
    x$logP=-log10(PValue)
    x$metric=x$logP/x$fcSign
    y<-x[,c("X", "metric")]
    write.table(y,file="expression.rnk",quote=FALSE,sep="\t", row.names=FALSE)
}


## Preparations to carry out Functional Enrichment analysis
GO_analysis <- function(fe_input_f) {
  # SET THE DESIRED ORGANISM (Human because INCISE patients are human)
  organism = "org.Hs.eg.db"
  BiocManager::install(organism, character.only = TRUE)
  library(organism, character.only = TRUE)
  
  # Reading data from EdgeR that will be used for functional enrichment analysis
  gene_list_f <- fe_input_f$logFC
  
  # Name the vector
  names(gene_list_f) <- fe_input_f$X
  
  # Omit possible NA values
  gene_list_f <- na.omit(gene_list_f)
  
  # Sort the list in decreasing order (this is necessary for clusterProfiler)
  gene_list_f = sort(gene_list_f, decreasing = TRUE)
  
  
  gse_f <- gseGO(geneList=gene_list_f, 
               ont ="ALL", 
               keyType = "SYMBOL",
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db,
               pAdjustMethod = "BH")
  
  require(DOSE)
  
  png("dotplot_FOCUSED_dysplVSnorm.png")
  dotplot(gse_f, showCategory=10, split=".sign") + facet_grid(.~.sign) 
  dev.off()
  
  png("emapplot_FOCUSED_advVsNA.png")
  emapplot(gse_f) 
  dev.off()
  
  png("ridgeplot_advVsNA.png")
  ridgeplot(gse_f) + labs(x = "enrichment distribution")
  dev.off()
  
  cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
  
  terms <- gse$Description[1:3]
  pmcplot(terms, 2012:2020, proportion=FALSE)
}

KEGG_analysis <- function(fe_input) {
  # Convert gene IDs for gseKEGG function
  # Some genes will be lost due to ID conversion
  ids_F<-bitr(fe_input_f$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
  dedup_ids_f = ids[!duplicated(ids_F[c("SYMBOL")]),]
  
  # Create a new dataframe df_KEGG which has only the genes which were successfully mapped using the bitr function above
  df_KEGG_f = fe_input_f[fe_input_f$X %in% dedup_ids_f$SYMBOL,]
  
  # Create a new column in df_KEGG with the corresponding ENTREZ IDs
  df_KEGG_f$Y = dedup_ids_f$ENTREZID
  
  # Create a vector of the gene universe
  kegg_gene_list_f <- df_KEGG_f$logFC
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list_f) <- df_KEGG_f$Y
  
  # omit any NA values 
  kegg_gene_list_f<-na.omit(kegg_gene_list_f)
  
  # sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list_f = sort(kegg_gene_list_f, decreasing = TRUE)
  
  kegg_organism = "hsa"
  kk2_F <- gseKEGG(geneList = kegg_gene_list_f,
                 organism = kegg_organism,
                 nPerm = 10000,
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 0.1,
                 pAdjustMethod = "BH",
                 keyType = "ncbi-geneid")
  
  
  dotplot(kk2_F, showCategory = 20, title = "Enriched KEGG Pathways" , split=".sign") + facet_grid(.~.sign)
  
  #Enrichment map plot of kegg results
  png("emapplot_KEGG_advVsNA.png")
  emapplot(kk2_F)
  dev.off()
  
  library(pathview)
  
  # Produce the native KEGG plot (PNG)
  dme <- pathview(gene.data=kegg_gene_list, pathway.id="05210", species = kegg_organism)
}


main <- function(){
  cat("Starting analysis, reading TempOSeq data files\n")
  ##Load INCISE data
  #PrunedData <- read.table(file = '../raw_data/INCISE_Batch1_PrunedExpression.tsv', sep = '\t', header = TRUE)
  FocusedData <- read.table(file = '../raw_data/BCL-CO0005b_gene_counts_Focused area_training_samples.tsv', sep = '\t', header = TRUE)
  Metadata <- read.table(file = './raw_data/INCISE_Batch1_Metadata.tsv', sep = '\t', header = TRUE)
  
  Metadata <- load_metadata_df(Metadata)
  FocusedData <- load_data_df(FocusedData)
  Targets_df <- create_targets_df(FocusedData, Metadata)
  
  #cat("Summary graphs of the data done")
  #fdata_graphs(Targets_df, Metadata, intersection, FocusedData)
  
  cat("DE analysis running. Please wait...")
  simpleDE_EdgeR(FocusedData, Targets_df)
  
  fe_input_f = read.csv("toptags_FOCUSED_Rec_AdvVSNAdv.csv", header = TRUE)
  
  cat("Done")
}

main()

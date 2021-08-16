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

##Load INCISE data
PrunedData <- read.table(file = 'INCISE_Batch1_PrunedExpression.tsv', sep = '\t', header = TRUE)
FocusedData <- read.table(file = 'BCL-CO0005b_gene_counts_Focused area_training_samples.tsv', sep = '\t', header = TRUE)
Metadata <- read.table(file = 'INCISE_Batch1_Metadata.tsv', sep = '\t', header = TRUE)

load_data_df <- function(FocusedData, Metadata) {
  
  #Rename columns of the dataframes so that they are named with INCISE code
  Metadata2 <- Metadata[,-1]
  rownames(Metadata2) <- Metadata[,1]
  
  FocusedData2 <- FocusedData[,-1]
  rownames(FocusedData2) <- FocusedData[,1]
  
  ##Remove rows that have all 0 (genes not expressed in any sample)
  FocusedData2 <- FocusedData2[rowSums(FocusedData2[])>10,]
  Quantile10 = rowQuantiles(as.matrix(FocusedData2), probs = 0.1)
  keep.exprs <- filterByExpr(FocusedData2, min.count=10)
  FocusedData2 = FocusedData2[keep.exprs,]

  ##Normalize TempOSeq raw counts to CPM (counts Per Million)
  SumCol = apply(FocusedData2, 2, sum)
  return(FocusedData2)
  return(Metadata2)
  return(Quantile10)
  return(SumCol)
  }

create_targets_df <- function(Metadata2, FocusedData2){
  
  ##Join FocusedData with Metadata
  #Identify which samples are present in the focused dataset
  # intersection <- rownames(Metadata2) %in% substr(colnames(FocusedData2), start = 1, stop = 7)
  intersection2 <- substr(colnames(FocusedData2), start = 1, stop = 7) %in% rownames(Metadata2)
  Metadata2[intersection,]
  
  ##Create a targets file for the FocusedData with information for each sample
  Sample <- colnames(FocusedData2)
  INC_Code <- substr(colnames(FocusedData2), start = 1, stop = 7)
  Annotation <- str_sub(colnames(FocusedData2), 9, -1)
  FUTURE_POLYP_OR_CRC <- Metadata2[INC_Code,]$FUTURE_POLYP_OR_CRC
  FUTURE_LESION <- Metadata2[INC_Code,]$FUTURE_LESION
  dysplasia <- Metadata2[INC_Code,]$dysplasia
  aden_tubular <- Metadata2[INC_Code,]$aden_tubular
  aden_villous <- Metadata2[INC_Code,]$aden_villous
  aden_tubulovillous <- Metadata2[INC_Code,]$aden_tubulovillous
  aden_serrated <- Metadata2[INC_Code,]$aden_serrated
  
  #Merge all data in a targets dataframe
  Targets_df <- data.frame(Sample, INC_Code, Annotation, FUTURE_POLYP_OR_CRC, FUTURE_LESION, dysplasia, 
                           aden_tubular, aden_villous, aden_tubulovillous, aden_serrated)
  
  #Add a factor for patients, will be used in DE analysis for blocking
  Targets_df <- transform(Targets_df, patient=as.numeric(factor(INC_Code)))

  #Join diff
  Targets_df$Annotation <- as.character(Targets_df$Annotation)
  Targets_df["Annotation"][Targets_df["Annotation"] == "Dyspl.Epi.1"] <- "Dysplasia"
  Targets_df["Annotation"][Targets_df["Annotation"] == "Dyspl.Epi.2"] <- "Dysplasia"
  Targets_df["Annotation"][Targets_df["Annotation"] == "Interface.1"] <- "Interface"
  Targets_df["Annotation"][Targets_df["Annotation"] == "Interface.2"] <- "Interface"
  Targets_df <- Targets_df[1:156,]
  
  return(Targets_df)
  return(intersection)
  }

fdata_graphs <- function(Targets_df, Metadata2, intersection, FocusedData2){
  
  ##Create simple graphs to get a view of the samples from the Focused dataset
  TableAnnot <- table(Targets_df$Annotation)
  lbls <- paste(names(TableAnnot), "\n", TableAnnot, sep="")
  pie(TableAnnot, labels = lbls,
      main="Pie Chart of number of samples for each polyp section")

  ##Barplot for polyp sections
  barplot(TableAnnot, main = "Proportion of dysplastic/normal sections in focused samples",
          ylab = "Number of samples")
  
  ##Barplot for polyp types
  barplot(c(count(Metadata2[intersection,]$aden_tubulovillous), count(Metadata2[intersection,]$aden_tubular),
            count(Metadata2[intersection,]$aden_serrated), count(Metadata2[intersection,]$aden_villous)), 
          names.arg = c("tubulovillous", "tubular", "serrated", "villous"), 
          main = "Proportion of polyp types in focused samples")
  
  ##Barplot for recurrence
  TableRecurrence <- table(Metadata2[intersection,]$FUTURE_POLYP_OR_CRC)
  barplot(TableRecurrence, main = "Recurrence of polyps in focused samples", names.arg = c("No", "Yes"))
  
  ##This is all to prepare for the PCA analysis
  ##Find way of doing this better
  
  Dyspl1 <- grep("*Dyspl.Epi.1",colnames(FocusedData2))
  names(Dyspl1) <- rep("red", length(Dyspl1))
  Dyspl2 <- grep("*Dyspl.Epi.2",colnames(FocusedData2))
  names(Dyspl2) <- rep("red", length(Dyspl2))
  
  Interface1 <- grep("*Interface.1",colnames(FocusedData2))
  names(Interface1) <- rep("blue", length(Interface1))
  Interface2 <- grep("*Interface.2",colnames(FocusedData2))
  names(Interface2) <- rep("blue", length(Interface2))
        
  Normal <- grep("*Normal.1",colnames(FocusedData2))
  names(Normal) <- rep("green", length(Normal))
  
  NotAnnotated <- grep("*Not.Annotated" ,colnames(FocusedData2))
  names(NotAnnotated) <- rep("yellow", length(NotAnnotated))
  
  
  Types <- c(Dyspl1, Dyspl2, Interface1, Interface2, Normal, NotAnnotated)
  Types <- sort(Types)
  Colours <- names(Types)
  
  ## Perform PCA
  
  pca <- prcomp(t(FocusedData2), scale=T)
  # Plot the PCA results
  #png("PCA_FocusedDataset.png")
  s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=Colours)
  s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
  text(s3d.coords$x, s3d.coords$y, labels = colnames(FocusedData2),pos = 3,offset = 0.5, col=Colours)
  #dev.off()
}

## PERFORM DE ANALYSIS WITH EDGER

simpleDE_EdgeR <- function(FocusedData2, Targets_df) {
  ##Simple DE analysis normal vs recurrence
  #Create a DGElist object from a table of counts
  DElist <- DGEList(counts=FocusedData2[,1:156], group=Targets_df$Annotation)
  #Estimate normalization factors
  DElist = calcNormFactors(DElist)
  
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
  plotSmear(de, de.tags = deg1)
  dev.off()
}

#Build a more complex DE analysis
#Multilevel analysis. Section 3.5 EdgeR User's Guide

multilvl_design_EdgeR <- function(Targets_df) {
  
    Patient <- factor(Targets_df$patient)
    Annotation <- factor(Targets_df$Annotation)
    Recurrence <- factor(Targets_df$FUTURE_POLYP_OR_CRC)
    
    design <- model.matrix(~Patient)
    
    #Define contrasts of interest
    Dysplasia.recurrent <- Targets_df$Annotation== "Dysplasia" & Targets_df$FUTURE_POLYP_OR_CRC == "1"
    Dysplasia.nonrecurrent <- Targets_df$Annotation== "Dysplasia" & Targets_df$FUTURE_POLYP_OR_CRC == "0"
    Normal.recurrent <- Targets_df$Annotation== "Normal.1" & Targets_df$FUTURE_POLYP_OR_CRC == "1"
    Normal.nonrecurrent <- Targets_df$Annotation== "Normal.1" & Targets_df$FUTURE_POLYP_OR_CRC == "0"
    
    design <- cbind(design, Dysplasia.recurrent, Dysplasia.nonrecurrent, Normal.recurrent, Normal.nonrecurrent)
    output(design)
  }

complexDE_EdgeR <- function(FocusedData2, design){}
    #Account for type of tissue and recurrence
    ComplexDE <- DGEList(counts=FocusedData2[,1:156])
    
    #Estimate normalizing factors
    ComplexDE <- calcNormFactors(ComplexDE)
    
    #Estimate dispersion values using CR-adjusted likelihood
    DElist2 = estimateGLMTrendedDisp(ComplexDE, design)
    DElist2 = estimateGLMTagwiseDisp(DElist2, design)
    fitcomplex <- glmQLFit(DElist2, design)
    
    #Fit a GLM to each feature
    drf <- glmQLFTest(fitcomplex, coef="Dysplasia.recurrent")
    
    ttcomp = topTags(drf, n = nrow(DElist2))
    head(ttcomp$table, 200)
    rn2 = rownames(ttcomp$table) 
    deg2 = rn2[ttcomp$table$FDR < .05]
    
    #Plot smear plot of results
    png("smear_normaldyspl_complex.png")
    plotSmear(drf, de.tags = deg2)
    dev.off()
    
    write.csv(ttcomp$table, file = "toptags_edgeR_recurrence.csv")

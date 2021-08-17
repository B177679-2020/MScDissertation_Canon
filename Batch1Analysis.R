# TRANSCRIPTOMICS ANALYSIS OF BATCH 1 DATASET
# Load in BATCH 1 dataset from INCISE project, link it to metadata, create a targets df
# Preliminary overview of dataset with proportions of polyp type, recurrence...
# PCA analysis.
# DE analysis with EdgeR
# Functional Enrichment Analysis with data from GO / KEGG


##Load libraries 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


##Creates list of packages that are needed to run the script. Checks if they are installed.
##If they are not, install the packages.
list.of.packages <- c("limma", "scatterplot3d", "stringr", "ggplot2", "edgeR", "tidyr", "dplyr", 
                      "org.Hs.eg.db", "clusterProfiler", "enrichplot", "EnrichmentBrowser")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

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
library(EnrichmentBrowser)
library(pathview)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)




load_metadata_df <- function(Metadata){
  #Rename columns of the dataframes so that they are named with INCISE code
  rownames(Metadata) <- Metadata[,1]
  return(Metadata)
}

load_data_df <- function(Data, Metadata) {
  
  names(Data)[1] = "Gene_ID"
  rownames(Data) <- Data[,1]

  
  ##Processing of sample names. TempOSeq data has names like "CD38_4647". There might be several probes for the
  #same gene with different numeric part. This lines of code aim to separate the gene names from the numeric code
  #and group the repeats adding up their counts. Also ensure the index column is called Gene_ID.
  Data <- separate(data = Data, col = Gene_ID, into = c("Gene_ID", "tag"), sep = "_")
  Data <- Data %>% group_by(Gene_ID) %>% summarise(across(where(is.numeric), sum))
  Data <- as.data.frame(Data)
  rownames(Data) <- Data[,1]
  Data <- Data[,-1]
  
  
  #GENE FILTERING BY EXPRESSION
  # Remove rows of genes that have less than 3 raw counts in more than 70% of the samples.
  # This is done to avoid noisy signal due to very low-expressed genes.

  keep.exprs <- filterByExpr(Data, min.count=3)
  Data = Data[keep.exprs,]
  
  ##Join Data with Metadata
  #Identify metadata samples which are present in the Batch dataset
  intersection <- colnames(Data) %in% rownames(Metadata)
  #Keep only data that has metadata associated
  Data <- Data[,c(intersection)]

  return(Data)
  
}


create_targets_df <- function(Data, Metadata){
  
  ##Create a targets file for the Data with information for each sample
  #This targets dataframe will be used for the DE analysis
  INC_Code <- colnames(Data)
  FUTURE_POLYP_OR_CRC <- Metadata[INC_Code,]$FUTURE_POLYP_OR_CRC
  FUTURE_LESION <- Metadata[INC_Code,]$FUTURE_LESION
  dysplasia <- Metadata[INC_Code,]$dysplasia
  aden_tubular <- Metadata[INC_Code,]$aden_tubular
  aden_villous <- Metadata[INC_Code,]$aden_villous
  aden_tubulovillous <- Metadata[INC_Code,]$aden_tubulovillous
  aden_serrated <- Metadata[INC_Code,]$aden_serrated
  Sex <- Metadata[INC_Code,]$Sex
  Dysplasia <- Metadata[INC_Code,]$dysplasia
  AdvAdenoma <- Metadata[INC_Code,]$AdvancedAdenomas
  
  #Merge all data in a targets dataframe
  Targets_df <- data.frame(INC_Code, FUTURE_POLYP_OR_CRC, FUTURE_LESION, dysplasia, 
                           aden_tubular, aden_villous, aden_tubulovillous, aden_serrated,
                           Dysplasia, Sex, AdvAdenoma)
  
  #Add a factor for patients, will be used in DE analysis for blocking
  Targets_df <- transform(Targets_df, patient=as.numeric(factor(INC_Code)))
  
  # Drop the samples that are not fully annotated. If there are NAs then some commands wont run.
  # Target_df <- drop_na(Targets_df)
  
  #Rename certain factors so that instead of 0 and 1 they carry informative info about sex, recurrence and polyp status
  levels(Targets_df$Sex) <- c(levels(Targets_df$Sex), "Male", "Female")
  Targets_df$Sex[Targets_df$Sex == "1"] <- "Male"
  Targets_df$Sex[Targets_df$Sex == "0"] <- "Female"
  
  levels(Targets_df$AdvAdenoma) <- c(levels(Targets_df$AdvAdenoma), "NoAdv", "Adv")
  Targets_df$AdvAdenoma[Targets_df$AdvAdenoma == "1"] <- "Adv"
  Targets_df$AdvAdenoma[Targets_df$AdvAdenoma == "0"] <- "NoAdv"
  
  levels(Targets_df$FUTURE_POLYP_OR_CRC) <- c(levels(Targets_df$FUTURE_POLYP_OR_CRC), "Rec", "NoRec")
  Targets_df$FUTURE_POLYP_OR_CRC[Targets_df$FUTURE_POLYP_OR_CRC == "1"] <- "Rec"
  Targets_df$FUTURE_POLYP_OR_CRC[Targets_df$FUTURE_POLYP_OR_CRC == "0"] <- "NoRec"
  
  return(Targets_df)
 }
  
data_graphs <- function(Targets_df, Metadata, intersection, Data){
    
    ##Barplot for polyp types
    barplot(c(as.numeric(length(which(Targets_df$aden_tubulovillous == 1))), 
              as.numeric(length(which(Targets_df$aden_tubular == 1))),
              as.numeric(length(which(Targets_df$aden_serrated == 1))), 
              as.numeric(length(which(Targets_df$aden_villous == 1)))), 
            names.arg = c("tubulovillous", "tubular", "serrated", "villous"), 
            main = "Number of polyps by polyp type")
    
    ##Barplot for patient recurrence
    TableRecurrence <- table(Targets_df$FUTURE_POLYP_OR_CRC)
    barplot(TableRecurrence, main = "Recurrent patients in batch 1 samples",
            names.arg = c("No", "Yes"))
    
    ##Barplot for patient sex
    TableSex <- table(Targets_df$Sex)
    barplot(TableSex, main = "Patients by sex in batch 1 samples",
            names.arg = c("Female", "Male"))
    
    ##Barplot for polyp status
    TableAdvanced <- table(Targets_df$AdvAdenoma)
    barplot(TableAdvanced, main = "Adenoma status in batch 1 samples",
            names.arg = c("Advanced", "Non advanced"))
    
    ## Perform PCA analysis of the data
    
    Recurrent_samples <- which(Targets_df$FUTURE_POLYP_OR_CRC == "Rec")
    names(Recurrent_samples) <- rep("red", length(Recurrent_samples))
    Nonrec_samples <- which(Targets_df$FUTURE_POLYP_OR_CRC == "NoRec")
    names(Nonrec_samples) <- rep("blue", length(Nonrec_samples))
    PCA_colors <- c(Recurrent_samples, Nonrec_samples)
    PCA_colors <- sort(PCA_colors)
    PCA_colors <- names(PCA_colors)
    
    pca <- prcomp(t(Data), scale=T)
    
    # Plot the PCA results
    png("PCA_Data.png")
    s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=PCA_colors,cex.lab = 1, cex.axis = 1, main = "PCA analysis of data by patient recurrence")
    s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
    legend("topleft", pch = 20, col=unique(PCA_colors), legend = c("FUTURE POLYP OR CRC", "NO FUTURE POLYP OR CRC"), bty='n', cex=1)
    
    dev.off()
    
    return(Targets_df)
}

#Multilevel analysis design with EdgeR. Section 3.5 EdgeR User's Guide.
#Create a design matrix that is fed to EdgeR.

multilvl_design_EdgeR <- function(Targets_df) {
  
  ##Define the factors of interest for the DE analysis
  Recurrence <- factor(Targets_df$FUTURE_POLYP_OR_CRC)
  Sex <- factor(Targets_df$Sex)
  Advanced <- factor(Targets_df$AdvAdenoma)
  
  
  #Define contrasts of interest
  Group <- factor(paste(Targets_df$Sex,
                        Targets_df$FUTURE_POLYP_OR_CRC, 
                        Targets_df$AdvAdenoma,
                        
                        sep="."))
  
  cbind(Targets_df, Group)
  
  ##Create a design matrix for DE analysis in EdgeR
  design <- model.matrix(~0+Group)
  colnames(design) <- levels(Group)
  return(design)
  
}

complexDE_EdgeR <- function(Data, design){
  #Account for type of tissue and recurrence
  cat("creating a dgelist object...")
  DEList <- DGEList(counts=Data)
  
  #Estimate normalizing factors
  cat("estimating normalizing factors...")
  ComplexDE <- calcNormFactors(DEList, method = "upperquartile")
  
  #Estimate dispersion values using Cox Reid-adjusted likelihood
  cat("Estimating data dispersion...")
  ComplexDE = estimateGLMTrendedDisp(ComplexDE, design)
  # DElist2 = estimateGLMTagwiseDisp(DElist2, design)
  fitcomplex <- glmQLFit(ComplexDE, design)
  
  #Establish the contrasts to study in DE analysis
  contrasts <- makeContrasts(Male.RecvsNonRec.Adv = Male.Rec.Adv - Male.NoRec.Adv,
                             Female.RecvsNonRec.Adv = Female.Rec.Adv - Female.NoRec.Adv,
                             Adv.RecvsNonadv.Rec = Male.Rec.Adv - Male.Rec.NoAdv,
                             Fem.rec.AdvvsNA = Female.Rec.Adv - Female.Rec.NoAdv,
                             levels=design)
  
  #CONTRAST 1: ADVANCED MALE, REC VS NO REC POLYP
  cat("Now comparing advanced polyps from male patients, recurrent vs. non recurrent patients")
  Contrast1 = glmQLFTest(fitcomplex, contrast = contrasts[, "Male.RecvsNonRec.Adv"])
  tt_Contrast1 = topTags(Contrast1, n = nrow(ComplexDE))
  head(tt_Contrast1$table, 20)
  rn_Contrast1 = rownames(tt_Contrast1$table)
  deg_Contrast1 = rn_Contrast1[tt_Contrast1$table$FDR < 0.05]
  png("Male.RecvsNoRec.Adv.png")
  plotMD(Contrast1, de.tags = deg_Contrast1, main = "Male advanced polyps, recurrent vs. non recurrent patients")
  dev.off()
  
  #CONTRAST 2: RECURRENT MALE, ADVANCED VS NO ADVANCED
  cat("now comparing recurrent male patients, advanced vs. non advanced polyps")
  Contrast2 = glmQLFTest(fitcomplex, contrast = contrasts[, "Adv.RecvsNonadv.Rec"])
  tt_Contrast2 = topTags(Contrast2, n = nrow(ComplexDE))
  head(tt_Contrast2$table, 20)
  rn_Contrast2 = rownames(tt_Contrast2$table)
  deg_Contrast2 = rn_Contrast2[tt_Contrast2$table$FDR < 0.05]
  
  write.csv(tt_Contrast2$table, file = "toptags_recurrence_advancedvsNA_MALE.csv")
  
  png("Male.Rec.AdvvsNoAdv.png")
  plotMD(Contrast2, de.tags = deg_Contrast2, main = "Male recurrent patients, advanced vs non advanced polyps")
  dev.off()
  
  #CONTRAST 3: FEMALE ADVANCED, RECURRENT VS. NON RECURRENT
  cat("now comparing advanced polyps from females, recurrent vs. non recurrent patients")
  Contrast3 = glmQLFTest(fitcomplex, contrast = contrasts[, "Female.RecvsNonRec.Adv"])
  tt_Contrast3 = topTags(Contrast3, n = nrow(ComplexDE))
  head(tt_Contrast3$table, 20)
  rn_Contrast3 = rownames(tt_Contrast3$table)
  deg_Contrast3 = rn_Contrast3[tt_Contrast3$table$FDR < 0.05]
  
  write.csv(tt_Contrast3$table, file = "toptags_advanced_RecVsNR_FEMALE.csv")
  
  png("Fem.Adv.RecvsNR.png")
  plotMD(Contrast3, de.tags = deg_Contrast3, main = "Female advanced polyps, recurrent vs. non recurrent patients")
  dev.off()
  
  #CONTRAST 4: RECURRENT FEMALE, ADVANCED VS NO ADVANCED
  cat("now comparing recurrent female patients, advanced vs. non advanced polyps")
  Contrast4 = glmQLFTest(fitcomplex, contrast = contrasts[, "Fem.rec.AdvvsNA"])
  tt_Contrast4 = topTags(Contrast4, n = nrow(ComplexDE))
  head(tt_Contrast4$table, 20)
  rn_Contrast4 = rownames(tt_Contrast4$table)
  deg_Contrast4 = rn_Contrast4[tt_Contrast4$table$FDR < 0.05]
  
  write.csv(tt_Contrast4$table, file = "toptags_recurrence_advancedvsNA_FEMALE.csv")
  
  png("Fem.Rec.AdvvsNE.png")
  plotMD(Contrast4, de.tags = deg_Contrast4, main = "Female recurrent patients, advanced vs. non advanced polyps")
  dev.off()
  
}

csv_to_rnk <- function(myfile){
  x<-read.table(myfile ,sep = ",", header=T)
  attach(x)
  x$fcSign=sign(logFC)
  x$logP=-log10(PValue)
  x$metric=x$logP/x$fcSign
  y<-x[,c("X", "metric")]
  write.table(y,file="expression.rnk",quote=FALSE,sep="\t", row.names=FALSE)
}
    
## Preparations to carry out Functional Enrichment analysis
GO_analysis <- function(fe_input) {
# SET THE DESIRED ORGANISM (Human because INCISE patients are human)

    
    # Reading data from EdgeR that will be used for functional enrichment analysis
    gene_list <- fe_input$logFC
    
    # Name the vector
    names(gene_list) <- fe_input$X
    
    # Omit possible NA values
    gene_list <- na.omit(gene_list)
    
    # Sort the list in decreasing order (this is necessary for clusterProfiler)
    gene_list = sort(gene_list, decreasing = TRUE)
    
    ##Perform GSEA analysis with gene sets taken from Gene Ontology (GO)
    gse <- gseGO(geneList=gene_list, 
                 ont ="ALL", 
                 keyType = "SYMBOL",
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.1, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db,
                 pAdjustMethod = "BH")
    
    ##Load DOSE database
    require(DOSE)
    
    ##Plot top 20 enriched GO terms split by sign (activated or deactivated)
    png("../graphs/dotplot_advVsNA.png")
    dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign) 
    dev.off()
    
    ##Plot network with relationship between enriched GO plots
    png("../graphs/emapplot_advVsNA.png")
    print(emapplot(gse, showCategory = 20, min_edge = 0.6))
    dev.off()
    
    png("../graphs/ridgeplot_advVsNA.png")
    print(ridgeplot(gse) + labs(x = "enrichment distribution"))
    dev.off()
    
    cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
    
    terms <- gse$Description[1:3]
    pmcplot(terms, 2012:2020, proportion=FALSE)
}


##FUNCTION TO PERFORM GSEA ANALYSIS WITH TERMS FROM THE KEGG DATABASE
KEGG_analysis <- function(fe_input) {
    # Convert gene IDs for gseKEGG function
    # Some genes will be lost due to ID conversion
    ids<-bitr(fe_input$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
    # remove duplicate IDS (here I use "SYMBOL", but it should be whatever was selected as keyType)
    dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
    
    # Create a new dataframe df_KEGG which has only the genes which were successfully mapped using the bitr function above
    df_KEGG = fe_input[fe_input$X %in% dedup_ids$SYMBOL,]
    
    # Create a new column in df_KEGG with the corresponding ENTREZ IDs
    df_KEGG$Y = dedup_ids$ENTREZID
    
    # Create a vector of the gene universe
    kegg_gene_list <- df_KEGG$logFC
    
    # Name vector with ENTREZ ids
    names(kegg_gene_list) <- df_KEGG$Y
    
    # omit any NA values 
    kegg_gene_list<-na.omit(kegg_gene_list)
    
    # sort the list in decreasing order (required for clusterProfiler)
    kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
    
    kegg_organism = "hsa"
    kk2 <- gseKEGG(geneList = kegg_gene_list,
                   organism = kegg_organism,
                   nPerm = 10000,
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = 0.20,
                   pAdjustMethod = "BH",
                   keyType = "ncbi-geneid")
    
    png("../graphs/dotplot_KEGG_AdvVsNA_try2.png")
    dotplot_kegg <- dotplot(kk2, showCategory = 20, title = "Enriched KEGG Pathways" , split=".sign") + 
      facet_grid(.~.sign)
    print(dotplot_kegg)
    dev.off()
    
    #Enrichment map plot of kegg results
    png("../graphs/emapplot_KEGG_advVsNA.png")
    print(emapplot(kk2, cex = 1.2))
    dev.off()
    
    # Produce the native KEGG plot (PNG)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="05210", species = kegg_organism)
}

#Run the whole analysis
main <- function(){
  cat("Starting analysis, reading TempOSeq data files")
  ##Load INCISE data
  Data <- read.table(file = '../raw_data/INCISE_Batch1_PrunedExpression.tsv', sep = '\t', header = TRUE)
  Metadata <- read.table(file = '../raw_data/INCISE_Batch1_Metadata.tsv', sep = '\t', header = TRUE)
  
  cat("Loading metadata")
  Metadata <- load_metadata_df(Metadata)
  cat("Loading data")
  Data <- load_data_df(Data, Metadata)
  # Creating targets dataframe
  cat("Creating targets dataframe")
  Targets_df <- create_targets_df(Data, Metadata)
  
  # cat("Summary graphs of the data done")
  # data_graphs(Targets_df, Metadata, intersection, Data)
  
  cat("DE analysis running. Please wait...")
  cat("Initializing design matrix...")
  design <- multilvl_design_EdgeR(Targets_df)
  
  complexDE_EdgeR(Data, design)
  
  #cat("Creating RNK files...")
  #csv_to_rnk("toptags_recurrenceVsNR.csv")
  
  fe_input = read.csv("toptags_recurrence_advancedvsNA_MALE.csv", header = TRUE)
   
  cat("Performing GO analysis...")
  GO <- GO_analysis(fe_input)
  
  cat("Performing KEGG analysis...")
  KEGG <- KEGG_analysis(fe_input)
  
  cat("DE and functional enrichment analysis done. All results and graphs are now in your working directory")
  
  cat("Done")
}

main()

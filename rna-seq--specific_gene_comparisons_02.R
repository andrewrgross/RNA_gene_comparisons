### RNA-seq--gene_comparisons 
### Andrew R Gross, 2016-08-25
### This script is intended to read in rna seq data, subset genes of interest, and compare them within the dataset
### And against references from Gtex

########################################################################
### Header
########################################################################

library(gplots)
library(RColorBrewer)
library(biomaRt)
library(DESeq2)
library(Hmisc)

ensembl = useMart(host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL')
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
#listMarts(host="www.ensembl.org")
#listDatasets(ensembl)
#filters = listFilters(ensembl)
#attributes = listAttributes(ensembl)

########################################################################
### Functions
########################################################################

addMedSD <- function(dataframe) {
  median <- apply(dataframe,1,median)
  sd <- apply(dataframe,1,sd)
  return(data.frame(dataframe,median,sd))  
}

sortByMed <- function(dataframe) {
  order <- order(dataframe$median,decreasing=TRUE)
  return(dataframe[order,])
}

convertIDs <- function(dataframe) {
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs
  return(dataframe)
}
addGene <- function(dataframe) {
  genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) {
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}
########################################################################
### Import Data
########################################################################

# Normalized refseq data in units of TPM
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

# Metadata
metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)
metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/HT_plus_reference_metadata.csv", row.names=1)

# Import genes of interest
genes.of.interest <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/Uthras_genes_of_interest.txt")

# Import housekeeping genes
housekeeping.genes <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/Housekeeping_genes.txt")

# Import references
references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

########################################################################
### Format
########################################################################

housekeeping.genes <- housekeeping.genes[c(1,2,3,4,7,8,9,10,11,13,14),]  # Subselect housekeeping genes to use

references <- convertIDs(references)  # Remove decimal values from Ensembl IDs

ht.reference <- references[c(1,17)]  # Subselect from the references for just hypothalamus
ht.reference <- ht.reference[order(ht.reference[2],decreasing = TRUE),]  # Sort hypothalamus from low to high

genes.of.interest <- ht.reference[1:20,]  # Define the genes of interest as top hypothalamic genes
genes.of.interest$Ensembl.ID <- row.names(genes.of.interest)

########################################################################
### Convert references to TPM
########################################################################

ht.reference$Brain...Hypothalamus <- ht.reference$Brain...Hypothalamus/sum(ht.reference$Brain...Hypothalamus)*1000000

########################################################################
### Subsample rows
########################################################################

gene.positions <- match(genes.of.interest$Ensembl.ID,row.names(TPMdata))  # Declare the row numbers in the rnaseq data which correspond to genes in or list of genes of interest

genes.df <- TPMdata[gene.positions,]  # Make a dataframe containing just the rows of RNAseq data corresponding to genes of interest

########################################################################
### Rename row names with genes
########################################################################

#genes.df <- addGene(genes.df)  # Add gene names to the rows
#genes.df$Gene <- genes.of.interest$Description
genes.df$Gene <- genes.of.interest$Gene
row.names(genes.df) <- genes.df$Gene  # Replace the row names with the gene names
genes.df <- genes.df[1:ncol(genes.df)-1]

########################################################################
### Select columns to plot
########################################################################

#genes.df$Reference <- ht.reference$Brain...Hypothalamus[1:nrow(genes.df)]
boxplot.df <- data.frame(t(genes.df))
metaData$Disease <- as.character(metaData$Disease)
metaData$Disease[grep('Adult',metaData$Source)] <- 'Adult HT'
metaData$Disease[grep('MN',metaData$Type)] <- 'iMN'
boxplot.df$disease<-metaData$Disease[match(row.names(boxplot.df),row.names(metaData))]
boxplot.df <- boxplot.df[c(2,3,4,19,20,1,11,14,15,16,17,18,7,8,9,10,12,5,6),]
boxplot.df$Sample <- factor(row.names(boxplot.df), levels=unique(row.names(boxplot.df)))

########################################################################
### Normalize
########################################################################

#for (column in 1:13) {
#  barplot.df[column] <- round(barplot.df[column]*1000/barplot.df$GAPDH,3)
#}

########################################################################
### Plot boxes
########################################################################
### Manual generation of single plot


gene = 6

g <- ggplot(data = boxplot.df, aes(x = disease, y = boxplot.df[gene], fill=disease)) + 
  geom_boxplot(varwidth = TRUE) +
  ggtitle(paste("Expression levels of ",names(boxplot.df[gene]))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_jitter(width = 0)

#plot.list[[length(plot.list)+1]] <- g
#png(paste0("plot_",gene,".png"))
g
#dev.off()
#dev.print(png,paste0("plot_",gene,".png"))
f <- ggplot(data = barplot.df, aes(x = Sample, y = barplot.df[gene], fill=disease)) + 
  geom_bar(stat = "identity") +
  ggtitle(paste("Expression levels of ",names(barplot.df[gene]))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

f


### Automatic generation of all plots

plot.list <- list()

for (gene in 1:13){
  print(gene)
  f <- ggplot(data = barplot.df, aes(x = Sample, y = barplot.df[gene], fill=disease)) + 
    geom_bar(stat = "identity") +
    ggtitle(paste("Expression levels of ",names(barplot.df[gene]))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  Sys.sleep(1)
  plot.list[[length(plot.list)+1]] <- f
  #print(barplot.df[gene][1,])
}

plot.list[[1]]
plot.list[[2]]
plot.list[[3]]
plot.list[[4]]
plot.list[[5]]
plot.list[[6]]
plot.list[[7]]
plot.list[[8]]
plot.list[[9]]
plot.list[[10]]
plot.list[[11]]
plot.list[[12]]
plot.list[[13]]




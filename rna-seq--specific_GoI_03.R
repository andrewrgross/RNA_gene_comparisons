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
library(grid)
library(gridExtra)

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
#metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/HT_plus_reference_metadata.csv", row.names=1)

# Import genes of interest
uthras.genes <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/Uthras_genes_of_interest.txt")

# Import housekeeping genes
housekeeping.genes <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/Housekeeping_genes.txt")

# Import hypothalamic genes at pSI 0.
ht.genes_0.01.df <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/ht.genes_0.01.csv")
ht.genes_0.005.df <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/ht.genes_0.005.csv")
ht.genes_0.0005.df <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/ht.genes_0.0005.csv")
ht.genes_0.0001.df <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/ht.genes_0.0001.csv")

# Import references
references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

########################################################################
### Format
########################################################################

#references <- convertIDs(references)  # Remove decimal values from Ensembl IDs
#ht.reference <- references[c(1,17)]  # Subselect from the references for just hypothalamus
#ht.reference <- ht.reference[order(ht.reference[2],decreasing = TRUE),]  # Sort hypothalamus from low to high
#genes.of.interest <- ht.reference[1:20,]  # Define the genes of interest as top hypothalamic genes
#genes.of.interest$Ensembl.ID <- row.names(genes.of.interest)

### Convert references to TPM
ht.reference$Brain...Hypothalamus <- ht.reference$Brain...Hypothalamus/sum(ht.reference$Brain...Hypothalamus)*1000000

metaData$Disease[grep('Adult',metaData$Source)] <- 'unknown'

### Reorder columns
TPMdata <- TPMdata[c(7,8,9,10,12,11,1,16,18,14,15,17,4,20,19,3,2,5,6)]


########################################################################
### Select comparison set
########################################################################

genes.of.interest <- uthras.genes

genes.of.interest <- housekeeping.genes

genes.of.interest <- ht.genes_0.0001.df

########################################################################
### Subsample rows
########################################################################

gene.positions <- match(genes.of.interest$Ensembl.ID,row.names(TPMdata))  # Declare the row numbers in the rnaseq data which correspond to genes in or list of genes of interest
genes.df <- TPMdata[gene.positions,]  # Make a dataframe containing just the rows of RNAseq data corresponding to genes of interest
### Rename the rows with genes
row.names(genes.df) <- genes.of.interest$Gene  # Replace the row names with the gene names

########################################################################
### Prepare BARplot data
########################################################################

#genes.df$Reference <- ht.reference$Brain...Hypothalamus[1:nrow(genes.df)]
barplot.df <- data.frame(t(genes.df))
barplot.df$disease<-factor(metaData$Disease[match(row.names(barplot.df),row.names(metaData))], levels = c('unknown','CTR','OBS','iMN'))
barplot.df$Sample <- factor(row.names(barplot.df), levels = names(genes.df))

########################################################################
### Prepare BOXplot data
########################################################################

boxplot.df<-barplot.df

########################################################################
### Normalize
########################################################################

#for (column in 1:13) {
#  barplot.df[column] <- round(barplot.df[column]*1000/barplot.df$GAPDH,3)
#}

########################################################################
### Plot
########################################################################
###### Bars

gene = 1

g <- ggplot(data = barplot.df, aes(x = Sample, y = barplot.df[gene], fill=disease)) + 
  geom_bar(stat = "identity") +
  ggtitle(paste("Expression levels of ",names(barplot.df[gene]))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

###### Boxes

h <- ggplot(data = boxplot.df, aes(x = disease, y = boxplot.df[gene], fill=disease)) + 
  geom_boxplot(varwidth = FALSE) +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle(paste("Expression levels of ",names(boxplot.df[gene]))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_jitter(width = 0)

grid.arrange(g,h, ncol = 2)

########################################################################
### Output
########################################################################





########################################################################
### Automatic generation of all plots
########################################################################

generate.bargraphs <- function(dataframe) {
  plot.list <- list()
  for(gene.column in 1:(ncol(dataframe)-2)){
    current.plot <- ggplot(data = dataframe, 
                           aes_string(x = 'Sample', y = names(dataframe)[gene.column], fill = "Disease")) +
      geom_bar(stat = 'identity') +
      labs(title = names(dataframe[gene.column])) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x = element_blank(),
            axis.title.y = element_blank(), legend.position = 'none',
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            panel.border = element_rect(color = "black", fill = NA)) +
      scale_fill_brewer(palette = 'Set1')
    plot.list[[length(plot.list)+1]] <- current.plot
  }
  print(paste("List contains",length(plot.list),"plots"))
  return(plot.list)
}

plot.list <- generate.bargraphs(barplot.df)

plot.list2[[2]]

########################################################################
### Format plots into tiles
########################################################################

(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
  plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],plot.list[[11]],plot.list[[12]],
  ncol = 4,top="Housekeeping genes"))


########################################################################
### Output tiled figure
########################################################################













####
plot.list <- list()

for (gene in 1:13){
  print(gene)
  #barplot.temp <- barplot.df[c(gene,ncol(boxplot.df)-1,ncol(boxplot.df))]
  #print(barplot.temp)
  barplot.temp <- barplot.df
  f = ggplot(data = barplot.temp, aes_string(x = "Sample", y = names(barplot.df)[gene], fill="disease")) + 
    geom_bar(stat = "identity") +
    ggtitle(names(barplot.df[gene])) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none")
  #Sys.sleep(1)
  print(f$labels$title)
  print(f$data[,1])
  plot.list[[length(plot.list)+1]] = f
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






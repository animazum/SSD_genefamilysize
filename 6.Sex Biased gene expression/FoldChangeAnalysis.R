
###### fold change Analysis ######
# Script created by Benjamin

rm(list=ls())
library(logr)

setwd("~/Dropbox/SSD/foldChange/")

## positive or negative
testtype<- "positive"
# brain, f_brain or many others
tissueOFinterest <- "f_brain"
## adult, fetus
mono<-"fetus"
### 
# 
PGLSbackground<-read.csv("~/Dropbox/SSD/foldChange/PGLS.background.genefams2geneID.BG.csv", row.names = 1)

expressioGenes<-read.csv("~/Dropbox/SSD/foldChange/BrainSpan_prenatalAges_MeanBrainStructure_WholeBrainExpression_Male_Female.csv")

tmp <- file.path(getwd(), paste("FoldchangeAnalysis.",mono,".",testtype,".associated.genes.",tissueOFinterest,".tissue.log", sep = ""))
lf <- log_open(tmp)


if (testtype == "negative" &&  mono == "adult") {
  log_print("Captain's log starts...")
  log_print("Using negatively selected genes", hide_notes = T)
  expressioGenes<-read.csv("~/Dropbox/SSD/foldChange/BrainSpan_plus18Ages_MeanBrainStructure_WholeBrainExpression_Male_Female.csv")
  
  # apply log transformation to gene expression data
  expressioGenes[2:9]<-log2(expressioGenes[2:9] + 1)
  
  # Create create average gene expression
  expressioGenes$average.male<- rowMeans(expressioGenes[2:5])
  expressioGenes$average.female<- rowMeans(expressioGenes[6:9])
  expressioGenes$fold.Change <- log2(expressioGenes$average.female/expressioGenes$average.male)
  
  # delete Inf and -Inf 
  expressioGenes<-expressioGenes[!is.infinite(expressioGenes$fold.Change),]
  expressioGenes<-expressioGenes[!is.na(expressioGenes$fold.Change),]
  
  ### remove the means with 0s, becase for the los2(fem/male) we would have some issues with infinite numbers
  background.expressioGenes<- expressioGenes[expressioGenes$Gene_stable_ID %in% PGLSbackground$Ensembl_gene_identifier,]
  
  NEGgenes<- read.csv("~/Dropbox/SSD/sex_biased_genexpression/gene_expression/PGLSoutput.genefams2geneID.NEGATIVE.csv")
  NEGgenes<- NEGgenes[NEGgenes$Ensembl_gene_identifier %in% background.expressioGenes$Gene_stable_ID,]
  
  # fantom <- read.csv("ADULT_UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv") #Read in the fantom data
  log_print("Using Adult selected genes", hide_notes = T)
} else if (testtype == "negative" &&  mono == "fetus") {
  log_print("Captain's log starts...")
  log_print("Using negatively selected genes", hide_notes = T)
  expressioGenes<-read.csv("~/Dropbox/SSD/foldChange/BrainSpan_prenatalAges_MeanBrainStructure_WholeBrainExpression_Male_Female.csv")
  
  # apply log transformation to gene expression data
  expressioGenes[2:16]<-log2(expressioGenes[2:16] + 1)
  
  # Create create average gene expression
  expressioGenes$average.male<- rowMeans(expressioGenes[2:8])
  expressioGenes$average.female<- rowMeans(expressioGenes[9:16])
  expressioGenes$fold.Change <- log2(expressioGenes$average.female/expressioGenes$average.male)
  
  # delete Inf and -Inf 
  expressioGenes<-expressioGenes[!is.infinite(expressioGenes$fold.Change),]
  expressioGenes<-expressioGenes[!is.na(expressioGenes$fold.Change),]
  ### remove the means with 0s, becase for the los2(fem/male) we would have some issues with infinite numbers
  background.expressioGenes<- expressioGenes[expressioGenes$Gene_stable_ID %in% PGLSbackground$Ensembl_gene_identifier,]
  
  NEGgenes<- read.csv("~/Dropbox/SSD/sex_biased_genexpression/gene_expression/PGLSoutput.genefams2geneID.NEGATIVE.csv")
  NEGgenes<- NEGgenes[NEGgenes$Ensembl_gene_identifier %in% background.expressioGenes$Gene_stable_ID,]
  
  log_print("Using fetous selected genes", hide_notes = T)
} else if (testtype == "positive" &&  mono == "adult") {
  log_print("Captain's log starts...")
  log_print("Using positive selected genes", hide_notes = T)
  expressioGenes<-read.csv("~/Dropbox/SSD/foldChange/BrainSpan_plus18Ages_MeanBrainStructure_WholeBrainExpression_Male_Female.csv")
  # apply log transformation to gene expression data
  expressioGenes[2:9]<-log2(expressioGenes[2:9] + 1)
  
  # Create create average gene expression
  expressioGenes$average.male<- rowMeans(expressioGenes[2:5])
  expressioGenes$average.female<- rowMeans(expressioGenes[6:9])
  expressioGenes$fold.Change <- log2(expressioGenes$average.female/expressioGenes$average.male)
  
  # delete Inf and -Inf 
  expressioGenes<-expressioGenes[!is.infinite(expressioGenes$fold.Change),]
  expressioGenes<-expressioGenes[!is.na(expressioGenes$fold.Change),]
  
  ### remove the means with 0s, becase for the los2(fem/male) we would have some issues with infinite numbers
  background.expressioGenes<- expressioGenes[expressioGenes$Gene_stable_ID %in% PGLSbackground$Ensembl_gene_identifier,]
  
  NEGgenes<- read.csv("~/Dropbox/SSD/sex_biased_genexpression/gene_expression/PGLSoutput.genefams2geneID.POSITIVE.csv")
  NEGgenes<- NEGgenes[NEGgenes$Ensembl_gene_identifier %in% background.expressioGenes$Gene_stable_ID,]
  
  # fantom <- read.csv("ADULT_UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv") #Read in the fantom data
  log_print("Using Adult selected genes", hide_notes = T)
} else if (testtype == "positive" &&  mono == "fetus") {
  log_print("Captain's log starts...")
  log_print("Using positive selected genes", hide_notes = T)
  expressioGenes<-read.csv("~/Dropbox/SSD/foldChange/BrainSpan_prenatalAges_MeanBrainStructure_WholeBrainExpression_Male_Female.csv")
  # apply log transformation to gene expression data
  expressioGenes[2:16]<-log2(expressioGenes[2:16] + 1)
  
  # Create create average gene expression
  expressioGenes$average.male<- rowMeans(expressioGenes[2:8])
  expressioGenes$average.female<- rowMeans(expressioGenes[9:16])
  expressioGenes$fold.Change <- log2(expressioGenes$average.female/expressioGenes$average.male)
  
  # delete Inf and -Inf 
  expressioGenes<-expressioGenes[!is.infinite(expressioGenes$fold.Change),]
  expressioGenes<-expressioGenes[!is.na(expressioGenes$fold.Change),]
  ### remove the means with 0s, becase for the log2(fem/male) we would have some issues with infinite numbers
  background.expressioGenes<- expressioGenes[expressioGenes$Gene_stable_ID %in% PGLSbackground$Ensembl_gene_identifier,]
  
  NEGgenes<- read.csv("~/Dropbox/SSD/sex_biased_genexpression/gene_expression/PGLSoutput.genefams2geneID.POSITIVE.csv")
  NEGgenes<- NEGgenes[NEGgenes$Ensembl_gene_identifier %in% background.expressioGenes$Gene_stable_ID,]
  
  # fantom <- read.csv("ADULT_UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv") #Read in the fantom data
  log_print("Using fetous selected genes", hide_notes = T)
}


### creating the index column to extract the genes from the back ground and use them in the analysis
rownames(background.expressioGenes)<- background.expressioGenes$Gene_stable_ID
rownames(NEGgenes)<-NEGgenes$Ensembl_gene_identifier
NEGgenes<-background.expressioGenes[(rownames(background.expressioGenes)) %in% rownames(NEGgenes),]
number_of_genes <- dim(NEGgenes)[1]


############## START HERE !!!!!!!!
# Convert all columns to numeric type
for (i in 2:ncol(NEGgenes)) {
  NEGgenes[, i] <- as.numeric(as.character(NEGgenes[, i]))
}

# Verify the changes
str(NEGgenes)  

##### Change specs
NEGgenes$pval<-0
for (i in 1:dim(NEGgenes)[1]) {
  
  # wilcoxtest<-wilcox.test(as.numeric(NEGgenes[i, c(2:16)], mu = 0)) # Prenatal
  wilcoxtest<-wilcox.test(as.numeric(NEGgenes[i, c(2:9)]), mu =0) # adult
  value_p <- wilcoxtest$p.value
  # NEGgenes[i, 20]<- as.numeric(value_p) # prenatal
  NEGgenes[i, "pval"]<- as.numeric(value_p) # adult
}

alpha <- 0.05  # Set your desired significance level

# False Discovery Rate (FDR) correction
library(stats)

# Calculate adjusted p-values using the FDR method
NEGgenes$FDR_p <- p.adjust(NEGgenes$pval, method = "fdr")

# Identify significant results after FDR correction
significant_results_fdr <- subset(NEGgenes,  NEGgenes$FDR_p < alpha)

max(significant_results_fdr$fold.Change)
min(significant_results_fdr$fold.Change)

sum(significant_results_fdr$fold.Change < -2)


# Creating a plot that matches the fold change analysis and the genes with their GO categories.

GOcategories<-read.csv("/Users/lisalisa/Library/CloudStorage/Dropbox/SSD/GO.terms.table.Process.mammals.csv")

PGLSbackground<-read.csv("~/Dropbox/SSD/foldChange/PGLS.background.genefams2geneID.BG.csv", row.names = 1)

PGLSbackground<-PGLSbackground[PGLSbackground$Ensembl_gene_identifier %in% rownames(NEGgenes),]
GOcategories<- GOcategories[GOcategories$GeneID %in% PGLSbackground$genefams2geneID.BG.GeneID,]


gene_data<- merge(x = PGLSbackground, y = GOcategories, by.x = "genefams2geneID.BG.GeneID", by.y = "GeneID", all = T)
gene_data<- merge(x = gene_data, y = NEGgenes, by.x = "Ensembl_gene_identifier", by.y = "Gene_stable_ID", all = T)

# Get the current date and time
current_date <- Sys.time()
# Format the date and time to your desired format
formatted_date <- format(current_date, "%Y-%m-%d_%H-%M-%S")
write.csv(gene_data, file = paste("GenelistVsGOcat", testtype, tissueOFinterest, mono, formatted_date,"csv", sep = "."))

gene_data<-gene_data[gene_data$FDR_p < 0.05,]

lista<- c("neural crest cell fate specification", "forebrain development", 
          "central nervous system development", "hypothalamus development", 
          "pituitary gland development", "central nervous system myelination", 
          "regulation of myelination", "neuron fate commitment",
          "spinal cord motor neuron cell fate specification", 
          "negative regulation of neuron differentiation", "neuron differentiation",
          "spinal cord motor neuron migration", "neuron differentiation"
          , "visceral motor neuron differentiation"
          , "negative regulation of neuron differentiation"
          , "neuron fate specification")

gene_data<-gene_data[gene_data$GO_term %in%  lista,]

library(ggplot2)

# Create the plot with x-axis as Fold Change and y-axis as GO Category
gene_plot <- ggplot(gene_data, aes(x = fold.Change, y = GO_term, color = FDR_p)) +
  geom_point(size = 3) +  # Scatter plot of fold change
  geom_text(aes(label = genefams2geneID.BG.GeneID), angle = 35, hjust = -0.2, size = 3) +  # Inclined GeneName labels
  labs(x = "Fold Change", y = "GO Category", color = "P-Value") +  # Axes labels and legend title
  theme_minimal() +  # Use a minimal theme (customize as needed)
  
  # Rotate y-axis labels for better readability if necessary
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  
  # Add dashed lines at x = -2 and x = 2
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black")

# Display the plot
print(gene_plot)


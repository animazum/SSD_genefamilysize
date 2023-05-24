
###### fold change Analysis ######
# Script created by Benjamin

rm(list=ls())
library(logr)

setwd("YOUR/DIRECTORY/foldChange/")

## positive or negative
testtype<- "negative"
## brain, f_brain or many others 
tissueOFinterest <- "brain"
## adult, fetus
mono<-"fetus"
### 1 up to how much time you want to lose. e.g.1000
perms= 1000

PGLSbackground<-read.csv("/PGLS.background.genefams2geneID.csv", row.names = 1)

expressioGenes<-read.csv("BrainSpan_SELECTED/GENES.csv")

tmp <- file.path(getwd(), paste("FoldchangeAnalysis.",mono,".",testtype,".associated.genes.",tissueOFinterest,".tissue.log", sep = ""))
lf <- log_open(tmp)

if (testtype == "negative" &&  mono == "adult") {
  log_print("Captain's log starts...")
  log_print("Using negatively selected genes", hide_notes = T)
  expressioGenes<-read.csv("BrainSpan_WholeBrainExpression.csv")
  ### remove the means with 0s, becase for the los2(fem/male) we would have some issues with infinite numbers
  expressioGenes<- expressioGenes[!expressioGenes$average.male == 0,]
  expressioGenes<- expressioGenes[!expressioGenes$average.female == 0,]
  background.expressioGenes<- expressioGenes[expressioGenes$Gene_stable_ID %in% PGLSbackground$Ensembl_gene_identifier,]
  
  NEGgenes<- read.csv("genefams2geneID.NEGATIVE.csv")
  NEGgenes<- NEGgenes[NEGgenes$Ensembl_gene_identifier %in% background.expressioGenes$Gene_stable_ID,]
  
  # fantom <- read.csv("ADULT_UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv") #Read in the fantom data
  log_print("Using Adult selected genes", hide_notes = T)
} else if (testtype == "negative" &&  mono == "fetus") {
  log_print("Captain's log starts...")
  log_print("Using negatively selected genes", hide_notes = T)
  expressioGenes<-read.csv("BrainSpan_WholeBrainExpression.csv")
  ### remove the means with 0s, becase for the los2(fem/male) we would have some issues with infinite numbers
  expressioGenes<- expressioGenes[!expressioGenes$average.male == 0,]
  expressioGenes<- expressioGenes[!expressioGenes$average.female == 0,]
  background.expressioGenes<- expressioGenes[expressioGenes$Gene_stable_ID %in% PGLSbackground$Ensembl_gene_identifier,]
  
  NEGgenes<- read.csv("genefams2geneID.NEGATIVE.csv")
  NEGgenes<- NEGgenes[NEGgenes$Ensembl_gene_identifier %in% background.expressioGenes$Gene_stable_ID,]
  
  log_print("Using fetous selected genes", hide_notes = T)
} else if (testtype == "positive" &&  mono == "adult") {
  log_print("Captain's log starts...")
  log_print("Using positive selected genes", hide_notes = T)
  expressioGenes<-read.csv("BrainSpan_WholeBrainExpression.csv")
  ### remove the means with 0s, becase for the los2(fem/male) we would have some issues with infinite numbers
  expressioGenes<- expressioGenes[!expressioGenes$average.male == 0,]
  expressioGenes<- expressioGenes[!expressioGenes$average.female == 0,]
  background.expressioGenes<- expressioGenes[expressioGenes$Gene_stable_ID %in% PGLSbackground$Ensembl_gene_identifier,]
  
  NEGgenes<- read.csv("genefams2geneID.POSITIVE.csv")
  NEGgenes<- NEGgenes[NEGgenes$Ensembl_gene_identifier %in% background.expressioGenes$Gene_stable_ID,]
  
  # fantom <- read.csv("ADULT_UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv") #Read in the fantom data
  log_print("Using Adult selected genes", hide_notes = T)
} else if (testtype == "positive" &&  mono == "fetus") {
  log_print("Captain's log starts...")
  log_print("Using positive selected genes", hide_notes = T)
  expressioGenes<-read.csv("BrainSpan_prenatalAges_WholeBrainExpression.csv")
  ### remove the means with 0s, becase for the log2(fem/male) we would have some issues with infinite numbers
  expressioGenes<- expressioGenes[!expressioGenes$average.male == 0,]
  expressioGenes<- expressioGenes[!expressioGenes$average.female == 0,]
  background.expressioGenes<- expressioGenes[expressioGenes$Gene_stable_ID %in% PGLSbackground$Ensembl_gene_identifier,]
  
  NEGgenes<- read.csv("genefams2geneID.POSITIVE.csv")
  NEGgenes<- NEGgenes[NEGgenes$Ensembl_gene_identifier %in% background.expressioGenes$Gene_stable_ID,]
  
  # fantom <- read.csv("ADULT_UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv") #Read in the fantom data
  log_print("Using fetous selected genes", hide_notes = T)
}

### creating the index column to extract the genes from the back ground and use them in the analysis
rownames(background.expressioGenes)<- background.expressioGenes$Gene_stable_ID
rownames(NEGgenes)<-NEGgenes$Ensembl_gene_identifier
NEGgenes<-background.expressioGenes[(rownames(background.expressioGenes)) %in% rownames(NEGgenes),]
number_of_genes <- dim(NEGgenes)[1]

#### permutations
fa1R<- data.frame(rownames(NEGgenes))
fa1RAvg<-data.frame()
BGlistR<-list()
for (i in 1:dim(background.expressioGenes)[1]) {
  
  if (i == dim(background.expressioGenes)[1]) {

    # log_print(cat("#############################\n##### average mean perms ####\n#############################"))
    log_print("#############################", hide_notes = T)
    log_print("##### average mean fold change perms ####", hide_notes = T)
    log_print("#############################", hide_notes = T)
    log_print("Done creating the Background fold change .___.", hide_notes = T)

    for (j in 1:perms) {
      Cerebro<-paste(tissueOFinterest,j, sep = "")
      die<-as.numeric(background.expressioGenes$fold.Change)
      fa1R[,Cerebro]<- as.numeric(sample(die, size = number_of_genes, replace = TRUE))
      if (j == perms) {
        
        fa1R[,1]<- NULL
        fa1RAvg<-data.frame(t(fa1R))
        fa1RAvg$totalfoldChangeAVG<-rowMeans(fa1RAvg)
        fa1RAvg<-subset(fa1RAvg, select = ("totalfoldChangeAVG"))
        log_print("Object fa1Ang and fa1R were created", hide_notes = F)

      }
    }
  }
}

#### real means
NEGgenesR<- mean(as.numeric(NEGgenes$fold.Change))



min(fa1RAvg)
max(fa1RAvg)
NEGgenesR
log_print(paste(tissueOFinterest," brain fold change has a mean bigger than random mean in :",
                data.frame(sum(NEGgenesR > fa1RAvg)), "/", dim(fa1RAvg)[1], sep = ""), hide_notes = T)
log_print(paste("Plot name: foldchangeAnalysis.",mono,".",testtype,".associated.genes.",tissueOFinterest,".tissue.pdf", sep = ""))

log_close()
writeLines(readLines(lf))

##### creating flot
densidad<-density(fa1RAvg$totalfoldChangeAVG)
pdf(paste("foldchangeAnalysis.",mono,".",testtype,".associated.genes.",tissueOFinterest,".tissue.pdf", sep = ""))
# Filled Density Plot
library(ggplot2)
p <- ggplot(fa1RAvg, aes(x=totalfoldChangeAVG)) + 
  geom_density() 
# Add mean line
p + scale_color_grey() + theme_classic() + geom_vline(aes(xintercept=NEGgenesR), 
                                                      color="red", linetype="dashed", size=1) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(densidad$y)+0.01)) + scale_x_continuous(expand = c(0, 0), limits = c(min(fa1RAvg)-0.2, max(fa1RAvg)+0.2))

dev.off()


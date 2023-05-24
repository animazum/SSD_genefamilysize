
###************************************************************************###
####  Script for tissue expression analysis using an average rank analysis. ####
###************************************************************************###

rm(list=ls())
library(logr)

setwd("YOUR/DATA/gene expression per tissue")

## positive or negative
testtype<- "negative"
## brain, f_brain or many others 
tissueOFinterest <- "f_brain"
## adult, fetus
mono<-"fetus"
### 1 up to how much time you want to lose. e.g.1000
perms= 1000

tmp <- file.path(getwd(), paste("RankingAnalysis.",mono,".",testtype,".associated.genes.",tissueOFinterest,".tissue.log", sep = ""))
lf <- log_open(tmp)

if (testtype == "negative") {
  log_print("Captain's log starts...")
  log_print("Using negatively selected genes", hide_notes = T)
  if (mono == "adult") {
    fantom <- read.csv("YOUR/FANTOM5/data.csv") #Read in the fantom data
    log_print("Using Adult selected genes", hide_notes = T)
  } else if (mono == "fetus") {
    fantomData<- read.csv("YOUR/FANTOM5/data.csv")
    feto<-grep(x = colnames(fantomData), pattern = "f_")
    fantom<- subset(fantomData, select = c(1,feto))
    log_print("Using fetous selected genes", hide_notes = T)
  }
  
  pgls_background <- read.csv("NEGATIVES/background_genes_for_SSD_fromPGLS")#Read in all the genes entered to PGLS
  
  pgls_assoc <- read.csv("NEGATIVES/associated_genes_for_SSD_fromPGLS")#Read in PGLS Associated genes
  
  pgls_notAssoc <- read.csv("NEGATIVES/associated_genes_for_SSD_fromPGLS")#Read in PGLS not associated genes
  
} else if (testtype == "positive") {
  log_print("Captain's log starts...")
  log_print("Using positively selected genes", hide_notes = T)
  if (mono == "adult") {
    fantom <- read.csv("YOUR/FANTOM5/data.csv") #Read in the fantom data
    log_print("Using Adult selected genes", hide_notes = T)
  } else if (mono == "fetus") {
    fantomData<- read.csv("YOUR/FANTOM5/data.csv")
    feto<-grep(x = colnames(fantomData), pattern = "f_")
    fantom<- subset(fantomData, select = c(1,feto))
    log_print("Using fetous selected genes", hide_notes = T)
  }
  
  pgls_background <- read.csv("POSITIVES/background_genes_for_SSD_fromPGLS")#Read in all the genes entered to PGLS
  
  pgls_assoc <- read.csv("POSITIVES/associated_genes_for_SSD_fromPGLS")#Read in PGLS Associated genes
  
  pgls_notAssoc <- read.csv("POSITIVES/not-associated_genes_for_SSD_fromPGLS")#Read in PGLS not associated genes
  
}

number_of_genes <- dim(pgls_assoc)[1] ##Change this
#Filter FANTOM data to PGLS background genes
fantom_background <- fantom[fantom$X %in% pgls_background$Gene.stable.ID ,]
fantom_background <- fantom_background[!duplicated(fantom_background$X) ,]
colnames(fantom_background)

row.names(fantom_background) <- fantom_background$X
fantom_background$X <- NULL

names <- row.names(fantom_background)
fantom_background$total<-(rowSums(fantom_background))
fantom_background<-fantom_background[!fantom_background$total == 0,]
fantom_background$total <-NULL

### creating the index column tio extract the genes from the back ground and use them in the analysis
rownames(pgls_assoc)<-pgls_assoc$Gene.stable.ID
NEGgenes<-fantom_background[(rownames(fantom_background)) %in% rownames(pgls_assoc),]
number_of_genes <- dim(NEGgenes)[1]

#### permutations
fa1R<- data.frame(rownames(NEGgenes))
fa1RAvg<-data.frame()
BGlistR<-list()
for (i in 1:dim(fantom_background)[1]) {
  
  genesR<- t(fantom_background[i,])
  GenenameR<- colnames(genesR)
  genesR[tissueOFinterest,]
  brainrank <- data.frame(sum(genesR[,1] < genesR[tissueOFinterest,]))
  rownames(brainrank)<- tissueOFinterest
  colnames(brainrank)<- GenenameR
  BGlistR[[i]]<- brainrank
  if (i == dim(fantom_background)[1]) {
    genesR<- do.call("cbind", BGlistR)
    genesRALLgenesvsBrain<<- t(genesR)
    genesR<<-rowMeans(genesR)
    
    # log_print(cat("#############################\n##### average mean perms ####\n#############################"))
    log_print("#############################", hide_notes = T)
    log_print("##### average mean perms ####", hide_notes = T)
    log_print("#############################", hide_notes = T)
    log_print("Done creating the Background means .___.", hide_notes = T)
    log_print("Object genesRALLgenesvsBrain and genesR were created", hide_notes = T)
    log_print(paste("The average rank for ", tissueOFinterest," is: ", genesR,sep = ""), hide_notes = T)
    # log_print("Object genesRALLgenesvsBrain and genesR were created")
    # log_print(paste("The average rank for ", tissueOFinterest," is: ", genesR, sep = ""))
    
    for (j in 1:perms) {
      Cerebro<-paste(tissueOFinterest,j, sep = "")
      # for (i in 1:dim(fantom_background)[1]) {
      # die<- c(0:(length(fantom_background)-1))
      die<-genesRALLgenesvsBrain[,tissueOFinterest]
      fa1R[,Cerebro]<- sample(die, size = number_of_genes, replace = TRUE)
      # }
      if (j == perms) {
        rownames(fa1R)<- fa1R[,1]
        fa1R[,1]<- NULL
        fa1RAvg<-data.frame(t(fa1R))
        fa1RAvg$totalBrainAVG<-rowMeans(fa1RAvg)
        fa1RAvg<-select(fa1RAvg, contains("totalBrainAVG"))
        log_print("Object fa1Ang and fa1R were created", hide_notes = F)
        # print("Object fa1Ang and fa1R were created", quote=F)
        
      }
    }
  }
}

#### real means
BGlist<- list()
for (i in 1:dim(NEGgenes)[1]) {
  NEGgenesR<- t(NEGgenes[i,])
  Genename<- colnames(NEGgenesR)
  NEGgenesR[tissueOFinterest,]
  brainrank <- data.frame(sum(NEGgenesR[,1] < NEGgenesR[tissueOFinterest,]))
  rownames(brainrank)<- tissueOFinterest
  colnames(brainrank)<- Genename
  BGlist[[i]]<- brainrank
  if (i == dim(NEGgenes)[1]) {
    NEGgenesR<- do.call("cbind", BGlist)
    NEGgenesALLgenesvsBrain<- t(NEGgenesR)
    NEGgenesR<-rowMeans(NEGgenesR)
    log_print("######################################", hide_notes = T)
    log_print("##### real means for focal tissue ####", hide_notes = T)
    log_print("######################################", hide_notes = T)
    log_print("Done creating the Real means .___.", hide_notes = T)
    log_print("Object NEGgenesALLgenesvsBrain and NEGgenesR were created", hide_notes = T)
    log_print(paste("The average rank for ", tissueOFinterest," is: ", NEGgenesR, sep = ""), hide_notes = T)
    # print("Done creating the Real means .___.", quote = F)
    # print("Object NEGgenesALLgenesvsBrain and NEGgenesR were created", quote =F)
    # print(paste("The average rank for ", tissueOFinterest," is: ", NEGgenesR, sep = ""), quote =F)
    
  }
}

min(fa1RAvg)
max(fa1RAvg)
NEGgenesR
log_print(paste(tissueOFinterest," tissue has an average rank bigger than random rank in :",
                data.frame(sum(NEGgenesR > fa1RAvg)), "/", dim(fa1RAvg)[1], sep = ""), hide_notes = T)
log_print(paste("Plot name: RankingAnalysis.",mono,".",testtype,".associated.genes.",tissueOFinterest,".tissue.pdf", sep = ""))

log_close()
writeLines(readLines(lf))

# hist(fa1RAvg$totalBrainAVG)

##### creating flot
densidad<-density(fa1RAvg$totalBrainAVG)
pdf(paste("RankingAnalysis.",mono,".",testtype,".associated.genes.",tissueOFinterest,".tissue.pdf", sep = ""))
# Filled Density Plot
p <- ggplot(fa1RAvg, aes(x=totalBrainAVG)) + 
  geom_density() 
# Add mean line
p + scale_color_grey() + theme_classic() + geom_vline(aes(xintercept=NEGgenesR), 
                                                      color="red", linetype="dashed", size=1) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(densidad$y)+0.01)) + scale_x_continuous(expand = c(0, 0), limits = c(min(fa1RAvg)-0.2, max(fa1RAvg)+0.2))

dev.off()


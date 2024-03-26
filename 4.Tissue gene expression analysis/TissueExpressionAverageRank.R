
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
    fantom <- read.csv("ADULT_UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv") #Read in the fantom data
    log_print("Using Adult selected genes", hide_notes = T)
  } else if (mono == "fetus") {
    fantomData<- read.csv("~/Dropbox/SSD/sex_biased_genexpression/phatom/UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv")
    feto<-grep(x = colnames(fantomData), pattern = "f_")
    fantom<- subset(fantomData, select = c(1,feto))
    log_print("Using fetous selected genes", hide_notes = T)
  }
  
  # pgls_background <- read.csv("NEGATIVES/8669background_genes_for_SSD_fromPGLS_SSD+BM.csv")#Read in all the genes entered to PGLS
  pgls_background <- read.csv("~/Dropbox/SSD/MANUSCRIPT/BG.DOS_negativePGLSgenes.csv")#Read in all the genes entered to PGLS
  
  pgls_assoc <- read.csv("NEGATIVES/738associated_genes_for_SSD_fromPGLS_SSD+BM.csv")#Read in PGLS Associated genes
  
  pgls_notAssoc <- read.csv("NEGATIVES/7931not-associated_genes_for_SSD_fromPGLS_SSD+BM.csv")#Read in PGLS not associated genes
  
} else if (testtype == "positive") {
  log_print("Captain's log starts...")
  log_print("Using positively selected genes", hide_notes = T)
  if (mono == "adult") {
    fantom <- read.csv("ADULT_UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv") #Read in the fantom data
    log_print("Using Adult selected genes", hide_notes = T)
  } else if (mono == "fetus") {
    fantomData<- read.csv("~/Dropbox/SSD/sex_biased_genexpression/phatom/UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv")
    feto<-grep(x = colnames(fantomData), pattern = "f_")
    fantom<- subset(fantomData, select = c(1,feto))
    log_print("Using fetous selected genes", hide_notes = T)
  }
  
  # pgls_background <- read.csv("POSITIVES/8669background_genes_for_SSD_fromPGLS_SSD+BM.csv")#Read in all the genes entered to PGLS
  pgls_background <- read.csv("~/Dropbox/SSD/MANUSCRIPT/BG.DOS_positivePGLSgenes.csv")#Read in all the genes entered to PGLS
  
  pgls_assoc <- read.csv("POSITIVES/588associated_genes_for_SSD_fromPGLS_SSD+BM.csv")#Read in PGLS Associated genes
  
  pgls_notAssoc <- read.csv("POSITIVES/8081not-associated_genes_for_SSD_fromPGLS_SSD+BM.csv")#Read in PGLS not associated genes
  
}

pgls_background<- data.frame(pgls_background[,3])
colnames(pgls_background)<- "Gene.stable.ID"

colnames(fantom)
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

### creating the index column to extract the genes from the back ground and use them in the analysis
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

##******* ************************************************###
##*  
# ... (Previous code remains unchanged until the permutation loop)

# Bootstrap Approach
B <- 5000  # Number of bootstrap samples

# Create an empty dataframe to store bootstrap results
bootstrap_results <- numeric(B)

fantom_background<-fantom_background[(rownames(fantom_background)) %in% rownames(pgls_assoc),]

for (b in 1:B) {
  # Generate a bootstrap sample (with replacement) from the observed data for the selected tissue column
  bootstrap_sample <- fantom_background[sample(nrow(fantom_background), replace = TRUE), tissueOFinterest]
  
  # Calculate the desired statistic or metric for the bootstrap sample (e.g., mean, median, sum)
  # For example, calculating the mean expression for the bootstrap sample
  bootstrap_sample_stat <- mean(bootstrap_sample)
  
  # Store the calculated statistic in the bootstrap_results vector
  bootstrap_results[b] <- bootstrap_sample_stat
}

# Calculate the desired statistic for the bootstrap results (e.g., mean, median)
bootstrap_statistic <- mean(bootstrap_results)  # Adjust this according to the desired metric

# ... (Further code to log results and create plots)

##******* ************************************************###

#### real means
BGlist<- list()
for (i in 1:dim(NEGgenes)[1]) {
  NEGgenesR<- t(NEGgenes[i,])
  ALLNEGgenesR<<-NEGgenesR
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
                data.frame(sum(NEGgenesR > fa1RAvg)), "/", dim(fa1RAvg)[1], ", p-value:",mean(fa1RAvg >= NEGgenesR),sep = ""), hide_notes = T)
log_print(paste("Plot name: RankingAnalysis.",mono,".",testtype,".associated.genes.",tissueOFinterest,".tissue.pdf", sep = ""))

log_close()
writeLines(readLines(lf))

##******* ************************************************###
# Compare observed means with bootstrap means
# Compare observed means with bootstrap means
p_value <- sum(bootstrap_results > NEGgenesR) / B # Calculate p-value
p_valueN <- sum(bootstrap_results < NEGgenesR) / B # Calculate p-value
cat("Empirical p-value:, right side", p_value, "\n")
cat("Empirical p-value, left side:", p_valueN, "\n")
# ... (Further code to log results and create plots)

##******* ************************************************###
# hist(fa1RAvg$totalBrainAVG)
bootstrap_results<- data.frame(bootstrap_results)
colnames(bootstrap_results)<- "totalBrainAVG"
##### creating flot

# Calculate density
densidad <- density(bootstrap_results$totalBrainAVG)

# Generate plot
p <- ggplot(bootstrap_results, aes(x = totalBrainAVG)) + 
  geom_density() +
  scale_color_grey() +
  theme_classic() +
  geom_vline(aes(xintercept = NEGgenesR), color = "red", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(densidad$y) + 0.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(min(bootstrap_results) - 20, max(bootstrap_results) + 10)) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
           label = paste("Empirical p-value, right side:", p_value)) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, 
           label = paste("Empirical p-value, left side:", p_valueN))

# Save plot to PDF
pdf(paste("REV6RankingAnalysis.", mono, ".", testtype, ".associated.genes.", tissueOFinterest, ".tissue.pdf", sep = ""))
print(p)
dev.off()

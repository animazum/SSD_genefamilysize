
## 03/24
#Script by Karina Diaz updated by Benjamin 
library(dplyr)
setwd("~/Dropbox/SSD/")
# setwd("YOUR/FILE/")
whereToSave<- getwd()
dir.create(file.path(whereToSave,'cellular_expression_Transformed'))
setwd(file.path(whereToSave,'cellular_expression_Transformed'))

rm(list=ls())
##Scatterplot for BrainSpan data. You can find the code to plot the mean expression or the Z-scores.

#load the data needed

expressioGenes<-read.csv("~/Dropbox/SSD/cellular_expression/Karina_brainspan/BrainSpan_Ages_MeanBrainStructure_WholeBrainExpression_Male_Female.csv")

###################################################################################################
######################### Zscore por etapa en las 4 listas distintas   ############################
###################################################################################################
#### ....................................
print("##########################")
print("######### Data in ########")
print("##########################")

### I used the same name from the GO enrichment output file to make the direct link between datasets
plotname<- "BothSEX_GO_SSD_Negative.GO.terms.table.Process.mammals.2Phenotype.fdr_fdr.BP.50.CorMammallistNoplatipus"

#####  Common backgound for the 4 list of genes
print("BG.DOS file should include all the output genes genes with R values calculated by your analysis. As the Z scores are calculated by the difference of means of your data and its standard deviation", quote =F)
print("BG.BOS file needs to have 2 columns, one with the Ensembl gene ID and the entrez gene ID", quote = F)

print("Creating the BG.DOS file...", quote=F)
## orthogroup or gene family data in linked to entrez IDs
NCBIgene2goFunc<-read.csv("~/Dropbox/SSD/PGLS/GO.terms.table.Process.mammals.csv", row.names = 1)
NCBIgene2goFunc<- NCBIgene2goFunc[,c(1,3)] %>% distinct()
## index file to link entrez id vs ensembl id
gene2ensembl<-data.table::fread("~/Dropbox/SSD/PGLS/gene2ensembl")
NCBIgene2goFunc<- merge(x = NCBIgene2goFunc, y = gene2ensembl[,c(2,3)], by.x = "GeneID", by.y = "GeneID")
# PGLS data in 
PGLS.result.1<-read.csv("~/Dropbox/SSD/PGLS.SSD_Av.Body.mass.124.Spp.5425.GeneFams.Rs_benjamini.2022-09-27.csv")

###******* Change me ************************************************###
BGRvalsCol <- 'R.t.value.SSD' ### the column were you R values are  ####
BGRvals <- 'Positive' #### can be Negative, Positive or ALL         ####
SEXO <- "ALL"         #### can be Male, Female or ALL               ####
###******************************* Change me ************************###

if (BGRvals == 'ALL') {
  PGLS.result.1 <- PGLS.result.1
  print("You are using a background set of ALL your data points", quote = F)
  plotname<- paste(BGRvals,'RsBackground.',plotname, sep = "")
  print(paste("Name in plot will be: ", plotname, sep = ""), quote = F)
} else if (BGRvals == 'Negative') {
  PGLS.result.1<- PGLS.result.1[PGLS.result.1[,BGRvalsCol] < 0, ]
  plotname<- paste(BGRvals,'RsBackground.',plotname, sep = "")
  print("You are using a background set of your NEGATIVE data points", quote = F)
  print(paste("Name in plot will be: ", plotname, sep = ""), quote = F)
} else if (BGRvals == 'Positive') {
  PGLS.result.1<- PGLS.result.1[PGLS.result.1[,BGRvalsCol] > 0, ]
  print("You are using a background set of your POSITIVE data points", quote = F)
  plotname<- paste(BGRvals,'RsBackground.',plotname, sep = "")
  print(paste("Name in plot will be: ", plotname, sep = ""), quote = F)
}

# creating table with ensembl id and entrez ids
BG.DOS<- merge(x = NCBIgene2goFunc, y = PGLS.result.1, by.x = "Orthogroup", by.y = "X")
BG.DOS <- as.data.frame(BG.DOS[,c(2:3)]) %>% distinct()
print(paste("Finished creating the BG.DOS file including ", dim(BG.DOS)[1], " genes from your analysis",sep = ""), quote=F)

##Keep only the genes of the Background from the BranSpan data.
DataBS<-expressioGenes[expressioGenes$Gene_stable_ID%in%BG.DOS$Ensembl_gene_identifier,] #14686   #### Modificacion 

print("Checking if you have duplicated genes")
if (sum(duplicated(DataBS$Gene_stable_ID)) == 0) {
  print("All good, you dont have any duplicated genes", quote = F)
} else { 
  print(paste("You have ",sum(duplicated(DataBS$Gene_stable_ID)),"duplicated genes"), quote = F)}

print(paste("you have ", dim(DataBS)[1], " genes with brain span data"), quote = F)

### ....................................
print('Creating table with data of the genes associated to longevity', quote =F)
print('Table should include be a csv file with 3 columns (entrez ID, Ensembl gene ID, R values)', quote = F)

if (BGRvals == 'Negative') {
  print('#******Negative selected Go categories *****###')
  # AGenes<- read.csv("/GO.ENRICHMENT.HM.SelectedGENES.csv", row.names = 1)
  AGenes<<- read.csv("~/Dropbox/SSD/PGLS.SSD_Av.Body.mass.124.Spp.5425.GeneFams.Rs_benjamini.2022-09-27/GO-HM.SelectedGENES.GO_SSD_Negative.GO.terms.table.Process.mammals.2Phenotype.fdr_fdr.BP.50All_Groups_BP50_fdr._Table.csv", row.names = 1)
} else if(BGRvals == 'Positive') {
  print('#******Positive selected Go categories *****###')
  # AGenes<- read.csv("/GO.ENRICHMENT.HM.SelectedGENES.csv", row.names = 1)
  AGenes<<- read.csv("~/Dropbox/SSD/PGLS.SSD_Av.Body.mass.124.Spp.5425.GeneFams.Rs_benjamini.2022-09-27/GO-HM.SelectedGENES.GO_SSD_Positive.GO.terms.table.Process.mammals.2Phenotype.fdr_fdr.BP.50All_Groups_BP50_fdr._Table.csv", row.names = 1)
}

AGenes<- AGenes[,c(2,25)]
AGenes<- as.data.frame(merge(x = gene2ensembl[,c(2,3)], y = AGenes ,by.x = "GeneID", by.y = "GeneID")) %>% distinct()
print(paste("Your AGenes object looks like this:"), quote = F)
print( head(AGenes, 3))

##Function to obtain the  Z-scores
FunZscore<-function(GENES, Rvalscol, PERMSnum){
  ##erase the id from the genes 
  print("Reading variables...")
  GENES<- GENES
  Rvalscol <- Rvalscol
  TablaLong<-DataBS[DataBS$Gene_stable_ID %in% GENES[,Rvalscol], -(1:2)]
  
  ## Apply log2 transformation to gene expression data
  TablaLong <- log2(TablaLong + 1)
  
  dim<-dim(TablaLong)[1]
  print("Filtering brain life span data vs your PGLS or analysis selected genes with r values")
  print(dim)
  
  print("calculationg mean of the data selected")
  meanLong<-apply(TablaLong, 2, mean)
  print(head(meanLong))
  meanRandom<-rep()
  
  print("Starting randomization", quote = F)
  print("Sampling the genes from DataBS by the number of your focal genes", quote =F)
  for (i in 1:PERMSnum){ ##change number
    #sample the genes and keep only the numeric info
    x<-DataBS[DataBS$Gene_stable_ID %in% (sample(DataBS$Gene_stable_ID, replace = FALSE, size = dim)), -(1:2)]
    x<-apply(x,2,mean)
    meanRandom<-rbind(meanRandom,x)
  }#for
  print("calculating the stats")
  meanRandom<-meanRandom[-1,]
  
  MeanOfMeans<<-apply(meanRandom, 2, mean)
  print("calc mean of random values")
  print(head(MeanOfMeans))
  
  sdRandom<-apply(meanRandom, 2, sd)
  print("SD or random values")
  print(head(sdRandom))
  
  zscore<-(meanLong-MeanOfMeans)/sdRandom
  print("calc. z scores")
  print(head(zscore))
  
  pval<-pnorm(-abs(zscore))
  print("calc. p vals")
  print(head(pval))
  
  temp<-rbind(meanLong,MeanOfMeans,sdRandom,zscore,pval)
  print("Final table")
  print(temp)
  
  return(temp)
}

#Run the function for every list of genes

Castillo<-FunZscore(GENES = AGenes, Rvalscol = "Ensembl_gene_identifier", PERMSnum = 1000)
Castillo[5,] < 0.05

Ages<-gsub(pattern = '.*?\\.(.*?)\\..*', replacement = '\\1', colnames(Castillo))
Ages <- sub("pcw", "/52", sub("mos", "/12", sub("yrs", "", Ages)))
convert_to_numeric <- function(x) {
  ifelse(grepl("/", x), {
    parts <- as.numeric(unlist(strsplit(x, "/")))
    parts[1] / parts[2]
  }, as.numeric(x))
}
# Apply the function to each string
Ages <- sapply(Ages, convert_to_numeric)

Sex<- gsub('.*\\.(.*)$', replacement = '\\1',colnames(Castillo))
Castillo <- rbind(Castillo, Ages)
Castillo <- rbind(Castillo, Sex)

if (SEXO == "ALL") {
  print("You are using female and male data")
} else if (SEXO == "Female") {
  print("You are using female data")
  # Find the column indices where the "SEX" row has "F"
  female_cols <- which(Castillo["Sex", ] == "F")
  Castillo <- Castillo[, female_cols]
} else if (SEXO == "Male") {
  print("You are using male data")
  # Find the column indices where the "SEX" row has "F"
  male_cols <- which(Castillo["Sex", ] == "M")
  Castillo <- Castillo[, male_cols]
}


write.csv(Castillo, paste("F.STATS.VALUES.",SEXO ,plotname,".csv", sep = ""))
write.csv(AGenes, paste("F.CellExp.GENES.", SEXO ,plotname, ".csv", sep = ""))


####################
######   PLOT ######
####################

#******************###
#Here I put as many options as I needed according to the numbers in my Y axis. 
#The function it is supposed to work with all data but if not, add or erase options according to your data and change the lenght of the Y axis.
#******************###

funPlot <- function(X, Y, COLOR, YLAB, MAIN) {
  # Calculate the range of Y values
  y_min <- min(Y)
  y_max <- max(Y)
  
  # Set the y-axis limits to include all data points
  y_range <- c(floor(y_min - abs(y_min * 0.2)), ceiling(y_max * 1.2))
  
  # Plot the data
  plot(x = log10(X), y = Y, type = "p", pch = 21, cex = 1.2, col = "black", main = MAIN, bg = COLOR,
       ylab = YLAB, xlab = "Time (years)", cex.lab = 1.2, axes = FALSE,
       cex.axis = 1.1, xlim = c(-1, 2), ylim = y_range)
  
  # Calculate the step size for the y-axis ticks
  y_step <- 1
  
  # Add y-axis ticks
  axis(2, at = seq(y_range[1], y_range[2], by = y_step), cex.axis = 1.1, tick = TRUE)
  
  # Add x-axis ticks
  axis(side = 1, at = c(-1, 0, 1, 2), labels = c("0.1", "1", "10", "100"), tick = TRUE, las = 1, cex.axis = 1.1, mgp = c(3, 1, 0))
  
  # Add vertical line
  abline(v = log10(0.75), col = "black", lty = "dashed", lwd = 1.2)
  
  # Add smooth line
  smooth_line <- smooth.spline(x = log10(X), y = Y, spar = 2)
  interpolated_y <- predict(smooth_line, x = log10(X))$y
  lines(log10(X), interpolated_y, lwd = 2.5, col = COLOR)
  
  # Calculate the slope of the smooth line
  slope <- coef(lm(interpolated_y ~ log10(X)))[[2]]
  
  # Display the slope value on the plot (at the top)
  text(x = max(log10(X)), y = max(Y), labels = paste("Slope:", round(slope, 2)), pos = 3)
}


###Plot and save plots as pdfs 
#Z scores
pdf(paste("F.Zscore_plot_age(years)_", SEXO, plotname,".pdf", sep = ""))
par(mfrow = c(2,2))
funPlot(X=as.numeric(Castillo[6,]), Y=as.numeric(Castillo[4,]), COLOR="green4", YLAB="Z score", MAIN="A")
dev.off()

##Expression
pdf(paste("F.Expression_plot_age(years)_", SEXO, plotname,".pdf", sep = ""))
par(mfrow = c(2,2))
funPlot(X=as.numeric(Castillo[6,]), Y=as.numeric(Castillo[1,]), COLOR="green4", YLAB="Expression", MAIN="A")
dev.off()




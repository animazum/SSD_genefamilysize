
##07/09/22
#Script by Karina Diaz updated by Benjamin 
library(dplyr)
setwd("YOUR/FILE/")
whereToSave<- getwd()
dir.create(file.path(whereToSave,'cellular_expression'))
setwd(file.path(whereToSave,'cellular_expression'))

rm(list=ls())
##Scatterplot for BrainSpan data. You can find the code to plot the mean expression or the Z-scores.

#load the data needed
#NewData1
load("~/Dropbox/SSD/cellular_expression/BrainSpan_Ages_MeanBrainStructure.RData")
#Df.edades
load("~/Dropbox/SSD/cellular_expression/Relacion_Id_Edades(years).RData")

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
NCBIgene2goFunc<-read.csv("GO.terms.table.Process.mammals.csv", row.names = 1)
NCBIgene2goFunc<- NCBIgene2goFunc[,c(1,3)] %>% distinct()
## index file to link entrez id vs ensembl id
gene2ensembl<-data.table::fread("~/Dropbox/SSD/PGLS/gene2ensembl")
NCBIgene2goFunc<- merge(x = NCBIgene2goFunc, y = gene2ensembl[,c(2,3)], by.x = "GeneID", by.y = "GeneID")
# PGLS data in 
PGLS.result.1<-read.csv("YOUR/PGLS/RESULT.csv")

###******* Change me ************************************************###
BGRvalsCol <- 'R.t.value.SSD' ### the column were you R values are  ####
BGRvals <- 'Negative' #### can be Negative, Positive or ALL         ####
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
DataBS<-NewData1[NewData1$Gene_stable_ID%in%BG.DOS$Ensembl_gene_identifier,] #14686   

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
  AGenes<- read.csv("/GO.ENRICHMENT.HM.SelectedGENES.csv", row.names = 1)
} else if(BGRvals == 'Positive') {
  print('#******Positive selected Go categories *****###')
  AGenes<- read.csv("/GO.ENRICHMENT.HM.SelectedGENES.csv", row.names = 1)
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
  TablaLong<-DataBS[DataBS[,1] %in% GENES[,Rvalscol],-(1:2)]
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
    x<-DataBS[DataBS[,1] %in% (sample(DataBS[,1],replace = F,size=dim)), -(1:2)]
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

Castillo<-FunZscore(GENES = AGenes, Rvalscol = "Ensembl_gene_identifier", PERMSnum = 10000)

write.csv(Castillo, paste("STATS.VALUES.",plotname,".csv", sep = ""))
write.csv(AGenes, paste("CellExp.GENES.", plotname, ".csv", sep = ""))
####################
######   PLOT ######
####################

#******************###
#Here I put as many options as I needed according to the numbers in my Y axis. 
#The function it is supposed to work with all data but if not, add or erase options according to your data and change the lenght of the Y axis.
#******************###

funPlot<-function(X,Y, COLOR, YLAB, MAIN){
  plot(x=log10(X),y=Y, type="p",pch=21, cex=1.2,  col="black", main=MAIN,bg=COLOR,
       ylab=YLAB,xlab="Time (years)", cex.lab=1.2,axes=F,
       cex.axis=1.1,xlim=c(-1,2), ylim = c(min(Y)-abs(min(Y*.2)),(max(Y)*1.2)))
  axis(2, at=seq(round(min(Y)-abs(min(Y*.2))) , round(max(Y*1.20)), 
                 by = round(abs((max(Y*1.2)) - (min(Y)-abs(min(Y*.2))))/5)), cex.axis=1.1, tick = TRUE)  ##y
  axis(side = 1, at = c(-1,0,1,2), labels = c("0.1","1","10","100"), tick = TRUE, las=1, cex.axis=1.1, mgp = c(3, 1, 0))##x
  abline(v=log10(.75),col="black",lty="dashed",lwd=1.2)##line
  lines(smooth.spline(x=log10(X),y=Y,spar=2),lwd=2.5,col=COLOR)##smooth ## Change to spar = 2 for straight line
}

###Plot and save plots as pdfs 

#Z scores
pdf(paste("3.Zscore_plot_age(years)_", plotname,".pdf", sep = ""))
par(mfrow = c(2,2))
funPlot(X=Df.edades$Age_years, Y=as.numeric(Castillo[4,]), COLOR="green4", YLAB="Z score", MAIN="A")
dev.off()

##Expression
pdf(paste("3.Expression_plot_age(years)_",plotname,".pdf", sep = ""))
par(mfrow = c(2,2))
funPlot(X=Df.edades$Age_years, Y=as.numeric(Castillo[1,]), COLOR="green4", YLAB="Expression", MAIN="A")
dev.off()






##### SELECTING GENES FOR TEMPORAL CELLULAR GENE EXPRESISION
library(dplyr)
#####  Common backgound for the 4 list of genes
print("BG.DOS file should include all the output genes genes with R values calculated by your analysis. As the Z scores are calculated by the difference of means of your data and its standard deviation", quote =F)
print("BG.BOS file needs to have 2 columns, one with the Ensembl gene ID and the entrez gene ID", quote = F)

print("Creating the BG.DOS file...", quote=F)
## orthogroup or gene family data in linked to entrez IDs
NCBIgene2goFunc<-read.csv("GO.terms.table.Process.mammals.csv", row.names = 1)
NCBIgene2goFunc<- NCBIgene2goFunc[,c(1,3)] %>% distinct()
## index file to link entrez id vs ensembl id
gene2ensembl<-data.table::fread("gene2ensembl")
NCBIgene2goFunc<- merge(x = NCBIgene2goFunc, y = gene2ensembl[,c(2,3)], by.x = "GeneID", by.y = "GeneID")

###******* Change me ************************************************###
PGLS.result.1<-read.csv("PGLS/RESULT.csv")
BGRvalsCol <- 'R.t.value.SSD' ### the column were you R values are  ####
BGpvalsCol <- 'p.adjusted.GFS.vs.SSD'
BGRvals <- 'Negative' #### can be Negative, Positive or ALL         ####
plotname<- ""
###******************************* Change me ************************###

if (BGRvals == 'ALL') {
  PGLS.result.1 <- PGLS.result.1
  print("You are using a background set of ALL your data points", quote = F)
  plotname<- paste(BGRvals,'RsBackground.',plotname, sep = "")
  print(paste("Name in plot will be: ", plotname, sep = ""), quote = F)
} else if (BGRvals == 'Negative') {
  PGLS.result.1<- PGLS.result.1[PGLS.result.1[,BGRvalsCol] < 0, ]
  PGLS.result.1<- PGLS.result.1[PGLS.result.1[,BGpvalsCol] < 0.05, ]
  plotname<- paste(BGRvals,'RsBackground.',plotname, sep = "")
  print("You are using a background set of your NEGATIVE data points", quote = F)
  print(paste("Name in plot will be: ", plotname, sep = ""), quote = F)
} else if (BGRvals == 'Positive') {
  PGLS.result.1<- PGLS.result.1[PGLS.result.1[,BGRvalsCol] > 0, ]
  PGLS.result.1<- PGLS.result.1[PGLS.result.1[,BGpvalsCol] < 0.05, ]
  print("You are using a background set of your POSITIVE data points", quote = F)
  plotname<- paste(BGRvals,'RsBackground.',plotname, sep = "")
  print(paste("Name in plot will be: ", plotname, sep = ""), quote = F)
}

# creating table with ensembl id and entrez ids
BG.DOS<- merge(x = NCBIgene2goFunc, y = PGLS.result.1, by.x = "Orthogroup", by.y = "X")
BG.DOS <- as.data.frame(BG.DOS[,c(2:3)]) %>% distinct()
print(paste("Finished creating the BG.DOS file including ", dim(BG.DOS)[1], " genes from your analysis",sep = ""), quote=F)

write.csv(BG.DOS, paste(BGRvals, "Genes.ForGeneExpression.csv", sep = ""))

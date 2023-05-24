

###******* ************************************************###
#### Creating the files for gene expression in tissue     ####
####    Creates the phantom table but with out            ####
####        multiple tissues with same name               ####
###********************************************************###


setwd(dir = "Your/FILE/sex_biased_genexpression/phatom/")
tissues<-read.csv("/FANTOM5/FILE/")

print("Creating the tissue table and index to map their ")
tissuescol<-colnames(tissues)
tissuescol<-gsub(pattern = "[.[:digit:]]", replacement = "", x = tissuescol)
tissuescol<- data.frame(table(tissuescol))
tissuesuni<- tissuescol[tissuescol$Freq == 1,]
tissuescol<- tissuescol[tissuescol$Freq >1,]


### average of multiple colums
tissues2<-list()
for (i in 1:dim(tissuescol)[1]) {
  termino<-paste("^",tissuescol[i,1], sep = "")
  tissues2[[i]] <- tissues[, grepl(termino, names(tissues))]
  tissues2[[i]] <- data.frame(rowMeans(tissues2[[i]]))
  colnames(tissues2[[i]])<- tissuescol[i,1]
}

tissues3<-tissues[,colnames(tissues) %in% tissuesuni$tissuescol]
tissues3$ref<-rownames(tissues3)
tissues4<-do.call("cbind", tissues2)
tissues4$ref<-rownames(tissues4)

tissues5<- merge(x = tissues3, y = tissues4, by = "ref")
tissues5$ref<-NULL
table(colnames(tissues5))


write.csv(tissues5, "UNIQUECOLs.rename.phatom data human.tissue.csv")




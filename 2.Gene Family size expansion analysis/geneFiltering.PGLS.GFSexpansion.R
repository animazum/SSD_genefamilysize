### Prepariong for PGLS 

## install older version of dplyr
# devtools::install_url("http://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.0.4.tar.gz")

## for pre pgls analyses
library("logr")
library("dplyr")
library("ggpubr")

## for pgls
library("ape")
library("nlme")
library("parallel")

## for the plot 
library("reshape")
library("ggplot2")

setwd("/Your/working/directory/")

print("change always the ulimit of the computer $ulimit -s 21000 !!!!!!!!!!!!")
rm(list=ls())

path1<- getwd()

### change here to set if you want to use MALE or FEMALE biased SSD
SSDbias <- "MALE"
### Change here depending in the phenotype of interest
trait<- "SSD"
### Change here to the type of correction that you want to use
Correcting<- "benjamini"

tmp <- file.path(getwd(), paste("GenefamilySize.",SSDbias,"biasSSD.",trait,".associated.genes.",Correcting,".correction.log", sep = ""))
lf <- log_open(tmp, traceback = F)

##### First we need to do the genome completion filtering. ####

## upload all the orthologs (gene families) per spp list. 
# Gene counts per orthogroup/Users/mona/Dropbox/SSD
# Orthogroups.GeneCount.tsv = orthofinder output
list_orthogroups_counts<- read.table(paste(path1,"/Orthogroups.GeneCount.tsv", sep = ''), header = T, row.names = 1)
log_print(paste("Species in the list:",dim(list_orthogroups_counts)[2]), hide_notes = T)
colnames(list_orthogroups_counts)<- gsub(pattern = ".LS", replacement = "", x = colnames(list_orthogroups_counts))

list_orthogroups_counts$Total<- NULL
### Single-copy genes, REMOVE!!! the last column that counts genes per orthogroup ###

## read the phenotype data file
Nombre.file<-c(paste(path1,"/SSD_mammalian_124sppWithOrtogroup_data_150721.csv", sep = ""))
traits<-read.csv(Nombre.file, stringsAsFactors = F)

### correct the name of some spp in the traits file
traits$Spp.Name<-gsub(pattern = "Balaenoptera_acutorostrata", replacement = "Balaenoptera_acutorostrata_scammoni", x = traits$Spp.Name)
traits$Spp.Name<-gsub(pattern = "Bison_bison", replacement = "Bison_bison_bison", x = traits$Spp.Name)
traits$Spp.Name<-gsub(pattern = "Ceratotherium_simum", replacement = "Ceratotherium_simum_simum", x = traits$Spp.Name)
traits$Spp.Name<-gsub(pattern = "Marmota_marmota", replacement = "Marmota_marmota_marmota", x = traits$Spp.Name)
traits$Spp.Name<-gsub(pattern = "Odocoileus_virginianus", replacement = "Odocoileus_virginianus_texanus", x = traits$Spp.Name)
traits$Spp.Name<-gsub(pattern = "Orycteropus_afer", replacement = "Orycteropus_afer_afer", x = traits$Spp.Name)
traits$Spp.Name<-gsub(pattern = "Saimiri_boliviensis", replacement = "Saimiri_boliviensis_boliviensis", x = traits$Spp.Name)
traits$Spp.Name<-gsub(pattern = "Gorilla_gorilla", replacement = "Gorilla_gorilla_gorilla", x = traits$Spp.Name)
traits<- traits[!traits$Spp.Name == "Ornithorhynchus_anatinus",] ### killing the Ornithorhynchus_anatinus

if (SSDbias == "FEMALE") {
  log_print("WE are removing all negative SSD values (Female biased SSD)")
  traits<-traits[traits$SSD >= 0,]
  log_print(sum(table(traits$SSD)))
} else if (SSDbias == "MALE") {
  log_print("We are not going to remove any SSD values (Male biased SSD)")
  log_print(table(traits$SSD))
  log_print(sum(table(traits$SSD)))
}

### read the tree that you used for orthofinder
tree <-read.tree(paste(path1,"/your/tree.tre", sep = ""))
## changing the name of the tree in case that you used the our "tutorial" in the Orthofinder section
tree$tip.label<-gsub(pattern = ".LS", replacement = "", x = tree$tip.label)

tree2<-drop.tip(tree,setdiff(tree$tip.label, traits$Spp.Name))

list_orthogroups_counts<-list_orthogroups_counts[, colnames(list_orthogroups_counts) %in% tree2$tip.label]

###############################################
# Open function for more detail of the process#
###############################################
Genome.completion.filtering<-function(FamData,Filt.Near.Uni, ssp.treshold){
  want.percent<- Filt.Near.Uni/100
  #percetout<<-Filt.Near.Uni
  #ssp.treshold.out<<- ssp.treshold
  ### Single-copy genes, REMOVE!!! the last column that counts genes per orthogroup ###
  NumCols1<-length(list_orthogroups_counts)
  
  # Filter the ortho counts by the single-copy orthologs. We are going to remove the orthogroups that have 0 or 1 number of genes per family 
  list1<- list_orthogroups_counts[1:NumCols1] %>% filter_all(all_vars(.< 2)) #  all orthogroups have 0 or 1 gene
  list1$Total<- rowSums(list1) # summatory of all the sigle-copy genes in new col. 
  ## Nearly universal genes calculation (filtering to get our core sigle-copy gene set)
  # NumCols1+1 (is the column with the total gene counts per orthogroup)
  
  # Why do I need to filter the orthogroups to keep the ones that have as many genes as the 90% or 95%  of the number of species?
  filtering.90.1<<- filter(list1, Total >= round(NumCols1*want.percent))    
  filtering.90.1<<- filter(list1, Total >= round(NumCols1*want.percent))    
  #filtering.90.1$Total<- NULL 
  # Notes:
  ## So we have 1339 orthogroups with at least a count of 117 genes (90% of the number of the spp).
  
  # keep the spp that have 90% of these single-copy genes orthogroups.
  # Find the treshold for the spp
  # All the spp most have at least the number of genes designated by the treshold. 
  treshold <- round(nrow(filtering.90.1)*(ssp.treshold/100))
  
  log_print(paste(round(nrow(filtering.90.1))))
  log_print(paste("All the spp most have at least the number of genes designated by the treshold:", treshold, "genes", sep = " "), hide_notes = T)
  
  # Count the number of single-copy genes per spp and remove those that does not have at least the number of genes dictated by the treshold.
  filtering.90.1<<- filtering.90.1[,colSums(filtering.90.1) >= treshold]
  filtering.90.1<- filtering.90.1[,colSums(filtering.90.1) >= treshold]
  # Note: apparently the 130 spp have at least 1205 single-copy genes. 
  
  log_print(paste( "Which spp or columns were removed? ", setdiff(colnames(list_orthogroups_counts),colnames(filtering.90.1))), hide_notes = T)
  
  ### Jackpot! the final table that will be used to run the PGLS with all the spp that we need. 
  gene_numbers<<-list_orthogroups_counts[colnames(list_orthogroups_counts) %in% colnames(filtering.90.1)]
  
  log_print(paste("You have", NumCols1, "spp that passed the genome completness filter", sep = " "),hide_notes = T )
  log_print("#### Jackpot! Now we have finished the filtering to know if the genomes used are complete ####", hide_notes = T)
  log_print("The dataframe gene_numbers was created", hide_notes = T)
  log_print("The filtering.90.1 object was created", hide_notes = T)
  write.csv(filtering.90.1, paste("nearly.universal.",Filt.Near.Uni,".percent.",ssp.treshold,"ssp.treshold",".csv", sep = ""))
  log_print(paste("File"," nearly.universal.",Filt.Near.Uni,".percent.",ssp.treshold,"ssp.treshold",".csv", " was created created", sep = ""))
  list.40_orthogroups_counts<-list_orthogroups_counts[colnames(list_orthogroups_counts) %in% colnames(filtering.90.1)]
  sppnumero<-dim(list.40_orthogroups_counts)[2]
  write.csv(list.40_orthogroups_counts, paste("Total.counts.",sppnumero,"spp.csv", sep = ""))
  log_print(paste("The file Total.counts.",sppnumero,"spp.csv was created", sep = ""))
}

#*************   function to filter orthogroups with gene counts   ********###
Genome.completion.filtering(FamData = list_orthogroups_counts, Filt.Near.Uni = 90, ssp.treshold = 90)
#************************###

colnames(x = traits)[5:8]

Trait.orthoGr.filtering<-function(traits, pheno.col.nums, spp.col.num, gene.numbers){
  phenotypes<- colnames(x = traits)[pheno.col.nums]
  for (i in 1:length(phenotypes)) {
    traits.filtered<<-traits[!is.na(traits[,phenotypes[i]]),]
  }
  spp.name<-colnames(x = traits)[spp.col.num]
  traits.filtered<<-traits.filtered[traits.filtered[,spp.name] %in% colnames(gene.numbers),]
  gene.numbers.filtered<-gene_numbers[,colnames(gene_numbers) %in% traits.filtered[,spp.name],]
  #gene.numbers.filtered<<-gene_numbers[colnames(gene_numbers) %in% traits.filtered$spp.name,]
  log_print("***********************************************************", hide_notes = T)
  log_print("****Spp missing in traits but present in gene.numbers:****", hide_notes = T)
  log_print("***********************************************************", hide_notes = T)
  log_print(setdiff(colnames(gene.numbers),traits.filtered[,spp.name]))
  log_print("***********************************************************", hide_notes = T)
  log_print("****Spp missing in gene.numbers but present in traits:****", hide_notes = T)
  log_print("***********************************************************", hide_notes = T)
  log_print(setdiff(traits$Spp.Name, traits.filtered$Spp.Name))
  if (dim(traits.filtered)[1] == 0) {
    log_print("The col. names from gene.number might be different to traits spp names", hide_notes = T)
  }
  log_print("****************************************************************************", hide_notes = T)
  log_print("Removing the orthogroups (gene families) that have 0 genes in 20% of the spp", hide_notes = T)
  log_print("****************************************************************************", hide_notes = T)
  # Keep only genes with a minimum number of missing values - in this case 20
  # As we have 40 species we want to keep rows with 20 or less zeros ??? ask araxi about this...
  ### change 20 to maximum zeros allowed e.g. 20/40 or 30/50 the diference should be 20...
  #Filter gene families to keep only those which have at least one gene in at least 6 species.
  
  #[huki] Remove any fams with more than 20% of zeros.
  zeros20perc <- ncol(gene.numbers.filtered)*0.2
  gene.numbers.filtered<-gene.numbers.filtered[rowSums(gene.numbers.filtered == 0) <= zeros20perc, ]
  
  log_print("************************************************************", hide_notes = T)
  log_print("Removing any orthogroupos (gene families) with variance of 0", hide_notes = T)
  log_print("************************************************************", hide_notes = T)
  #[huki]Remove any fams with variance of zeros.
  gene.numbers.filtered<-gene.numbers.filtered[apply(gene.numbers.filtered, 1, var) > 0,]
  
  #[huki]Remove any fams where the max number of genes in a row is 1
  # to avoid gene families present in only one sp. 
  log_print("****************************************************************************", hide_notes = T)
  log_print("Removing orthogroups (gene families) of 1 gene or less present in only 1 spp", hide_notes = T)
  log_print("*****************************************************************************", hide_notes = T)
  gene.numbers.filtered<<-gene.numbers.filtered[apply(gene.numbers.filtered, 1, max) > 2,]
  log_print(paste("The file ", paste("Filtered.GeneNum.",phenotypes,".csv"," was created", sep = "")))
  
  
  write.csv(gene.numbers.filtered, paste("Filtered.GeneNum.",phenotypes,".csv", sep = ""))
}

#*************   change the phenotype column number, the species list number   ********###
Trait.orthoGr.filtering(traits = traits, pheno.col.nums = 7, spp.col.num = 3, gene.numbers = gene_numbers)
#************************###

write.csv(gene.numbers.filtered, paste("Filtered.GeneNum.",dim(gene.numbers.filtered)[2],"spp","SSD",".csv", sep = ""))

for (i in 1:dim(traits.filtered)[1]) {
  traits.filtered$Av.Body.mass[i]<-as.numeric(log10((traits.filtered$Bodymass_male..g.[i] + traits.filtered$Bodymass_female..g.[i])/2))
}

### change this section based on your phenotype
traits.filtered[,"Av.Body.mass"]<- as.numeric(traits.filtered[,"Av.Body.mass"])

#*************   IT IS IMPORTANT TO HAVE THE SAME SPP ORDER BETWEEN TRAITS AND GENE NUMBERSS!!!!!!!   ********###
gene.numbers.filtered <- gene.numbers.filtered[,traits.filtered[,3]]
dim(traits.filtered)

# total of gene counts per spp
sumgenes<-rowSums(t(gene.numbers.filtered))


############################################
######## PGLS     PGLS      PGLS ###########
############################################

library(ape)
library(caper)
library(nlme)
library(parallel)

################
##### PGLS #####
### function ###
################

### calculation average bodymass, the results are in absolut values.

PGLS.GF.size<- function(traits, pheno.col.nums, tree, spp.col.num, gene.numbers, No.variables, where.Save.it){
  phenotypes<<- colnames(x = traits)[c(spp.col.num, pheno.col.nums)]
  log_print("#########################", hide_notes = T)
  log_print("##### PGLS PGLS PGLS ####", hide_notes = T)
  log_print("#########################", hide_notes = T)
  if (No.variables == 1){
    ###################################
    log_print("#### PGLS ONE Variable ######", hide_notes = T)
    ###################################
    out<<-as.data.frame(t(apply(X = gene.numbers,1,FUN = function(x){
      
      traits2<<-cbind(traits[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
      names(traits2)[dim(traits2)[2]]<<-"GFS"
      #next bit removes Genus.Species with NAs in events from tree, only necessary for losses, deletion, not for GFS
      tree<- drop.tip(tree,names(x)[is.na(x)]) 
      
      ## creating the formula for the model for 1 variable
      formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
      formilin2<-as.formula(paste("~",phenotypes[1], sep = ""))
      pglsModel<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      
      if (inherits(pglsModel, "try-error")) 
        return(c(NA,NA))
      else
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2], anova(pglsModel)$`p-value`[2], 
                 summary(pglsModel)$tTable[1,1], summary(pglsModel)$tTable[2,1], summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2], 
                 summary(pglsModel)$tTable[1,3], summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4]))
    })))
    colnames(out)<<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[2] , sep=" "), paste("p_value", phenotypes[2], sep=" "), 
                      "Coefficient Intercept", paste("Coefficient", phenotypes[2], sep=" "), "Coef SE Intercept", paste("Coef SE", phenotypes[2], sep=" "), 
                      "Coef tval Intercept", paste("Coef tval", phenotypes[2], sep=" "), "Coef pval Intercept", paste("Coef pval", phenotypes[2], sep=" "))
    
  } 
  else if (No.variables == 2) {
    ###################################
    log_print("#### PGLS TWO Variable ######", hide_notes = T)
    ###################################
    out<<-as.data.frame(t(apply(X = gene.numbers,MARGIN = 1,FUN = function(x){
      pheno1<<-phenotypes[3]
      pheno2<<-phenotypes[2]
      traits2<<-cbind(traits[,c("Spp.Name", pheno2, pheno1)],as.numeric(x)) #add phenotypes as needed
      # traits2<<-cbind(traits[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
      names(traits2)[dim(traits2)[2]]<<-"GFS"
      #next bit removes Genus.Species with NAs in events from tree, only necessary for losses, deletion, not for GFS
      tree<- drop.tip(tree,names(x)[is.na(x)])
      ## creating the formula for the model for 2 variables
      # formilin<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
      formilin<<-as.formula(paste("GFS", paste(phenotypes[c(3,2)], collapse = " + "), sep = " ~ "))
      #log_print((paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ ")))
      # formilin2<-as.formula(paste("~",colnames(traits)[spp.col.num], sep = ""))
      formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
      #log_print((paste("~",colnames(traits)[spp.col.num], sep = "")))
      # pglsModel<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      pglsModel<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      
      if (inherits(pglsModel, "try-error"))
        return(c(NA,NA))
      else
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2], anova(pglsModel)$`p-value`[2],
                 anova(pglsModel)$`F-value`[3], anova(pglsModel)$`p-value`[3], summary(pglsModel)$tTable[1,1],
                 summary(pglsModel)$tTable[2,1], summary(pglsModel)$tTable[3,1], summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2],
                 summary(pglsModel)$tTable[3,2], summary(pglsModel)$tTable[1,3], summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[3,3],
                 summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4], summary(pglsModel)$tTable[3,4]))
    })))
    colnames(out)<<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[3] , sep=" "), paste("p_value", phenotypes[3], sep=" "),
                      paste("F_value", phenotypes[2], sep=" "), paste("p_value", phenotypes[2], sep=" "), "Coefficient Intercept", paste("Coefficient", phenotypes[3], sep=" "),
                      paste("Coefficient", phenotypes[2], sep=" "),  "Coef SE Intercept", paste("Coef SE", phenotypes[3], sep=" "), paste("Coef SE", phenotypes[2], sep=" "),
                      "Coef tval Intercept", paste("Coef tval", phenotypes[3], sep=" "), paste("Coef tval", phenotypes[2], sep=" "), "Coef pval Intercept",
                      paste("Coef pval", phenotypes[3], sep=" "), paste("Coef pval", phenotypes[2], sep=" "))
  }
  else if (No.variables > 2){
    ###################################
    ## PGLS more than two Variables ###
    ###################################
    log_print("At the moment the function only works with upto 2 variables, sorry :P", hide_notes = T)
  }
  for (i in 1:No.variables) {
    ###########################
    ### R^2 calculation ###
    ###########################
    log_print("Calc. Rs", hide_notes = T)
    sps<-nrow(traits2)
    DeFr <- sps - (1 + No.variables)
    Coef.tval1<-toString(paste("Coef tval", phenotypes[i+1], sep=" "))
    R1 <- toString(paste("R t value", phenotypes[i+1], sep=" "))
    R_t_value1 <- out[[Coef.tval1]] / (sqrt(((out[[Coef.tval1]])^2) + DeFr)) # tvalue / square root of(tvalue^2 + DF))
    out[[R1]] <<- R_t_value1
    ### benjamini correction
    log_print("Calc. benjamini correction", hide_notes = T)
    out[[paste("p.adjusted GFS vs.", phenotypes[i+1], sep = "")]] <<- p.adjust(out[,paste("p_value", phenotypes[i+1], sep=" ")], method = "fdr")
  }
  
  nameout<- paste("PGLS", paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits2)[1], "Spp", dim(out)[1], "GeneFams", "Rs_benjamini", Sys.Date(), "csv", sep = ".")
  log_print(paste("Location and name of your file:", paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  log_print(paste("object created:", paste( "out", sep = "")), hide_notes = T)
  write.csv(out, paste(where.Save.it, nameout, sep = "/"))
}

### PGLS function
dirSave<-path1

log_print(paste("You have ",dim(gene.numbers.filtered)[1], " gene families to input the PGLS", sep = ""), hide_notes = T)

#*************   change the phenotype column number, the species list number, the number of phenotypes to use   ********###
PGLS.GF.size(traits = traits.filtered, pheno.col.nums = 7:8, No.variables = 2, 
             tree = tree2, spp.col.num = 3, gene.numbers = gene.numbers.filtered, 
             where.Save.it = dirSave)
#*********************###

###################################
#### PERMS PERMS PERMS PERMS ###### Permutations
###################################

PGLS.GF.size.perms<-function(traits.cl, pheno.col.nums, cores,tree.cl, spp.col.num, gene.numbers.cl, where.Save.it, No.perms){
  log_print("###########################", hide_notes = T)
  log_print("##### Perms perms PGLS ####", hide_notes = T)
  log_print("###########################", hide_notes = T)
  log_print("calc. cores", hide_notes = T)
  nucleos<-detectCores()
  nucleos<-nucleos + cores - nucleos
  cl <- makeCluster(nucleos) 
  log_print(paste(nucleos + cores - nucleos," cores detected and ready to roll"), hide_notes = T)
  
  #gene.numbers.cl<- gene.numbers.cl
  #traits.cl<- traits.cl
  #tree.cl<- tree.cl
  phenotypes<<- colnames(x = traits.cl)[c(spp.col.num, pheno.col.nums)]
  # traits.cl<<-traits.cl
  #No.perms <- No.perms
  
  ##### this way  clusterExports stops using .GlobalEnv and uses the eviropnment inside the function. 
  clusterExport(cl = cl , varlist = c("gene.numbers.cl", "traits.cl", "gls", "tree.cl", "corBrownian", "phenotypes", "No.perms"), envir = environment())
  
  log_print("Starting perms PGLS", hide_notes = T)
  
  if (length(pheno.col.nums) == 2) {
    perms<<- parLapply(cl, 1:No.perms, function(i,...){
      traits.cl<-traits.cl
      y<-sample(c(1:length(gene.numbers.cl)))
      gene_numbers2 <- gene.numbers.cl[,y]
      outperm1<-apply(gene_numbers2, MARGIN = 1,function(x){
        traits2<<-cbind(traits.cl[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c("Spp.Name" , "Av.Body.mass" , "SSD")],x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c(cat(paste0(dQuote(phenotypes[1:length(phenotypes)], F), collapse=" , ")))],x)) #add phenotypes as needed
        # The dQuote adds double quotes, the paste0 inserts the commas and cat shows the result without escaping special characters.
        #list2env(traits2, envir = .GlobalEnv)## Assign them to the global environment
        names(traits2)[dim(traits2)[2]]<-"GFS"
        
        # formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        formilin<<-as.formula(paste("GFS", paste(phenotypes[c(3,2)], collapse = " + "), sep = " ~ "))
        # formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        # formilin2<<-as.formula(paste("~",colnames(traits.cl)[spp.col.num], sep = ""))
        formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
        
        pglsModel<-try(gls(formilin, correlation = corBrownian(phy = tree.cl, form= formilin2), method = "ML", data = traits2))  #add phenotypes as needed
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2], 
                 anova(pglsModel)$`p-value`[2], anova(pglsModel)$`F-value`[3], anova(pglsModel)$`p-value`[3], 
                 summary(pglsModel)$tTable[1,1], summary(pglsModel)$tTable[2,1], summary(pglsModel)$tTable[3,1], 
                 summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2], summary(pglsModel)$tTable[3,2], 
                 summary(pglsModel)$tTable[1,3], summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[3,3], 
                 summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4], summary(pglsModel)$tTable[3,4]))
      })
      outperm1<-as.data.frame(t(outperm1))
      colnames(outperm1)<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[2] , sep=" "), paste("p_value", phenotypes[2], sep=" "), 
                            paste("F_value", phenotypes[3], sep=" "), paste("p_value", phenotypes[3], sep=" "), "Coefficient Intercept", 
                            paste("Coefficient", phenotypes[2], sep=" "), paste("Coefficient", phenotypes[3], sep=" "),  "Coef SE Intercept", 
                            paste("Coef SE", phenotypes[2], sep=" "), paste("Coef SE", phenotypes[3], sep=" "), "Coef tval Intercept", 
                            paste("Coef tval", phenotypes[2], sep=" "), paste("Coef tval", phenotypes[3], sep=" "), "Coef pval Intercept",
                            paste("Coef pval", phenotypes[2], sep=" "), paste("Coef pval", phenotypes[3], sep=" "))
      
      pglsList<- list(outperm1)
      return(outperm1)
      
    })
  } else if (length(pheno.col.nums) == 1) {
    perms<<- parLapply(cl, 1:No.perms, function(i,...){
      traits.cl<-traits.cl
      y<-sample(c(1:length(gene.numbers.cl)))
      gene_numbers2 <- gene.numbers.cl[,y]
      outperm1<-apply(gene_numbers2, MARGIN = 1,function(x){
        traits2<<-cbind(traits.cl[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c("Spp.Name" , "SSD" , "Av.Body.mass")],x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c(cat(paste0(dQuote(phenotypes[1:length(phenotypes)], F), collapse=" , ")))],x)) #add phenotypes as needed
        # The dQuote adds double quotes, the paste0 inserts the commas and cat shows the result without escaping special characters.
        #list2env(traits2, envir = .GlobalEnv)## Assign them to the global environment
        names(traits2)[dim(traits2)[2]]<-"GFS"
        
        formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        # formilin2<<-as.formula(paste("~",colnames(traits.cl)[spp.col.num], sep = ""))
        formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
        
        
        # pglsModel<-pgls(formilin, correlation = corBrownian(1, phy = tree.cl, form = formilin2), method = "ML", data = traits2) #add phenotypes as needed
        # pglsModel<-gls(formilin, correlation = corBrownian(1, phy = tree.cl, form = formilin2), method = "ML", data = traits2, ) #add phenotypes as needed
        
        pglsModel<-try(gls(formilin, correlation = corBrownian(phy = tree.cl, form= formilin2), method = "ML", data = traits2))  #add phenotypes as needed
        return(c(anova(pglsModel)$`F-value`[1], anova(pglsModel)$`p-value`[1], anova(pglsModel)$`F-value`[2],
                 anova(pglsModel)$`p-value`[2], summary(pglsModel)$tTable[1,1], summary(pglsModel)$tTable[2,1],
                 summary(pglsModel)$tTable[1,2], summary(pglsModel)$tTable[2,2], summary(pglsModel)$tTable[1,3],
                 summary(pglsModel)$tTable[2,3], summary(pglsModel)$tTable[1,4], summary(pglsModel)$tTable[2,4]))
      })
      outperm1<-as.data.frame(t(outperm1))
      colnames(outperm1)<-c("F_value Intercept", "p_value Intercept", paste("F_value", phenotypes[2] , sep=" "),
                            paste("p_value", phenotypes[2], sep=" "), "Coefficient Intercept",
                            paste("Coefficient", phenotypes[2], sep=" "), "Coef SE Intercept",
                            paste("Coef SE", phenotypes[2], sep=" "), "Coef tval Intercept",
                            paste("Coef tval", phenotypes[2], sep=" "), "Coef pval Intercept",
                            paste("Coef pval", phenotypes[2], sep=" "))
      
      pglsList<- list(outperm1)
      return(outperm1)
      
    })
  }
  
  log_print("Perms PGLS done", hide_notes = T)
  
  perms<<- as.data.frame(do.call(rbind, perms), rownames=TRUE)
  log_print(colnames(perms))
  
  for (i in 1:length(pheno.col.nums)) {
    
    log_print("creating variables for stats", hide_notes = T)
    
    #traits2<-as.data.frame(cbind(traits.cl[,c(phenotypes[1:length(phenotypes)])],x)) #add phenotypes as needed
    
    # perms<- bind_rows(perms,.id = NULL)
    sps<-nrow(traits.cl)
    DeFr<- nrow(traits.cl) - (1 + (length(phenotypes)-1)) ## the -1 removes the col.name Spp.Name that its not a phenotype.
    
    head(traits.cl)
    #### the 62 are the numbers of spp minus 2 = degrees of freedom
    log_print("tval and R values", hide_notes = T)
    
    Coef.tval1<<-toString(paste("Coef tval", phenotypes[i+1], sep=" "))
    
    head(perms)
    colnames(perms[[1]])
    rownames(perms)
    R_t_value <<- ((perms[,Coef.tval1]))/sqrt((((perms[,Coef.tval1])^2) + DeFr )) # square root of (tvalue^2 / (tvalue^2 + DF))
    log_print(head(R_t_value), hide_notes = T)
    
    R1 <- paste("R t value", phenotypes[i+1], sep=" ")
    perms[R1]<<- R_t_value
    log_print(" benjamini correction", hide_notes = T)
    
    ## benjamini correction
    perms[paste("p.adjusted Benjamini GFS vs ", phenotypes[i+1], sep = "")] <<- p.adjust(perms[,paste("p_value", phenotypes[i+1], sep=" ")], method = "fdr")
    #permsRs<<-permsRs
  }
  log_print("Creating output file", hide_notes = T)
  nameout<- paste("PGLS_perms", paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits.cl)[1], "Spp", dim(perms)[1], "Perm_GeneFams", "Rs_benjamini", Sys.Date(), "csv", sep = ".")
  log_print(paste("Location and name of your file:", paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  write.csv(perms, paste(where.Save.it, nameout, sep = "/"))
  perms<<- perms
  
  closeAllConnections()
}

#*************   change the phenotype column number, the species list number, the number of phenotypes to use   ********###
PGLS.GF.size.perms(traits.cl = traits.filtered, cores = 20, pheno.col.nums = 7:8, tree.cl = tree2, spp.col.num = 3, 
                   gene.numbers.cl = testgf, where.Save.it = dirSave, No.perms = 2)
#*********************###

colnames(perms)

head(perms,3)

#########################
##### plot plot plot ####
#########################
DensityPlots<- function(out, outcol ,perms , permscol, traits.cl, spp.col.num, pheno.col.nums, No.perms, where.Save.it, smoothing, Opacity, xaxisMin, XaxisMax, Phenos){
  
  log_print("#########################", hide_notes = T)
  log_print("##### plot plot plot ####", hide_notes = T)
  log_print("#########################", hide_notes = T)
  
  log_print("It is important to set genefamilies as rownames for out and perms objects", hide_notes = T)  
  
  phenotypes<<- colnames(x = traits.cl)[c(spp.col.num, pheno.col.nums)]
  log_print("Phenos in...", hide_notes = T)
  
  R1<-paste("R.t.value", phenotypes[1+1], sep=".")
  P1<-paste("R.t.value", phenotypes[1+1], sep=".")
  
  log_print("Change R1 if you there is an issue here", hide_notes = T)
  log_print(R1, hide_notes = T)
  log_print(P1, hide_notes = T)
  
  C1<-data.frame(rownames(out),out[,outcol])
  C2<-data.frame(rownames(perms),perms[,permscol])
  
  dat<-merge(x = C1, y = C2, by.x = "rownames.out.", by.y = "rownames.perms.", all.y = T)
  dat$rownames.out.<-NULL
  
  log_print("Setting columns [][]", hide_notes = T)
  #for two variables pĺot
  colnames(dat) <- c(paste("Real R²", phenotypes[1+1], sep=" "),paste("Perms R²", phenotypes[1+1], sep=" "))
  dat2<-melt(dat)
  
  log_print("Name of the plot", hide_notes = T)

    nameout.plot<- paste("PGLS_perms", Phenos,"Phenotype",No.perms,"perms",paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits.cl)[1], "Spp", dim(perms)[1], "Perm_GeneFams", "Rs_benjamini", Sys.Date(), "pdf", sep = ".")
  log_print(nameout.plot, hide_notes = T)
  pdf(paste(where.Save.it, nameout.plot, sep = "/")) ### not yet possible to save it from

  print("Building plot")
  plot1<-ggplot(dat2, aes(x = value, fill = variable)) + geom_density(alpha= Opacity, na.rm = T, adjust = smoothing, ) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) + xlim(xaxisMin, XaxisMax)
  
  print("Saving plot")
  print(plot1)
  dev.off()
}

DensityPlots(out = out, outcol = "R t value SSD", perms = perms, permscol= "R t value SSD", 
             traits.cl = traits.filtered, spp.col.num = 3, pheno.col.nums = 7, 
             No.perms = 1000, Opacity = 0.1 , smoothing = 1.5, xaxisMin= -0.8, XaxisMax = 0.8, 
             Phenos = "SSD", where.Save.it = dirSave)

closeAllConnections()

colnames(out)

log_close()
writeLines(readLines(lf))

##################################
################################## plot of the perms vs real Rs
library("reshape")
library("ggplot2")

# rm(list=ls())

#This is the file with wich you want to compare to the infile
R1<-paste("R t value", phenotypes[1+1], sep=" ")
#permsfile<-toString(paste("R t value", phenotypes[1+1], sep=" "))
dim(out)
dim(perms)

dat<-data.frame(out, perms)
#for two variables pĺot
colnames(dat) <- c(paste("Real R²", phenotypes[1+1], sep=" "),paste("Perms R²", phenotypes[1+1], sep=" "))

dat2<-melt(dat)


pdf("hist1000_Longevity.pdf")
ggplot(dat2, aes(x=value, fill = variable)) +
  geom_density(alpha=0.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

cor(x = traits.filtered$SSD, y = traits.filtered$Av.Body.mass)
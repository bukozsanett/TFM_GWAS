
setwd("/Volumes/ws_helix/HELIX_analyses/GWAS_chemical_exposure/db")

library(manta)
library(devtools)


Phtalates<-read.csv("HELIX_phtalates_cov_plink.txt", sep="\t",header = TRUE)


#transform to a data frame
Phtalates<-as.data.frame(Phtalates)

head(Phtalates)

#get the chemicals in a separate data frame and concert to a matrix for manta
chem<-Phtalates[,39:43]

write.table(chem, file = "chem.txt", sep = "\t")

chem<-as.matrix(Phtalates[,39:43])

#get the variables I need for manta

var<-data.frame("PC1"=Phtalates$PC1,"PC2"=Phtalates$PC2,"PC3"=Phtalates$PC3,"PC4"=Phtalates$PC4,"PC5"=Phtalates$PC5,"PC6"=Phtalates$PC6,"PC7"=Phtalates$PC7,"PC8"=Phtalates$PC8,"PC9"=Phtalates$PC9,"PC10"=Phtalates$PC10,"sex"=Phtalates$e3_sex,"age"=Phtalates$hs_child_age_days_None)

write.table(var, file = "covar.txt", sep = "\t")

manta(chem~.,data = var)


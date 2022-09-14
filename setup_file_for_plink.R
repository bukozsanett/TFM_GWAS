setwd("/Volumes/ws_helix")
getwd()
#Reading data 

#reading covariates and checking out the dataset
meta_EU<-read.table("/Volumes/ws_helix/HELIX_preproc/gwas/Final_data/HELIX_GWAS_FINAL_meta_updPCsEUR.txt",sep="\t", header=T)
head(meta_EU)
#check the dimension
dim(meta_EU)
#1155   36


#Let's read the phtalates cov

phtalates<-read.csv("/Volumes/ws_helix/HELIX_analyses/GWAS_chemical_exposure/db/HELIX_phtalates_cov_20220330.txt",sep="\t", header=T)
head(phtalates)
dim(phtalates)
#1301   92

table(is.na(phtalates$hs_mehhp_cadj_Log2))


#Let's calculate the reatios and add to the phtalate data.frame

#MEHHP/MEHP
#hs_mehhp_cadj_Log2/hs_mehp_cadj_Log2

phtalates$MEHHP_MEHP<-as.numeric(as.numeric(phtalates$hs_mehhp_cadj_Log2)/as.numeric(phtalates$hs_mehp_cadj_Log2))

table(is.na(phtalates$MEHHP_MEHP))


#MEOHP/MEHP 
#hs_meohp_cadj_Log2/hs_mehp_cadj_Log2

phtalates$MEOHP_MEHP<-as.numeric(as.numeric(phtalates$hs_meohp_cadj_Log2)/as.numeric(phtalates$hs_mehp_cadj_Log2))

table(is.na(phtalates$MEOHP_MEHP))
#MECPP/MEHP

#hs_mecpp_cadj_Log2/hs_mehp_cadj_Log2

phtalates$MECPP_MEHP<-as.numeric(as.numeric(phtalates$hs_mecpp_cadj_Log2)/as.numeric(phtalates$hs_mehp_cadj_Log2))


#check if we actually added to our dataset

head(phtalates)
dim(phtalates)

#1301   95


#merge the phtalates and the meta_EU dataframe together

df<-merge(meta_EU,phtalates,by.x="HelixID")
head(df)



#for the plink we need the FID first and after the IID so we have to change the order
#FID=IID=HelixID

colnames(df)[1:3]<-c("FID","IID","HelixID")
head(df)


#49-0,
#with plink to categorical variable can not be NA we have to change it to NONE
df[50][is.na(df[50])] = "NONE"
df[78][is.na(df[78])] = "NONE"
df[79][is.na(df[79])] = "NONE"
df[83][is.na(df[83])] = "NONE"
df[84][is.na(df[84])] = "NONE"
df[127][is.na(df[127])] = "NONE"


#let's save it in a txt

write.table(df,"/Volumes/ws_helix/HELIX_analyses/GWAS_chemical_exposure/db/HELIX_phtalates_cov_plink.txt",quote=F, row.names=F, col.names = T,sep="\t")

#continue in bash and run the plink

eu<-read.csv("/Volumes/ws_helix/HELIX_analyses/GWAS_chemical_exposure/db/HELIX_phtalates_cov_plink.txt",sep="\t", header=T)

table(is.na(eu))

head(eu)

class(eu$hs_mecpp_cadj_Log2)
table(is.na(eu$hs_mecpp_cadj_Log2))


class(eu$hs_mehp_cadj_Log2)
table(is.na(eu$hs_mehp_cadj_Log2))


class(eu$MEHHP_MEHP)

table(is.na(eu$MEHHP_MEHP))


dim(eu)

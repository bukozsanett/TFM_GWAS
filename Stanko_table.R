setwd("/Volumes/ws_helix/HELIX_analyses/GWAS_chemical_exposure/results/phthalates/MEHP/data")
getwd()
setwd("/Users/zsanika/Documents/Omics_Data_Analyst/Final_Master_Project/PHTALATES/plink_data")
library(data.table)

# MEHP
dd<-fread("MEHP.hs_mehp_cadj_Log2.glm.linear")
dim(dd)
head(dd)

#CYP2C9
dd1<-dd[dd$ID=="rs1057910",]
dd1
dd1<-dd[dd$ID=="rs1799853",]
dd1

#CYP2C19
dd1<-dd[dd$ID=="rs4244285",]
dd1
dd1<-dd[dd$ID=="rs12248560",]
dd1

#CYP2D6
dd1<-dd[dd$ID=="rs3892097",]
dd1

#UGT2B15
dd1<-dd[dd$ID=="rs1902023",]
dd1

#UGT1A7	
dd1<-dd[dd$ID=="rs11692021",]
dd1

# MEHHP
cc<-fread("MEHHP.hs_mehhp_cadj_Log2.glm.linear")
dim(cc)
head(cc)

#CYP2C9
cc1<-cc[cc$ID=="rs1057910",]
cc1
cc1<-cc[cc$ID=="rs1799853",]
cc1

#CYP2C19
dd1<-cc[cc$ID=="rs4244285",]
dd1
dd1<-cc[cc$ID=="rs12248560",]
dd1

#CYP2D6
dd1<-cc[cc$ID=="rs3892097",]
dd1

#UGT2B15
dd1<-cc[cc$ID=="rs1902023",]
dd1

#UGT1A7	
dd1<-cc[cc$ID=="rs11692021",]
dd1



#####MEOHP

dd<-fread("MEoHP.hs_meohp_cadj_Log2.glm.linear")
dim(dd)
head(dd)

#CYP2C9
dd1<-dd[dd$ID=="rs1057910",]
dd1
dd1<-dd[dd$ID=="rs1799853",]
dd1

#CYP2C19
dd1<-dd[dd$ID=="rs4244285",]
dd1
dd1<-dd[dd$ID=="rs12248560",]
dd1

#CYP2D6
dd1<-dd[dd$ID=="rs3892097",]
dd1

#UGT2B15
dd1<-dd[dd$ID=="rs1902023",]
dd1

#UGT1A7	
dd1<-dd[dd$ID=="rs11692021",]
dd1


##MECPP

dd<-fread("MECPP.hs_mecpp_cadj_Log2.glm.linear")
dim(dd)
head(dd)

#CYP2C9
dd1<-dd[dd$ID=="rs1057910",]
dd1
dd1<-dd[dd$ID=="rs1799853",]
dd1

#CYP2C19
dd1<-dd[dd$ID=="rs4244285",]
dd1
dd1<-dd[dd$ID=="rs12248560",]
dd1

#CYP2D6
dd1<-dd[dd$ID=="rs3892097",]
dd1

#UGT2B15
dd1<-dd[dd$ID=="rs1902023",]
dd1

#UGT1A7	
dd1<-dd[dd$ID=="rs11692021",]
dd1


###MEOHP/MEHP

dd<-fread("MEOHP_MEHP.MEOHP_MEHP.glm.linear")
dim(dd)
head(dd)

#CYP2C9
dd1<-dd[dd$ID=="rs1057910",]
dd1
dd1<-dd[dd$ID=="rs1799853",]
dd1

#CYP2C19
dd1<-dd[dd$ID=="rs4244285",]
dd1
dd1<-dd[dd$ID=="rs12248560",]
dd1

#CYP2D6
dd1<-dd[dd$ID=="rs3892097",]
dd1

#UGT2B15
dd1<-dd[dd$ID=="rs1902023",]
dd1

#UGT1A7	
dd1<-dd[dd$ID=="rs11692021",]
dd1

#MECPP-MEHP

dd<-fread("MECPP_MEHP.MECPP_MEHP.glm.linear")
dim(dd)
head(dd)

#CYP2C9
dd1<-dd[dd$ID=="rs1057910",]
dd1
dd1<-dd[dd$ID=="rs1799853",]
dd1

#CYP2C19
dd1<-dd[dd$ID=="rs4244285",]
dd1
dd1<-dd[dd$ID=="rs12248560",]
dd1

#CYP2D6
dd1<-dd[dd$ID=="rs3892097",]
dd1

#UGT2B15
dd1<-dd[dd$ID=="rs1902023",]
dd1

#UGT1A7	
dd1<-dd[dd$ID=="rs11692021",]
dd1

getwd()



#MBZP

dd<-fread("MBZP.hs_mbzp_cadj_Log2.glm.linear")
dim(dd)
head(dd)

#CYP2C9
dd1<-dd[dd$ID=="rs1057910",]
dd1
dd1<-dd[dd$ID=="rs1799853",]
dd1

#CYP2C19
dd1<-dd[dd$ID=="rs4244285",]
dd1
dd1<-dd[dd$ID=="rs12248560",]
dd1

#CYP2D6
dd1<-dd[dd$ID=="rs3892097",]
dd1

#UGT2B15
dd1<-dd[dd$ID=="rs1902023",]
dd1

#UGT1A7	
dd1<-dd[dd$ID=="rs11692021",]
dd1

#MIBP

dd<-fread("MIBP.hs_mibp_cadj_Log2.glm.linear")
dim(dd)
head(dd)

#CYP2C9
dd1<-dd[dd$ID=="rs1057910",]
dd1
dd1<-dd[dd$ID=="rs1799853",]
dd1

#CYP2C19
dd1<-dd[dd$ID=="rs4244285",]
dd1
dd1<-dd[dd$ID=="rs12248560",]
dd1

#CYP2D6
dd1<-dd[dd$ID=="rs3892097",]
dd1

#UGT2B15
dd1<-dd[dd$ID=="rs1902023",]
dd1

#UGT1A7	
dd1<-dd[dd$ID=="rs11692021",]
dd1


#MNBP

dd<-fread("MNBP.hs_mnbp_cadj_Log2.glm.linear")
dim(dd)
head(dd)

#CYP2C9
dd1<-dd[dd$ID=="rs1057910",]
dd1
dd1<-dd[dd$ID=="rs1799853",]
dd1

#CYP2C19
dd1<-dd[dd$ID=="rs4244285",]
dd1
dd1<-dd[dd$ID=="rs12248560",]
dd1

#CYP2D6
dd1<-dd[dd$ID=="rs3892097",]
dd1

#UGT2B15
dd1<-dd[dd$ID=="rs1902023",]
dd1

#UGT1A7	
dd1<-dd[dd$ID=="rs11692021",]
dd1


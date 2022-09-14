setwd("/Volumes/ws_helix/HELIX_analyses/GWAS_chemical_exposure/results/phthalates/MECPP/data")
library(data.table)
library(qqman)
library(GWASTools)
library("xlsx")

#MEHP
MEHP<-fread("MEHP.hs_mehp_cadj_Log2.glm.linear")

#choose the right test
MEHP<-MEHP[MEHP$TEST=="ADD",]
head(MEHP)

dim(MEHP)

summary(MEHP$P)

#convert to a data.frame
MEHP<-as.data.frame(MEHP)

#Manhattan plot with all data

png("Manhattan_MEHP.png")
manhattanPlot(MEHP$P,MEHP$`#CHROM`,main="MEHP",signif = 5e-8)
dev.off()

#calculate lambda
pvalue<-MEHP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MEHP.png")
qq(MEHP$P,main="MEHP lambda=1.000107")
dev.off()

#let's use a a limit P<1e-8

MEHP_bonf<-subset(MEHP, P<1e-8)

head(MEHP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MEHP_04<-subset(MEHP, P<1e-4)

dim(MEHP_04)

#Manhattan plot with the subsetted data

png("Manhattan_MEHP_04.png")
manhattanPlot(MEHP_04$P,MEHP_04$`#CHROM`,main="MEHP limit P<1e-4",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEHP_04$P)

#qq plot with all samples
png("qq_MEHP_04.png")
qq(MEHP_04$P,main="MEHP limit P<1e-4,lambda=1.390586e-04")
dev.off()

#let's use an other limit P<1e-5

MEHP_05<-subset(MEHP, P<1e-5)

dim(MEHP_05)

#Manhattan plot with the subsetted data

png("Manhattan_MEHP_05.png")
manhattanPlot(MEHP_05$P,MEHP_05$`#CHROM`,main="MEHP limit P<1e-5",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEHP_05$P)

#qq plot with all samples
png("qq_MEHP_05.png")
qq(MEHP_04$P,main="MEHP limit P<1e-5,lambda=1.578456e-05")
dev.off()

#get an excel file with the top SNPs

MEHP_topSNPs<-subset(MEHP, P<1e-5)

write.xlsx(MEHP_topSNPs, file = "MEHP_topSNPs.xlsx")

##############MEHHP###########################

MEHHP<-fread("MEHHP.hs_mehhp_cadj_Log2.glm.linear")

#choose the right test
MEHHP<-MEHHP[MEHHP$TEST=="ADD",]
head(MEHHP)

dim(MEHHP)

summary(MEHHP$P)

#convert to a data.frame
MEHHP<-as.data.frame(MEHHP)

#Manhattan plot with all data

png("Manhattan_MEHHP.png")
manhattanPlot(MEHHP$P,MEHHP$`#CHROM`,main="MEHHP",signif = 5e-8)
dev.off()

#calculate lambda

pvalue<-MEHHP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MEHHP.png")
qq(MEHHP$P,main="MEHHP lambda=1.015766")
dev.off()

#let's use a a limit P<1e-8

MEHHP_bonf<-subset(MEHHP, P<1e-8)

head(MEHHP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MEHHP_04<-subset(MEHHP, P<1e-4)

dim(MEHHP_04)

#Manhattan plot with the subsetted data

png("Manhattan_MEHHP_04.png")
manhattanPlot(MEHHP_04$P,MEHHP_04$`#CHROM`,main="MEHHP limit P<1e-4",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEHHP_04$P)

#qq plot with all samples
png("qq_MEHHP_04.png")
qq(MEHHP_04$P,main="MEHHP limit P<1e-4,lambda=1.307197e-04")
dev.off()

#let's use an other limit P<1e-5

MEHHP_05<-subset(MEHHP, P<1e-5)

dim(MEHHP_05)

#Manhattan plot with the subsetted data

png("Manhattan_MEHHP_05.png")
manhattanPlot(MEHHP_05$P,MEHHP_05$`#CHROM`,main="MEHHP limit P<1e-5",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEHHP_05$P)

#qq plot with all samples
png("qq_MEHHP_05.png")
qq(MEHHP_05$P,main="MEHHP limit P<1e-5,lambda=2.83327e-05")
dev.off()

#get an excel file with the top SNPs

MEHHP_topSNPs<-subset(MEHHP, P<1e-5)

write.xlsx(MEHHP_topSNPs, file = "MEHHP_topSNPs.xlsx")


################MEOHP

MEOHP<-fread("MEoHP.hs_meohp_cadj_Log2.glm.linear")

#choose the right test
MEOHP<-MEOHP[MEOHP$TEST=="ADD",]
head(MEOHP)

dim(MEOHP)

summary(MEOHP$P)

#convert to a data.frame
MEOHP<-as.data.frame(MEOHP)

#Manhattan plot with all data

png("Manhattan_MEOHP.png")
manhattanPlot(MEOHP$P,MEOHP$`#CHROM`,main="MEOHP",signif = 5e-8)
dev.off()

#calculate lambda
pvalue<-MEOHP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MEOHP.png")
qq(MEOHP$P,main="MEOHP lambda=1.013677")
dev.off()

#let's use a a limit P<1e-8

MEOHP_bonf<-subset(MEOHP, P<1e-8)

head(MEOHP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MEOHP_04<-subset(MEOHP, P<1e-4)

dim(MEOHP_04)

#Manhattan plot with the subsetted data

png("Manhattan_MEOHP_04.png")
manhattanPlot(MEOHP_04$P,MEOHP_04$`#CHROM`,main="MEOHP limit P<1e-4",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEOHP_04$P)

#qq plot with all samples
png("qq_MEOHP_04.png")
qq(MEOHP_04$P,main="MEOHP limit P<1e-4,lambda=1.361906e-04")
dev.off()

#let's use an other limit P<1e-5

MEOHP_05<-subset(MEOHP, P<1e-5)

dim(MEOHP_05)

#Manhattan plot with the subsetted data

png("Manhattan_MEOHP_05.png")
manhattanPlot(MEOHP_05$P,MEOHP_05$`#CHROM`,main="MEOHP limit P<1e-5",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEOHP_05$P)

#qq plot with all samples
png("qq_MEOHP_05.png")
qq(MEOHP_05$P,main="MEHHP limit P<1e-5,lambda=2.304437e-05")
dev.off()

#get an excel file with the top SNPs

MEOHP_topSNPs<-subset(MEOHP, P<1e-5)

write.xlsx(MEOHP_topSNPs, file = "MEOHP_topSNPs.xlsx")

###################MEHHP_MEHP



MEHHP_MEHP<-fread("MEHHP_MEHP.MEHHP_MEHP.glm.linear")

#choose the right test
MEHHP_MEHP<-MEHHP_MEHP[MEHHP_MEHP$TEST=="ADD",]
head(MEHHP_MEHP)

dim(MEHHP_MEHP)

summary(MEHHP_MEHP$P)

#convert to a data.frame
MEHHP_MEHP<-as.data.frame(MEHHP_MEHP)

#Manhattan plot with all data

png("Manhattan_MEHHP_MEHP.png")
manhattanPlot(MEHHP_MEHP$P,MEHHP_MEHP$`#CHROM`,main="MEHHP_MEHP",signif = 5e-8)
dev.off()

#calculate lambda
pvalue<-MEHHP_MEHP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MEHHP_MEHP.png")
qq(MEHHP_MEHP$P,main="MEHHP_MEHP lambda=0.9876547")
dev.off()

#let's use a a limit P<1e-8

MEHHP_MEHP_bonf<-subset(MEHHP_MEHP, P<1e-8)

head(MEHHP_MEHP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MEHHP_MEHP_04<-subset(MEHHP_MEHP, P<1e-4)

dim(MEHHP_MEHP_04)

#Manhattan plot with the subsetted data

png("Manhattan_MEHHP_MEHP_04.png")
manhattanPlot(MEHHP_MEHP_04$P,MEHHP_MEHP_04$`#CHROM`,main="MEHHP_MEHP limit P<1e-4",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEHHP_MEHP_04$P)

#qq plot with all samples
png("qq_MEHHP_MEHP_04.png")
qq(MEOHP_04$P,main="MEHHP_MEHP limit P<1e-4,lambda=2.151993e-04")
dev.off()

#let's use an other limit P<1e-5

MEHHP_MEHP_05<-subset(MEHHP_MEHP, P<1e-5)

dim(MEHHP_MEHP_05)

#Manhattan plot with the subsetted data

png("Manhattan_MEHHP_MEHP_05.png")
manhattanPlot(MEHHP_MEHP_05$P,MEHHP_MEHP_05$`#CHROM`,main="MEHHP_MEHP limit P<1e-5",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEHHP_MEHP_05$P)

#qq plot with all samples
png("qq_MEHHP_MEHP_05.png")
qq(MEHHP_MEHP_05$P,main="MEHHP_MEHP limit P<1e-5,lambda=1.988264e-05")
dev.off()

#get an excel file with the top SNPs

MEHHP_MEHP_topSNPs<-subset(MEHHP_MEHP, P<1e-5)

write.xlsx(MEHHP_MEHP_topSNPs, file = "MEHHP_MEHP_topSNPs.xlsx")

#################MECPP
MECPP<-fread("MECPP.hs_mecpp_cadj_Log2.glm.linear")

#choose the right test
MECPP<-MECPP[MECPP$TEST=="ADD",]
head(MECPP)

dim(MECPP)

summary(MECPP$P)

#convert to a data.frame
MECPP<-as.data.frame(MECPP)

#Manhattan plot with all data

png("Manhattan_MECPP.png")
manhattanPlot(MECPP$P,MECPP$`#CHROM`,main="MECPP",signif = 5e-8)
dev.off()

#calculate lambda
pvalue<-MECPP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MECPP.png")
qq(MECPP$P,main="MECPP lambda=1.020796")
dev.off()

#let's use a a limit P<1e-8

MECPP_bonf<-subset(MECPP, P<1e-8)

head(MEHHP_MEHP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MECPP_04<-subset(MECPP, P<1e-4)

dim(MECPP_04)

#Manhattan plot with the subsetted data

png("Manhattan_MECPP_04.png")
manhattanPlot(MECPP_04$P,MECPP_04$`#CHROM`,main="MECPP limit P<1e-4",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MECPP_04$P)

#qq plot with all samples
png("qq_MECPP_04.png")
qq(MECPP_04$P,main="MECPP limit P<1e-4,lambda=1.999622e-05")
dev.off()

#let's use an other limit P<1e-5

MECPP_05<-subset(MECPP, P<1e-5)

dim(MECPP_05)

#Manhattan plot with the subsetted data

png("Manhattan_MECPP_05.png")
manhattanPlot(MECPP_05$P,MECPP_05$`#CHROM`,main="MECPP limit P<1e-5",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MECPP_05$P)

#qq plot with all samples
png("qq_MECPP_05.png")
qq(MECPP_05$P,main="MECPP limit P<1e-5,lambda=2.19478e-06")
dev.off()

#get an excel file with the top SNPs

MECPP_topSNPs<-subset(MECPP, P<0.3e-6)

write.xlsx(MECPP_topSNPs, file = "MMECPP_topSNPs.xlsx")


#####################MECP_MEHP

MECPP_MEHP<-fread("MECPP_MEHP.MECPP_MEHP.glm.linear")

#choose the right test
MECPP_MEHP<-MECPP_MEHP[MECPP_MEHP$TEST=="ADD",]
head(MECPP_MEHP)

dim(MECPP_MEHP)

summary(MECPP_MEHP$P)

#convert to a data.frame
MECPP_MEHP<-as.data.frame(MECPP_MEHP)

#Manhattan plot with all data

png("Manhattan_MECPP_MEHP.png")
manhattanPlot(MECPP_MEHP$P,MECPP_MEHP$`#CHROM`,main="MECPP_MEHP",signif = 5e-8)
dev.off()

#calculate lambda
pvalue<-MECPP_MEHP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MECPP_MEHP.png")
qq(MECPP_MEHP$P,main="MECPP_MEHP lambda=0.9881681")
dev.off()

#let's use a a limit P<1e-8

MECPP_MEHP_bonf<-subset(MECPP_MEHP, P<1e-8)

head(MECPP_MEHP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MECPP_MEHP_04<-subset(MECPP_MEHP, P<1e-4)

dim(MECPP_MEHP_04)

#Manhattan plot with the subsetted data

png("Manhattan_MECPP_MEHP_04.png")
manhattanPlot(MECPP_MEHP_04$P,MECPP_MEHP_04$`#CHROM`,main="MECPP_MEHP limit P<1e-4",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MECPP_MEHP_04$P)

#qq plot with all samples
png("qq_MECPP_MEHP_04.png")
qq(MECPP_MEHP_04$P,main="MECPP_MEHP limit P<1e-4,lambda=2.208451e-04")
dev.off()

#let's use an other limit P<1e-5

MECPP_MEHP_05<-subset(MECPP_MEHP, P<1e-5)

dim(MECPP_MEHP_05)

#Manhattan plot with the subsetted data

png("Manhattan_MECPP_MEHP_05.png")
manhattanPlot(MECPP_MEHP_05$P,MECPP_MEHP_05$`#CHROM`,main="MECPP_MEHP limit P<1e-5",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MECPP_MEHP_05$P)

#qq plot with all samples
png("qq_MECPP_MEHP_05.png")
qq(MECPP_MEHP_05$P,main="MECPP_MEHP limit P<1e-5,lambda=1.897041e-05")
dev.off()

#get an excel file with the top SNPs

MECPP_MEHP_topSNPs<-subset(MECPP_MEHP, P<1e-5)

write.xlsx(MECPP_MEHP_topSNPs, file = "MECPP_MEHP_topSNPs.xlsx")


#####################MEOHP_MEHP

MEOHP_MEHP<-fread("MEOHP_MEHP.MEOHP_MEHP.glm.linear")

#choose the right test
MEOHP_MEHP<-MEOHP_MEHP[MEOHP_MEHP$TEST=="ADD",]
head(MEOHP_MEHP)

dim(MEOHP_MEHP)

summary(MEOHP_MEHP$P)

#convert to a data.frame
MEOHP_MEHP<-as.data.frame(MEOHP_MEHP)

#Manhattan plot with all data

png("Manhattan_MEOHP_MEHP.png")
manhattanPlot(MEOHP_MEHP$P,MEOHP_MEHP$`#CHROM`,main="MEOHP_MEHP",signif = 5e-8)
dev.off()

#calculate lambda
pvalue<-MEOHP_MEHP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MEOHP_MEHP.png")
qq(MEOHP_MEHP$P,main="MEOHP_MEHP lambda=0.9859866")
dev.off()

#let's use a a limit P<1e-8

MEOHP_MEHP_bonf<-subset(MEOHP_MEHP, P<1e-8)

head(MEOHP_MEHP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MEOHP_MEHP_04<-subset(MEOHP_MEHP, P<1e-4)

dim(MEOHP_MEHP_04)

#Manhattan plot with the subsetted data

png("Manhattan_MEOHP_MEHP_04.png")
manhattanPlot(MEOHP_MEHP_04$P,MEOHP_MEHP_04$`#CHROM`,main="MEOHP_MEHP limit P<1e-4",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEOHP_MEHP_04$P)

#qq plot with all samples
png("qq_MEOHP_MEHP_04.png")
qq(MEOHP_MEHP_04$P,main="MEOHP_MEHP limit P<1e-4,lambda=2.179256e-04")
dev.off()

#let's use an other limit P<1e-5

MEOHP_MEHP_05<-subset(MEOHP_MEHP, P<1e-5)

dim(MEOHP_MEHP_05)

#Manhattan plot with the subsetted data

png("Manhattan_MEOHP_MEHP_05.png")
manhattanPlot(MEOHP_MEHP_05$P,MEOHP_MEHP_05$`#CHROM`,main="MEOHP_MEHP limit P<1e-5",signif = 5E-8)
dev.off()

#calculate lambda
qq.chisq(MEOHP_MEHP_05$P)

#qq plot with all samples
png("qq_MEOHP_MEHP_05.png")
qq(MEOHP_MEHP_05$P,main="MEOHP_MEHP limit P<1e-5,lambda=1.426311e-05")
dev.off()

#get an excel file with the top SNPs

MEOHP_MEHP_topSNPs<-subset(MEOHP_MEHP, P<1e-5)

write.xlsx(MEOHP_MEHP_topSNPs, file = "MEOHP_MEHP_topSNPs.xlsx")


#MBZP
MBZP<-fread("MBZP.hs_mbzp_cadj_Log2.glm.linear")

#choose the right test
MBZP<-MBZP[MBZP$TEST=="ADD",]
head(MBZP)

dim(MBZP)

summary(MBZP$P)

#convert to a data.frame
MBZP<-as.data.frame(MBZP)

#Manhattan plot with all data

png("Manhattan_MBZP.png")
manhattanPlot(MBZP$P,MBZP$`#CHROM`,main="MBZP",signif = 1e-5)
dev.off()

#calculate lambda
pvalue<-MBZP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MBZP.png")
qq(MBZP$P,main="MBZP lambda=0.9953749")
dev.off()

#let's use a a limit P<1e-8

MBZP_bonf<-subset(MBZP, P<1e-8)

head(MBZP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MBZP_04<-subset(MBZP, P<1e-4)

dim(MBZP_04)


#let's use an other limit P<1e-5

MBZP_05<-subset(MBZP, P<1e-5)

dim(MBZP_05)



#get an excel file with the top SNPs

MBZP_topSNPs<-subset(MBZP, P<1e-5)

dim(MBZP_topSNPs)

write.xlsx(MBZP_topSNPs, file = "MBZP_topSNPs.xlsx")


#MEP
MEP<-fread("MEP.hs_mep_cadj_Log2.glm.linear")

#choose the right test
MEP<-MEP[MEP$TEST=="ADD",]
head(MEP)

dim(MEP)

summary(MEP$P)

#convert to a data.frame
MEP<-as.data.frame(MEP)

#Manhattan plot with all data

png("Manhattan_MEP.png")
manhattanPlot(MEP$P,MEP$`#CHROM`,main="MEP",signif = 1e-5)
dev.off()

#calculate lambda
pvalue<-MEP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MEP.png")
qq(MEP$P,main="MEP lambda=1.01265")
dev.off()

#let's use a a limit P<1e-8

MEP_bonf<-subset(MEP, P<1e-8)

head(MEP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MEP_04<-subset(MEP, P<1e-4)

dim(MEP_04)


#let's use an other limit P<1e-5

MEP_05<-subset(MEP, P<1e-5)

dim(MEP_05)



#get an excel file with the top SNPs

MEP_topSNPs<-subset(MEP, P<1e-5)

dim(MEP_topSNPs)

write.xlsx(MEP_topSNPs, file = "MEP_topSNPs.xlsx")




#MIBP
MIBP<-fread("MIBP.hs_mibp_cadj_Log2.glm.linear")

#choose the right test
MIBP<-MIBP[MIBP$TEST=="ADD",]
head(MIBP)

dim(MIBP)

summary(MIBP$P)

#convert to a data.frame
MIBP<-as.data.frame(MIBP)

#Manhattan plot with all data

png("Manhattan_MIBP.png")
manhattanPlot(MIBP$P,MIBP$`#CHROM`,main="MIBP",signif = 1e-5)
dev.off()

#calculate lambda
pvalue<-MIBP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MIBP.png")
qq(MIBP$P,main="MIBP lambda=1.013187")
dev.off()

#let's use a a limit P<1e-8

MIBP_bonf<-subset(MIBP, P<1e-8)

head(MIBP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MIBP_04<-subset(MIBP, P<1e-4)

dim(MIBP_04)


#let's use an other limit P<1e-5

MIBP_05<-subset(MIBP, P<1e-5)

dim(MIBP_05)



#get an excel file with the top SNPs

MIBP_topSNPs<-subset(MIBP, P<1e-5)

dim(MIBP_topSNPs)

write.xlsx(MIBP_topSNPs, file = "MIBP_topSNPs.xlsx")



#MNBP
MNBP<-fread("MNBP.hs_mnbp_cadj_Log2.glm.linear")

#choose the right test
MNBP<-MNBP[MNBP$TEST=="ADD",]
head(MNBP)

dim(MNBP)

summary(MNBP$P)

#convert to a data.frame
MNBP<-as.data.frame(MNBP)

#Manhattan plot with all data

png("Manhattan_MNBP.png")
manhattanPlot(MNBP$P,MNBP$`#CHROM`,main="MNBP",signif = 1e-5)
dev.off()

#calculate lambda
pvalue<-MNBP$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

#qq plot with all samples
png("qq_MNBP.png")
qq(MNBP$P,main="MNBP lambda=1.006838")
dev.off()

#let's use a a limit P<1e-8

MNBP_bonf<-subset(MNBP, P<1e-8)

head(MNBP_bonf)
#with this limit I didn't get back anything so I have to use an other limit

#1e<4limit

MNBP_04<-subset(MNBP, P<1e-4)

dim(MNBP_04)


#let's use an other limit P<1e-5

MNBP_05<-subset(MNBP, P<1e-5)

dim(MNBP_05)



#get an excel file with the top SNPs

MNBP_topSNPs<-subset(MNBP, P<1e-5)

dim(MNBP_topSNPs)

write.xlsx(MNBP_topSNPs, file = "MNBP_topSNPs.xlsx")







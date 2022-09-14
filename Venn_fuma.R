setwd("/Volumes/ws_helix/HELIX_analyses/GWAS_chemical_exposure/results/phthalates/FUMA_result/FUMA_MEOHP_MEHP")

library("venn")

MECPP<-read.csv("genes.txt",sep="\t")

MECPP_MEHP<-read.csv("genes.txt",sep="\t")

MEHHP<-read.csv("genes.txt",sep="\t")


MEHHP_MEHP<-read.csv("genes.txt",sep="\t")

MEHP<-read.csv("genes.txt",sep="\t")

MEOHP<-read.csv("genes.txt",sep="\t")

MEOHP_MEHP<-read.csv("genes.txt",sep="\t")

input<-list(MECPP=MECPP$symbol,MECPP_MEHP=MECPP_MEHP$symbol,MEHHP=MEHHP$symbol,MEHHP_MEHP=MEHHP_MEHP$symbol,MEHP=MEHP$symbol,MEOHP=MEOHP$symbol,MEOHP_MEHP=MEOHP_MEHP$symbol)

png("Venn_FUMA.png")
venn(input,ilabels = TRUE,zcolor = "style")
dev.off()

#only chemicals

input2<-list(MECPP=MECPP$symbol,MEHHP=MEHHP$symbol,MEHP=MEHP$symbol,MEOHP=MEOHP$symbol)

png("Venn_FUMA_only_chem.png")
venn(input2,ilabels = TRUE,zcolor = "style")
dev.off()


#only ratios

input3<-list(MECPP_MEHP=MECPP_MEHP$symbol,MEHHP_MEHP=MEHHP_MEHP$symbol,MEOHP_MEHP=MEOHP_MEHP$symbol)

png("Venn_FUMA_only_ratio.png")
venn(input3,ilabels = TRUE,zcolor = "style")
dev.off()

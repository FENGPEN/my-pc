setwd("C:\\Users\\suinlee\\Desktop\\泌尿")

# load R package
# load R package
source("CIBERSORT.R")
library(sva)
library(ConsensusClusterPlus)
library(pheatmap)
library(estimate)
library(survival)
library(survminer)
library(corrplot)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Boruta)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggalluvial)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(maftools)

# customized function
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

# set color
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
red    <- "#f25f5c"
purple <- "#531f7a"
grey   <- "#8693ab"
orange <- "#fa7921"
white  <- "#f2d7ee"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
lightblue <- "#B2EBFF"
darkblue  <- "#1d00ff"
cherry    <- "#700353"
lightgrey <- "#dcddde"
nake <- "#F8C364"
gold <- "#ECE700"
cyan <- "#00B3D0"
sun  <- "#E53435"
peach  <- "#E43889"
violet <- "#89439B"
soil   <- "#EC7D21"
lightgreen <- "#54B642"
darkblue   <- "#21498D"
darkgreen  <- "#009047"
brown      <- "#874118"
seagreen   <- "#008B8A"
jco <- c("#2874C5","#EABF00","#5FC1C2","#C6524A","#80A7DE")



ciber.res <- read.table("exp.txt", sep = "\t", row.names = 1, header = T, check.names = F)
head(ciber.res)


ciber.res <- ciber.res
input_data <- cbind.data.frame(ciber.res


indata <- t(scale(input_data))
tmp <- cor(t(indata))
pdf("Figure 1C correlation plot of immune cells.pdf",width = 10,height = 10)
corrplot(tmp, 
         method = "pie",
         type = "lower",
         diag = FALSE,
         addCoef.col = "black",
         tl.col = "black",
         tl.cex = 0.8,
         tl.srt = 90,
         number.cex = 0.8)
invisible(dev.off())

dev.off()
dev.off()
dev.off()

cc <- ConsensusClusterPlus(d = indata, 
                           maxK = 3, 
                           reps = 100, 
                           pItem = 0.8,
                           pFeature = 1, 
                           clusterAlg = "km", 
                           innerLinkage = "ward.D2", 
                           finalLinkage = "ward.D2", 
                           distance = "pearson", 
                           seed = 2020103,
                           title = "ConsensusCluster",
                           plot = "pdf") 
ICIcluster <- cc[[3]]$consensusClass

tmp <- data.frame(expr = as.numeric(combined.expr.combat["PDCD1", names(ICIcluster)]),
                  ICIcluster = factor(ICIcluster))
my_comparisons <- list(c("1","2"),
                       c("2","3"),
                       c("1","3"))
dp <- ggplot(tmp, aes(x=ICIcluster, y=expr, fill=ICIcluster)) +
  scale_fill_manual(values = c(darkred,jco[1],jco[2])) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width = 0.1, fill = "white",outlier.color = NA) + 
  labs(title="",x="ICIcluster", y = "PD-1 expression") + theme_classic() +
  theme(legend.position = "none") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  stat_compare_means(method = "kruskal.test", label.y = 0)
ggsave("Figure 1E violin plot of PDCD1.pdf", width = 4,height = 4)

tmp <- data.frame(expr = as.numeric(combined.expr.combat["CTLA4", names(ICIcluster)]),
                  ICIcluster = factor(ICIcluster)) # cannot match this gene
my_comparisons <- list(c("1","2"),
                       c("2","3"),
                       c("1","3"))
dp <- ggplot(tmp, aes(x=ICIcluster, y=expr, fill=ICIcluster)) +
  scale_fill_manual(values = c(darkred,jco[1],jco[2])) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width = 0.1, fill = "white",outlier.color = NA) + 
  labs(title="",x="ICIcluster", y = "CTLA4 expression") + theme_classic() +
  theme(legend.position = "none") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
  stat_compare_means(method = "kruskal.test", label.y = 0)
ggsave("Figure 1F violin plot of CTLA4.pdf", width = 4,height = 4)

#--------------#

library(MOVICS)
library(MCPcounter)

combined.expr.combat<-read.table("output_combined_expr.txt", sep = "\t", row.names = 1, header = T, check.names = F)
tcga.clin<- read.table("ICIcluster in pooled cohorts.txt", sep = "\t", row.names = 1, header = T, check.names = F)
# create pseudo MOVICS object
tcga.clin$samID <- rownames(tcga.clin)
tcga.clin$clust <- ifelse(tcga.clin$ICIcluster == "A",1,ifelse(tcga.clin$ICIcluster == "B",2,3))
cmoic.brca <- list("clust.res" = tcga.clin,
                   "mo.method" = "Lee")
library(MOVICS)
library(MCPcounter)
# get marker for expression and methylation
runDEA(dea.method = "limma",
       expr       =combined.expr.combat, # normalized expression data
       moic.res   = cmoic.brca,
       prefix     = "TCGA-BRCA")



annCol <- tcga.clin[,c("ICIcluster"), drop = FALSE]


GSET.FILE <-
  system.file("extdata", "新建文本文档.gmt", package = "MOVICS", mustWork = TRUE)


gsva.res <-
  runGSVA(moic.res = cmoic.brca,
          norm.expr = combined.expr.combat,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method = "gsva", # method to calculate single sample enrichment score
          annCol = annCol,
          
          
          fig.path = getwd(),
          fig.name = " immune GENE SETS OF INTEREST HEATMAP",
          height = 10,
          width = 8)

GSET.FILE <-
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)


gsva.res <-
  runGSVA(moic.res = cmoic.brca,
          norm.expr = combined.expr.combat,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method = "gsva", # method to calculate single sample enrichment score
          annCol = annCol,
          
          
          fig.path = getwd(),
          fig.name = " gene sets of interest SETS OF INTEREST HEATMAP",
          height = 10,
          width = 8)


GSET.FILE <-
  system.file("extdata", "代谢.gmt", package = "MOVICS", mustWork = TRUE)


gsva.res <-
  runGSVA(moic.res = cmoic.brca,
          norm.expr = combined.expr.combat,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method = "gsva", # method to calculate single sample enrichment score
          annCol = annCol,
          
          
          fig.path = getwd(),
          fig.name = " 11111immune GENE SETS OF INTEREST HEATMAP",
          height = 30,
          width = 8)
library(data.table)
maf<- fread("data_mutations_extended.txt")
mutect.dataframe <- function(x){
  # delete rows of Silent
  cut_id <- x$Variant_Classification %in% c("Silent")
  x <- x[!cut_id,]
  somatic_sum <- x %>% group_by(Tumor_Sample_Barcode) %>% summarise(TCGA_sum = n())
}

label <- c("Tumor_Sample_Barcode","Hugo_Symbol","NCBI_Build","Chromosome","Start_Position","End_Position","Strand",
           "Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","HGVSp_Short")

maf <- read_tsv("data_mutations_extended.txt", comment = "#")
maf$Tumor_Sample_Barcode <- paste0(maf$Tumor_Sample_Barcode,"A")

comsam <- intersect(rownames(tcga.clin),maf$Tumor_Sample_Barcode)

# 计算TMB
maf <- maf[which(maf$Tumor_Sample_Barcode %in% comsam),label]



tmb.brca <-
  compTMB(moic.res = cmoic.brca,
          maf = maf,
          rmDup = TRUE, # remove duplicated variants per sample
          rmFLAGS = FALSE, # keep FLAGS mutations
          exome.size = 38, # estimated exome size
          test.method = "nonparametric", # statistical testing method
          fig.name = "DISTRIBUTION OF TMB AND TITV")

segment<-fread("TCGA-BLCA.cnv.tsv")
colnames(segment) <- c("sample","chrom","start","end","value")

fga.brca <-
  compFGA(moic.res = cmoic.brca,
          segment = segment,
          iscopynumber = FALSE, # this is a segmented copy number file
          cnathreshold = 0.2, # threshold to determine CNA gain or loss
          test.method = "nonparametric", # statistical testing method
          fig.name = "BARPLOT OF FGA")

drug.brca <-
  compDrugsen(moic.res = cmoic.brca,
              norm.expr = combined.expr,
              drugs = c("Cisplatin", "Paclitaxel"), # a vector of drug name
              tissueType = "urogenital_system", # choose specific tissue type
              test.method = "nonparametric", # statistical testing method
              prefix = "BOXVIOLIN OF ESTIMATED IC50")
head(drug.brca$Cisplatin)


## mutation
tcga.mut <- read.table("geneMutation.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tcga.mut[tcga.mut == "Wild"] <- 0
tcga.mut[tcga.mut == "Mutation"] <- 1
tmp <- dimnames(tcga.mut)
tcga.mut <- sapply(tcga.mut, as.numeric)
dimnames(tcga.mut) <- tmp
tcga.mut <- as.data.frame(tcga.mut)
## filter features
sort(rowSums(tcga.mut),decreasing = T) # check mutation frequency
elite.mut <- getElites(dat       = tcga.mut,
                       method    = "freq", # must set as 'freq'
                       #elite.num = 80, # note: in this scenario elite.num refers to frequency of mutation
                       elite.pct = 0.1) 

mut.marker <- compMut(moic.res   = cmoic.brca,
                      mut.matrix   = tcga.mut, # binary somatic mutation matrix
                      doWord       = F, # generate table in .docx format
                      doPlot       = F) # draw OncoPrint
mut.marker <- mut.marker[which(as.numeric(mut.marker$pvalue) < 0.05),]
write.table(mut.marker, "significant marker for mutation data.txt",sep = "\t",row.names = F,col.names = T,quote = F)




setwd("C:\\Users\\suinlee\\Desktop\\charac\\预后\\新建文件夹\\新建文件夹\\新建文件夹\\新建文件夹")
rt<-read.delim2("新建文本文档 (2).txt",header = T,as.is = T)
bt<-read.delim2("新建文本文档.txt",header = T,as.is = T)
merge<-merge(rt,bt,by="ID")
write.table(merge, "merge.txt", quote=F, row.names=T,col.names = T,sep = "\t")



setwd("C:\\Users\\suinlee\\Desktop\\糖\\8\\1\\差异免疫细胞浸润\\免疫预测")
# load R packages
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(circlize)
library(gridBase)

# set colors
q.col <- "#E1B52E"
g.col <- "#6AA1AB"
c.col <- "#9AAF57"
m.col <- "#81552E"
# Subtype prediction in IMvigor210 cohort #
train.expr <- read.table("output_combined_expr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
train.subt <- read.table("cluster in pooled cohorts.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

comgene<-intersect(colnames(train.expr),rownames(train.subt))
train.expr<-train.expr[,comgene]
train.subt <-train.subt[comgene,]


# include original clinical information as `clust.res` and a string value for `mo.method` to a list
pseudo.moic.res                 <- list("clust.res" = train.subt,
                                        "mo.method" = "COMBINED")

# make pseudo samID
pseudo.moic.res$clust.res$samID <- rownames(pseudo.moic.res$clust.res)

# make pseudo clust using a mapping relationship
pseudo.moic.res$clust.res$clust <- sapply(pseudo.moic.res$clust.res$cluster,
                                          switch,
                                         
                                          "A" = 1,
                                          "B" = 2,
                                          "C" = 3
                                         
                                         
                                         
                                         
) # relabel Mixed as 4
# run DEA with limma
runDEA(dea.method = "limma",
       expr       = train.expr, # normalized expression data
       moic.res   = pseudo.moic.res,
       overwt     = TRUE,
       prefix     = "cluster")

# get marker template
marker.up <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "cluster", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 1, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       norm.expr     = train.expr,
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE) # no heatmap

dev.off()
dev.off()
dev.off()
# run NTP in training cohort by using up-regulated biomarkers

library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(circlize)
library(gridBase)
library(MOVICS)
library(survival)
library(survminer)

# set colors
q.col <- "#E1B52E"
g.col <- "#6AA1AB"
c.col <- "#9AAF57"
m.col <- "#81552E"
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
red    <- "#f25f5c"
purple <- "#531f7a"
grey   <- "#8693ab"
orange <- "#fa7921"
white  <- "#f2d7ee"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
lightblue <- "#B2EBFF"
darkblue  <- "#1d00ff"
cherry    <- "#700353"
lightgrey <- "#dcddde"
nake <- "#F8C364"
gold <- "#ECE700"
cyan <- "#00B3D0"
sun  <- "#E53435"
peach  <- "#E43889"
violet <- "#89439B"
soil   <- "#EC7D21"
lightgreen <- "#54B642"
darkblue   <- "#21498D"
darkgreen  <- "#009047"
brown      <- "#874118"
seagreen   <- "#008B8A"
jco <- c("#2874C5","#EABF00","#5FC1C2","#C6524A","#80A7DE")
train.ntp.pred <- runNTP(expr       = train.expr,
                         templates  = marker.up$templates, # the template has been already prepared in runMarker()
                         scaleFlag  = TRUE, # scale input data (by default)
                         centerFlag = TRUE, # center input data (by default)
                         doPlot     = TRUE, # to generate heatmap
                         fig.name   = "NTP HEATMAP FOR TRAINING OF COMBINED DATA") 

runKappa(subt1     = pseudo.moic.res$clust.res$clust,
         subt2     = train.ntp.pred$clust.res$clust,
         subt1.lab = "TRUE",
         subt2.lab = "NTP",
         fig.name  = "CONSISTENCY HEATMAP FOR TRAINING between TRUE and NTP OF COMBINED DATA")

# load Imvigor data
test.expr <- read.table("cds.tpm.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
test.subt <- read.table("cds.pData.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)



test.ntp.pred <- runNTP(expr       = test.expr,
                        templates  = marker.up$templates, # the template has been already prepared in runMarker()
                        scaleFlag  = TRUE, # scale input data (by default)
                        centerFlag = TRUE, # center input data (by default)
                        doPlot     = TRUE, # to generate heatmap
                        fig.name   = "NTP HEATMAP FOR IMvigor OF cluster") 
tmp <- test.ntp.pred$ntp.res
write.table(tmp,file = "NTP RESULTS OF IMvigor OF cluster.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
test.subt$Subtype <- sapply(tmp$prediction,
                            switch,
                            "CS1" = "A", 
                            "CS2" = "B",
                            "CS3" = "C"
                            
                            
) 
test.subt$Subtype <- factor(test.subt$Subtype,levels = c("A","B","C"))
fitd <- survdiff(Surv(os, censOS) ~ Subtype,
                 data      = test.subt,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(os, censOS)~ Subtype,
               data      = test.subt,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
names(fit$strata) <- gsub("Subtype=","",names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = TRUE,
                risk.table.col    = "strata",
                palette           = c(q.col, g.col, c.col, m.col,blue,yellow),
                data              = test.subt,
                size              = 1,
                break.time.by     = 2,
                legend.title      = "",
                pval              = TRUE,
                surv.median.line  = "hv",
                xlab              = "Time (month)",
                ylab              = "Survival probability",
                risk.table.y.text = FALSE)
p
pdf.options(reset = TRUE, onefile = FALSE)
pdf("Figure 5B km curve of IMvigor predicted subtype.pdf", width = 5.5, height = 6)
print(p)
dev.off()



fisher.test(table(test.subt$binaryResponse,test.subt$Subtype)) # 0.04567
indata <- as.data.frame.array(table(test.subt$binaryResponse,test.subt$Subtype))
indata$A<- indata$A/colSums(indata)[1]
indata$B <- indata$B/colSums(indata)[2]
indata$C <- indata$C/colSums(indata)[3]


pdf("Figure 5C barplot of subtype regarding immunotherapy in Imvigor210 cohort.pdf", width = 6.5,height = 6)
par(bty="o", mgp = c(2,0.5,0), mar = c(4.1,4.1,2.1,4.1),tcl=-.25, font.main=3, las = 1)
a <- barplot(as.matrix(indata),
             col = c(darkblue,darkred),
             border = NA,
             xlab = "Subtype",
             ylab = "Percentage weight")
box(bty = "l")
legend(x = par("usr")[2]-0.15, y = 1.1,
       xpd = T,
       legend = c("SD/PD","CR/PR"), fill = c(darkred,darkblue),bty = "n",
       y.intersp = 0.8,
       x.intersp = 0.2)
text(a[1],indata[1,1]/2,labels = paste0(round(indata[1,1],3)*100,"%"))
text(a[1],indata[1,1] + (indata[2,1])/2,labels = paste0(round(indata[2,1],3)*100,"%"))
text(a[2],indata[1,2]/2,labels = paste0(round(indata[1,2],3)*100,"%"))
text(a[2],indata[1,2] + (indata[2,2])/2,labels = paste0(round(indata[2,2],3)*100,"%"))
text(a[3],indata[1,3]/2,labels = paste0(round(indata[1,3],3)*100,"%"))
text(a[3],indata[1,3] + (indata[2,3])/2,labels = paste0(round(indata[2,3],3)*100,"%"))
text(a[4],indata[1,4]/2,labels = paste0(round(indata[1,4],3)*100,"%"))
text(a[4],indata[1,4] + (indata[2,4])/2,labels = paste0(round(indata[2,4],3)*100,"%"))
invisible(dev.off())
dev.off()
dev.off()
dev.off()







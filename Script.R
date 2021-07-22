###################################################
# SUPPLEMENTARY_FIGURE_1                          #
###################################################

# LIBRARY
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(factoextra)

# DATA
data <- read.csv("Corr_Year.csv", header=TRUE, sep = ";", check.names=FALSE)  
dim (data)  # 173 X 15 

# CORRELATION 
year.cor <- data[, 2:ncol(data)]
year.cor <- as.matrix(year.cor)

corr_matrix <- cor(year.cor, use = "complete.obs", method="spearman")

cor_year.test <- cor.mtest(corr_matrix, conf.level = .95) 

# PLOT
par(mar=c(1,1,1,1)+0.1)
corrplot(corr_matrix,
         is.corr=TRUE,
         type = "upper",
         method = "ellipse",
         col =colorRampPalette(c("#389638","white","#ff6347"))(200),
         tl.col="black", 
         tl.srt=45,
         tl.cex =11/ncol(corr_matrix),
         p.mat = cor_year.test$p, 
         sig.level = c(.001, .01, .05),
         insig = "label_sig",
         pch.col = "gray13", pch.cex = 1.4,
         number.cex=0.7,
         diag=TRUE, tl.pos="lt")


corrplot(corr_matrix, 
         type="lower",
         method="ellipse",
         addCoef.col = "gray13",
         tl.col = "black",
         col =colorRampPalette(c("#389638","white","#ff6347"))(200),
         insig="blank", 
         p.mat = cor_year.test$p,
         cl.pos="n",
         number.cex=10/ncol(corr_matrix), 
         add=T, diag=FALSE, tl.pos="n")



###################################################
# SUPPLEMENTARY_FIGURE_4                          #
###################################################

# LIBRARY
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(factoextra)

# DATA
data_trait <- read.csv("Corr_Trait.csv", header=TRUE, sep = ";", check.names=FALSE)  
dim(data_trait) #173 8
str(data_trait)

colnames(data_trait)[2] <- "a"
colnames(data_trait)[3] <- "DCh"
colnames(data_trait)[4] <- "Hst"
colnames(data_trait)[5] <- "Stl"
colnames(data_trait)[6] <- "SHL"
colnames(data_trait)[7] <- "SHT"
colnames(data_trait)[8] <- "Fs"


group1 <- read.csv("PCA.csv", header=TRUE, sep = ";", check.names=FALSE)  

data_1 <- merge(data_trait, group1, by="ACCESSION")

data_1 <- data_1[order(data_1$DCh),]
row.names(data_1) <- data_1$ACCESSION  
data_1_temp <- na.omit(data_1)


# PCA
pca.1 <-prcomp(data_1_temp[,c(2:8)], center = TRUE,scale. = TRUE)

fviz_pca_biplot(pca.1, label="var", col.var = "black",
                geom.var = c("arrow", "text"),
                labelsize= 4,
                arrowsize= 1,
                pointsize=3,
                palette = c("dodgerblue4", "#389638", "firebrick2"),
                habillage=data_1_temp$Group, 
                addEllipses=TRUE, 
                ellipse.alpha =0.2, 
                ellipse.level=0.95)

###################################################
# FIGURE_8                                        #
###################################################
library(ggpubr)
data <- read.csv("WxBy_SHL_Fs.txt")
library("ggpubr")
ggboxplot(data, x = "geno", y = "SHL", 
          color = "geno", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Di2Di2", "Di2di2", "di2di2"),
          ylab = "SHL", xlab = "geno")
res.aov <- aov(SHL ~ geno, data = data)
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )

ggboxplot(data, x = "geno", y = "Fs", 
          color = "geno", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Di2Di2", "Di2di2", "di2di2"),
          ylab = "Fs", xlab = "geno")
res.aov <- aov(Fs ~ geno, data = data)
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )

###################################################
# FIGURE_6                                        #
###################################################
# BxO_SHL_Fs
library(ggpubr)
data <- read.csv("BxO_SHL_Fs.txt")
ggboxplot(data, x = "haplo", y = "SHL", 
          order = c("H1H2", "H2"),
          ylab = "SHL", xlab = "geno")
with(data, shapiro.test(SHL[haplo == "H1H2"]))
with(data, shapiro.test(SHL[haplo == "H2"]))
res.ftest <- var.test(SHL ~ haplo, data = data)
res.ftest
res <- t.test(SHL ~ haplo, data = data, var.equal = TRUE)
res

ggboxplot(data, x = "haplo", y = "Fs", 
          order = c("H1H2", "H2"),
          ylab = "SHL", xlab = "geno")
with(data, shapiro.test(Fs[haplo == "H1H2"]))
with(data, shapiro.test(Fs[haplo == "H2"]))
res.ftest <- var.test(Fs ~ haplo, data = data)
res.ftest
res <- t.test(Fs ~ haplo, data = data, var.equal = TRUE)
res

#CxELF2_SHL_Fs
library("ggpubr")
data <- read.csv("CxELF2_SHL_Fs.txt")
ggboxplot(data, x = "haplo", y = "SHL", 
          color = "haplo", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("H1", "H1H2", "H2"),
          ylab = "SHL", xlab = "haplo")
res.aov <- aov(SHL ~ haplo, data = data)
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )

ggboxplot(data, x = "haplo", y = "Fs", 
          color = "haplo", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("H1", "H1H2", "H2"),
          ylab = "Fs", xlab = "haplo")
res.aov <- aov(Fs ~ haplo, data = data)
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )

#panel
library("ggpubr")
data <- read.csv("panel_SHL.txt")
ggboxplot(data, x = "haplo", y = "SHL", 
          color = "haplo", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("H1", "H1H2", "H2"),
          ylab = "SHL", xlab = "haplo")
res.aov <- aov(SHL ~ haplo, data = data)
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals)

data <- read.csv("panel_Fs.txt")
ggboxplot(data, x = "haplo", y = "Fs", 
          color = "haplo", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("H1", "H1H2", "H2"),
          ylab = "Fs", xlab = "haplo")
res.aov <- aov(Fs ~ haplo, data = data)
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )
kruskal.test(SHL ~ haplo, data = data)
pairwise.wilcox.test(data$SHL, data$haplo,
                 p.adjust.method = "BH")
kruskal.test(Fs ~ haplo, data = data)
pairwise.wilcox.test(data$Fs, data$haplo,
                 p.adjust.method = "BH")

###################################################
# SUPPL. FIGURE_8                                 #
###################################################
library("ggpubr")
data <- read.csv("panel_MD_Fs.txt")
ggboxplot(data, x = "geno", y = "Fs", 
          order = c("AA", "AG", "GG"),
          ylab = "Fs", xlab = "Peach_AO_0424020")
res.aov <- aov(SHL ~ geno, data = data)
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals)

ggboxplot(data, x = "geno", y = "MD", 
          order = c("AA", "AG", "GG"),
          ylab = "MD", xlab = "Peach_AO_0424020")
# Compute the analysis of variance
res.aov <- aov(MD ~ geno, data = data)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals)

###################################################
# FIGURE_7B                                       #
###################################################
#panel_LocusG_vs_SHL_vs_Fs
library("ggpubr")
setwd("/home/marco/Desktop/shape/Finale_April2021/Revision_July2021/R_notebooks")
data <- read.csv("panel_G_SHL_Fs_1.txt")
#SHL
ggboxplot(data, x = "LocusG", y = "SHL", 
          order = c("P", "N"),
          ylab = "SHL", xlab = "LocusG")
with(data, shapiro.test(SHL[LocusG == "P"]))
with(data, shapiro.test(SHL[LocusG == "N"]))
res.ftest <- var.test(SHL ~ LocusG, data = data)
res.ftest
res <- t.test(SHL ~ LocusG, data = data, var.equal = TRUE)
res
#Fs
ggboxplot(data, x = "LocusG", y = "Fs", 
          order = c("P", "N"),
          ylab = "Fs", xlab = "LocusG")
with(data, shapiro.test(Fs[LocusG == "P"]))
with(data, shapiro.test(Fs[LocusG == "N"]))
res.ftest <- var.test(Fs ~ LocusG, data = data)
res.ftest
res <- t.test(Fs ~ LocusG, data = data, var.equal = TRUE)
res
#SHLvsHaplo
library("ggpubr")
setwd("/home/marco/Desktop/shape/Finale_April2021/Revision_July2021/R_notebooks")
data <- read.csv("panel_G_SHL_Fs_2.txt")
ggboxplot(data, x = "haplo", y = "SHL", color = "LocusG", palette = c("#00AFBB", "#E7B800"))
res.aov2 <- aov(SHL ~ LocusG + haplo, data = data)
summary(res.aov2)
res.aov3 <- aov(SHL ~ LocusG * haplo, data = data)
res.aov3 <- aov(SHL ~ LocusG + haplo + LocusG:haplo, data = data)
summary(res.aov3)
TukeyHSD(res.aov3, which = "haplo")
TukeyHSD(res.aov3, which = "LocusG")
pairwise.t.test(data$SHL, data$haplo, p.adjust.method = "BH")
library(car)
my_anova <- aov(SHL ~ LocusG * haplo, data = data)
Anova(my_anova, type = "III")
#FsvsHaplo
library("ggpubr")
setwd("/home/marco/Desktop/shape/Finale_April2021/Revision_July2021/R_notebooks")
data <- read.csv("panel_G_SHL_Fs_2.txt")
ggboxplot(data, x = "haplo", y = "Fs", color = "LocusG", palette = c("#00AFBB", "#E7B800"))
res.aov2 <- aov(Fs ~ LocusG + haplo, data = data)
summary(res.aov2)
res.aov3 <- aov(Fs ~ LocusG * haplo, data = data)
res.aov3 <- aov(Fs ~ LocusG + haplo + LocusG:haplo, data = data)
summary(res.aov3)
TukeyHSD(res.aov3, which = "haplo")
TukeyHSD(res.aov3, which = "LocusG")
pairwise.t.test(data$Fs, data$haplo, p.adjust.method = "BH")
library(car)
my_anova <- aov(Fs ~ LocusG * haplo, data = data)
Anova(my_anova, type = "III")
#SHL_vs_Haplo_LocusG
library("ggpubr")
setwd("/home/marco/Desktop/shape/Finale_April2021/Revision_July2021/R_notebooks")
data <- read.csv("panel_G_SHL_Fs_3.txt")
# Compute the analysis of variance
res.aov <- aov(SHL ~ haploG, data = data)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals)
kruskal.test(SHL ~ haploG, data = data)
pairwise.wilcox.test(data$SHL, data$haploG, p.adjust.method = "BH")
#Fs_vs_Haplo_LocusG
library("ggpubr")
setwd("/home/marco/Desktop/shape/Finale_April2021/Revision_July2021/R_notebooks")
data <- read.csv("panel_G_SHL_Fs_3.txt")
# Compute the analysis of variance
res.aov <- aov(Fa ~ haploG, data = data)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals)
kruskal.test(Fs ~ haploG, data = data)
pairwise.wilcox.test(data$Fs, data$haploG, p.adjust.method = "BH")


# BiocManager::install("pasilla")
# BiocManager::install("DESeq")
# BiocManager::install("edgeR")
# BiocManager::install("GO.db")
# BiocManager::install("org.Hs.eg.db")

 # Note: using the devel versions of both packages!
library(DESeq) # version 1.9.11
library(edgeR) # version 2.99.8
library(VennDiagram)
library(GO.db)
library(org.Hs.eg.db)
library(ggplot2) 

# Read in data ------------------------------------------------------------

datafile = read.table("input_data/test_gene_expression_matrix.txt", header=T, row.names = 1) # making column 1 row names
head(datafile)
summary(datafile)

## Read in the data making the row names the first column
counttable <- datafile #creating a count table for later use
head(counttable)

## Make metadata data.frame
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("bud", "bud", "bud", "flower", "flower", "flower"))
meta$condition <- relevel(meta$condition, ref="bud")
meta

## Independent filtering?
# keep_cpm <- rowSums(cpm(counttable)>2) >=2
# keep_quantile <- rowSums(counttable)>quantile(rowSums(counttable), probs=.5)
# addmargins(table(keep_cpm, keep_quantile))
# counttable <- counttable[keep_cpm, ]


# edgeR -------------------------------------------------------------------

## Make design matrix
condition <- relevel(factor(meta$condition), ref="bud")
#libType <- factor(meta$libType)
edesign <- model.matrix(~condition)

## Make new DGEList, normalize by library size, and estimate dispersion allowing possible trend with average count size
e <- DGEList(counts=counttable)
e <- calcNormFactors(e)
e <- estimateGLMCommonDisp(e, edesign)
e <- estimateGLMTrendedDisp(e, edesign) 
e <- estimateGLMTagwiseDisp(e, edesign)

## MDS Plot
plotMDS(e, main="edgeR MDS Plot")

## Biological coefficient of variation plot
plotBCV(e, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance")

## Fit the model, testing the coefficient for the treated vs untreated comparison
efit <- glmFit(e, edesign)
efit <- glmLRT(efit, coef="conditionflower")

## Make a table of results
etable <- topTags(efit, n=nrow(e))$table
etable <- etable[order(etable$FDR), ]
head(etable)


# Looking at which factors were significant
addmargins(table(sig_edgeR=etable$FDR<0.05))
summary(decideTests(efit)) # about 4k down regulated; about 4k up regulated @ 0.05FDR
# PLot Log-fold change against log-counts per million with DE genes highlighted
plotMD(efit)
abline(h=c(-1,1), col="blue")

# Want to pull out the names of the genes that were up regulated and down regulated
test <- as.data.frame(decideTests(efit))
test$gene <- as.factor(row.names(test))
summary(test)
up.reg <- test[test$conditionflower > 0,]
down.reg <- test[test$conditionflower < 0,]

summary(up.reg)
summary(down.reg)

# creating a manhattan plot to show which genes would be up regulated and which genes are down regulated in this scenario
rna.dat <- datafile

rna.dat$gene <- as.factor(row.names(rna.dat))
summary(rna.dat)

# separatign bud phase from flower phase so that we can take the mean of each of these
bud.dat <- rna.dat[,c("gene", "Bud_Rep1", "Bud_Rep2", "Bud_Rep3")]
flower.dat <- rna.dat[,c("gene", "Day2_Rep1", "Day2_Rep2", "Day2_Rep3")]

summary(bud.dat)
bud.mean <- data.frame(gene=bud.dat$gene,
                       value = rowMeans(bud.dat[,c("Bud_Rep1", "Bud_Rep2", "Bud_Rep3")]),
                       type = as.factor("bud"))
summary(bud.mean)
bud.mean$log2.value <- ifelse(bud.mean$value == 0, log2(bud.mean$value + 1e-6), log2(bud.mean$value)) # adding a small value to all zero values to avoid infinity error

flower.mean <- data.frame(gene=flower.dat$gene,
                       value = rowMeans(flower.dat[,c("Day2_Rep1", "Day2_Rep2", "Day2_Rep3")]),
                       type = as.factor("flower"))
flower.mean$log2.value <- ifelse(flower.mean$value == 0, log2(flower.mean$value + 1e-6), log2(flower.mean$value)) # adding a small value to all zero values to avoid infinity error
summary(flower.mean)

# stacking together for graphing purposes

rna.mean.stack <- rbind(bud.mean, flower.mean)

# adding color for the up/down regulation of genes
rna.mean.stack$reg <- as.factor(ifelse(rna.mean.stack$gene %in% up.reg$gene, 1,
                             ifelse(rna.mean.stack$gene %in% down.reg$gene,-1,0)) )

summary(rna.mean.stack)

gn.names <- as.factor(matrix(unlist(strsplit(paste(rna.mean.stack$gene), "_")), nrow=nrow(rna.mean.stack), byrow=T)[,2])
rna.mean.stack$pos <- matrix(unlist(strsplit(paste(gn.names), "-")), nrow=nrow(rna.mean.stack), byrow=T)[,1]
rna.mean.stack$pos <- as.numeric(rna.mean.stack$pos)

ggplot(data=rna.mean.stack[rna.mean.stack$type=="flower",]) + # facet_grid(type~.) +
  geom_point(aes(x=pos, y=-log2.value, col=reg)) +
  scale_color_manual(values=c("blue", "grey50", "red"))

ggplot(data=rna.mean.stack) +  facet_grid(type~.) +
  geom_point(aes(x=pos, y=-log2.value, col=reg)) +
  scale_color_manual(values=c("blue", "grey50", "red"))

# Pulling out the band of values at the top that are colored red in the bud
# Meaning they were turned off in the bud but were then turned on in the flower

rna.mean.stack.bud <- rna.mean.stack[rna.mean.stack$type=="bud",]
summary(rna.mean.stack.bud)
rna.mean.stack.bud$gn.sig <- as.factor(ifelse(rna.mean.stack.bud$log2.value < -10 & rna.mean.stack.bud$reg == 1,"Y","N"))
summary(rna.mean.stack.bud)

gn.on <- ggplot(rna.mean.stack.bud[rna.mean.stack.bud$gn.sig=="Y",]) +
  geom_hline(aes(, yintercept=1),lwd=1000, col="grey65") +
  geom_vline(aes(xintercept=pos, y=1), col="purple") +
  labs(x="Position", y="", fill="", title="Turned On", caption = "n = 99") +
  theme(axis.line=element_line(color="black"), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),  
        panel.background=element_blank(), 
        axis.text.x=element_text(angle=0, color="black", size=16, vjust= 0.5, face="bold"), 
        axis.text.y=element_blank(), 
        strip.text=element_text(face="bold", size=22),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="top",
        legend.key.size = unit(0.75, "cm"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22),
        legend.key = element_rect(fill = "white")) + 
  guides(fill = guide_colourbar(barwidth = 15, barheight = 3,title="")) +
  theme(axis.title.y= element_text(size=24, face="bold")) +
  theme(axis.title.x= element_text(size=24, face="bold")) +
  scale_x_continuous(limits = c(0,max(rna.mean.stack.bud$pos)),breaks  = seq(0,max(rna.mean.stack.bud$pos), by=10000))

# The center of the distribution is just a flat out mess, but there are blips of blue around the 20 mark that merit investigation
# Pulling out the band of values at the top that are colored blue in the flower meaning they were turned on in teh bud and are not turned off
rna.mean.stack.flower <- rna.mean.stack[rna.mean.stack$type=="flower",]
summary(rna.mean.stack.flower)
rna.mean.stack.flower$gn.sig <- as.factor(ifelse(rna.mean.stack.flower$log2.value < -10 & rna.mean.stack.flower$reg == -1,"Y","N"))

summary(rna.mean.stack.flower)

gn.off <- ggplot(rna.mean.stack.flower[rna.mean.stack.flower$gn.sig=="Y",]) +
  geom_hline(aes(, yintercept=1),lwd=1000, col="grey65") +
  geom_vline(aes(xintercept=pos, y=1), col="forestgreen") +
  labs(x="Position", y="", fill="", title = "Turned Off", caption ="n = 108") +
  theme(axis.line=element_line(color="black"), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),  
        panel.background=element_blank(), 
        axis.text.x=element_text(angle=0, color="black", size=16, vjust= 0.5, face="bold"), 
        axis.text.y=element_blank(), 
        strip.text=element_text(face="bold", size=22),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="top",
        legend.key.size = unit(0.75, "cm"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22),
        legend.key = element_rect(fill = "white")) + 
  guides(fill = guide_colourbar(barwidth = 15, barheight = 3,title="")) +
  theme(axis.title.y= element_text(size=24, face="bold")) +
  theme(axis.title.x= element_text(size=24, face="bold")) +
  scale_x_continuous(limits = c(0,max(rna.mean.stack.bud$pos)),breaks  = seq(0,max(rna.mean.stack.bud$pos), by=10000))


library(cowplot)

png("figures/gene_onoff.png", height= 8, width = 15, res=300, unit = "in")
  plot_grid(gn.on,gn.off,align = c("v","h"), nrow = 2)
dev.off()



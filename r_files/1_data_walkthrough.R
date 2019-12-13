# Ball Horticulture Data Analysis Test
# Data Walkthrough
library(ggplot2)
# Load in data file

rna.dat <- read.table("input_data/test_gene_expression_matrix.txt", header=T)
summary(rna.dat)
dim(rna.dat)

# Quick scatter plots
ggplot(data=rna.dat) +
  geom_point(aes(x=log2(Bud_Rep1), y=log2(Day2_Rep1)))+
  labs(title="Rep1")

ggplot(data=rna.dat) +
  geom_point(aes(x=log2(Bud_Rep2), y=log2(Day2_Rep2))) +
  labs(title="Rep2")

ggplot(data=rna.dat) +
  geom_point(aes(x=log2(Bud_Rep3), y=log2(Day2_Rep3))) +
  labs(title="Rep3")

# Data looks to be structured where individual genes were sampled for across different buds
# Samples are paired with samples being taken during bud stage and then samples being taken 2 days after flowers open

# Making a new dataframe to make the graphing a bit easier
# Buds
## rep1
bud1 <- data.frame(gene = rna.dat$gene,
                   value = rna.dat$Bud_Rep1)
bud1$stage <- as.factor("bud")  
bud1$rep <- as.factor("rep1")

# Buds
## rep2
bud2 <- data.frame(gene = rna.dat$gene,
                   value = rna.dat$Bud_Rep2)
bud2$stage <- as.factor("bud")  
bud2$rep <- as.factor("rep2")

# Buds
## rep3
bud3 <- data.frame(gene = rna.dat$Name,
                   value = rna.dat$Bud_Rep3)
bud3$stage <- as.factor("bud")  
bud3$rep <- as.factor("rep3")

# Flowers
## rep1
flr1 <- data.frame(gene = rna.dat$Name,
                   value = rna.dat$Day2_Rep1)
flr1$stage <- as.factor("flower")  
flr1$rep <- as.factor("rep1")

## rep2
flr2 <- data.frame(gene = rna.dat$Name,
                   value = rna.dat$Day2_Rep2)
flr2$stage <- as.factor("flower")  
flr2$rep <- as.factor("rep2")

## rep1
flr3 <- data.frame(gene = rna.dat$Name,
                   value = rna.dat$Day2_Rep3)
flr3$stage <- as.factor("flower")  
flr3$rep <- as.factor("rep3")


rna.stack <- rbind(bud1, bud2, bud3, flr1, flr2, flr3)

summary(rna.stack)

gn.names <- as.factor(matrix(unlist(strsplit(paste(rna.stack$gene), "_")), nrow=nrow(rna.stack), byrow=T)[,2])
rna.stack$pos <- matrix(unlist(strsplit(paste(gn.names), "-")), nrow=nrow(rna.stack), byrow=T)[,1]
rna.stack$pos <- as.numeric(rna.stack$pos)

png("figures/manhattan_test.png", height=8, width=15, res=300, unit="in")
  ggplot(data = rna.stack[paste0(rna.stack$gene, rna.stack$rep) %in% paste0(names.express$gene, names.express$rep),]) + facet_grid(rep~stage) +
    geom_point(aes(x=pos, y=log2(value)))
dev.off()


# Creating a stacked dataframe of ratios to look at gene expression
# rep1
rna.exp1 <- data.frame(gene = rna.dat$Name,
                          value = rna.dat$Day2_Rep1/rna.dat$Bud_Rep1)

rna.exp1$rep <- as.factor("rep1") 
rna.exp1$type <- as.factor("ratio")

# getting the names from replicate 1 that did show changes
rep1.name.chng <- rna.exp1[!rna.exp1$value==Inf & !rna.exp1$value==1 & !is.na(rna.exp1$value),"gene"]
summary(rep1.name.chng)
length(rep1.name.chng)
# rep2
rna.exp2 <- data.frame(gene = rna.dat$Name,
                       value = rna.dat$Day2_Rep2/rna.dat$Bud_Rep2)

rna.exp2$rep <- as.factor("rep2") 
rna.exp2$type <- as.factor("ratio")

# getting the names from replicate 2 that did show changes
rep2.name.chng <- rna.exp2[!rna.exp2$value==Inf & !rna.exp2$value==1 & !is.na(rna.exp2$value),"gene"]
summary(rep2.name.chng)
length(rep2.name.chng)

# rep3
rna.exp3 <- data.frame(gene = rna.dat$Name,
                       value = rna.dat$Day2_Rep3/rna.dat$Bud_Rep3)

rna.exp3$rep <- as.factor("rep3") 
rna.exp3$type <- as.factor("ratio")

# getting the names from replicate 2 that did show changes
rep3.name.chng <- rna.exp3[!rna.exp3$value==Inf & !rna.exp3$value==1 & !is.na(rna.exp3$value),"gene"]
summary(rep3.name.chng)
length(rep3.name.chng)


# making data frame of names that actually showed changes between life stages
change.names1 <- data.frame(gene = rep1.name.chng,
                           rep = as.factor("rep1"))

change.names2 <- data.frame(gene = rep2.name.chng,
                            rep = as.factor("rep2"))

change.names3 <- data.frame(gene = rep3.name.chng,
                            rep = as.factor("rep3"))

names.express <- rbind(change.names1, change.names2, change.names3)

# Stacking expression data

rna.expression <- rbind(rna.exp1, rna.exp2, rna.exp3)
summary(rna.expression)

# Splitting the name apart as I think these are teh gene positions
gn.names <- matrix(unlist(strsplit(paste(rna.expression$gene), "_")), nrow=nrow(rna.expression), byrow=T)[,2]
rna.expression$pos <- matrix(unlist(strsplit(paste(gn.names), "-")), nrow=nrow(rna.expression), byrow=T)[,1]
rna.expression$pos <- as.numeric(rna.expression$pos)

png("figures/ratio_manhattan_test.png", height=8, width=15, res=300, unit="in")
  ggplot(data=rna.expression[rna.expression$type=="ratio" & paste0(rna.expression$gene, rna.expression$rep) %in% paste0(names.express$gene, names.express$rep),]) + facet_grid(rep~.) +
    geom_point(aes(x=pos, y=log2(value)))
dev.off()




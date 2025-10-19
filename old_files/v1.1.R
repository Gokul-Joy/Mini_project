install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install("GEOquery")
library(Biobase)
library(GEOquery)
#-----------------
gse<-getGEO(filename="pancreatic.txt")
expression_data<-exprs(gse)
#-----------------
dim(expression_data)
colSums(is.na((expression_data)))

#-----------------
pheno<-pData(gse)
head(pheno$`tissue:ch1`)
unique(pheno$`tissue:ch1`)

#-----------------
labels<-ifelse(pheno$`tissue:ch1`=='Pancreatic tumor',1,0)
table(labels)

#-----------------
summary(expression_data)
hist(expression_data)
hist(expression_data, breaks = 100, main = "Distribution of Expression Values", xlab = "Expression")
range(expression_data)

rm(list = ls())
library(tidyverse)
library(openxlsx)
## get the work path

path_initial<- getwd()


rawdata_path <- paste0(path_initial,"/1.Study raw data")
cleandata_path <- paste0(path_initial,"/2. clean data")


#Intersection of related genes between NSCLC and COVID-19

## NSCLC related genes
setwd(cleandata_path)
load("TCGA-NSCLC_diffgene.Rdata")

setwd(rawdata_path)

NSCLC_genes <- read.xlsx("NSCLC related genes from public databases .xlsx")


## NSCLC related DEGs
TCGA_NSCLC_diffgene %>% 
  filter(abs(log2FoldChange)>1&padj<0.05) %>% 
  pull(gene_id) %>% 
  unique()->NSCLC_degs_all   

## intersection of the NSCLC DEGs from TCGA  and NSCLC 
## related genes from public databases 
intersect(NSCLC_degs_all,
          unique(NSCLC_genes$gene)) %>% 
  unique()->
  NSCLC_genes_inter



## COVID 19 related genes
setwd(cleandata_path)
load("GSE147507_diffgene_new.Rdata")
load("GSE157103_diffgene.Rdata")
load("GSE166190_diffgene.Rdata")


setwd(rawdata_path)

COVID_genes <- read.xlsx("COVID-19 related genes from public databases.xlsx")

## unite the COVID-19 DEGs from GEO 
c(

  GSE147507_diffgene_new %>% 
    filter(abs(log2FoldChange)>1&padj<0.05) %>% 
    pull(gene_id),
  GSE157103_diffgene %>% 
    filter(abs(log2FoldChange)>1&padj<0.05) %>% 
    pull(gene_id),
  GSE166190_diffgene %>% 
    filter(abs(log2FoldChange)>1&padj<0.05) %>% 
    pull(gene_id)
) %>% 
  unique()->COVID_degs_all

## intersection of the COVID-19 DEGs from GEO  and COVID-19 
## related genes from public databases 
intersect(COVID_degs_all,unique(COVID_genes$Gene))->
  COVID_genes_all



intersect(unique(NSCLC_genes_inter),
          unique(COVID_genes_all)) %>% 
  unique()->interact_genes

### count
interact_genes %>% n_distinct()

#[1] 122

## The interact_genes could be applied for PPI network construction
## and KEGG pathway/ GO analysis
### write out the interact_genes

setwd(cleandata_path)
data.frame(interact_genes=interact_genes) %>% 
  write.xlsx("interact_genes.xlsx")
  
###


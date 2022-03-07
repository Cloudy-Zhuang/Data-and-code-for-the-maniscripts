# KEGG and GO enrichment for the genes in the biggest sub-network
# the genes in the biggest sub-network was provided in the 2.clean data file


library(clusterProfiler)
library(tidyverse)

## get the work path
path_initial<- getwd()
cleandata_path <- paste0(path_initial,"/2. clean data")

## set work place

setwd(cleandata_path)


read.csv("genes of the biggest sub-network.csv") %>% 
  pull(name) %>% 
  unique()->gene

print(gene)

mygene <- gene

gene = bitr(mygene, fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")



## KEGG
EGG <- enrichKEGG(gene= gene$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

EGG <- setReadable(EGG ,OrgDb="org.Hs.eg.db",keyType ="ENTREZID" )
cluster_EGG <- setReadable(EGG ,OrgDb="org.Hs.eg.db",keyType ="ENTREZID" )


## GO
go <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
cluster_go<- setReadable(go ,OrgDb="org.Hs.eg.db",keyType ="ENTREZID" )


## save data

save(cluster_EGG,cluster_go,
     file="sub-network_KEGG_GO.Rdata")

## KEGG GO enrichment
library(clusterProfiler)
library(tidyverse)
library(openxlsx)

path_initial<- getwd()
cleandata_path <- paste0(path_initial,"/2. clean data")


# read data and analyze---------------------------------------------------------------
setwd(cleandata_path)
mygene <- read.xlsx("interact_genes.xlsx")

mygene <- mygene$interact_genes 

gene = bitr(mygene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")



## KEGG enrichment
EGG <- enrichKEGG(gene= gene$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

#set the result readable
EGG <- setReadable(EGG ,OrgDb="org.Hs.eg.db",keyType ="ENTREZID" )
diseases_EGG <- setReadable(EGG ,OrgDb="org.Hs.eg.db",keyType ="ENTREZID" )

## GO enrichment
go <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
diseases_GO  <- setReadable(go ,OrgDb="org.Hs.eg.db",keyType ="ENTREZID" )


go_res <- diseases_GO@result




## save resluts

save(diseases_EGG,
      diseases_GO,
     file = "GO_KEGG_result.Rdata")



rm(list = ls())

## get the work path

path_initial<- getwd()

GEOdataset_path <- paste0(path_initial,"/1.Study raw data/COVID-19 GEO dataset")

## load R package
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(data.table)


#Computing DEGs for every dataset
# # GSE147507 -------------------------------------------------------------
setwd(GEOdataset_path)

load("GSE147507.Rdata")

# expression data
GSE147507_count <-fread("GSE147507_RawReadCounts_Human.tsv")
GSE147507_count %>% 
  column_to_rownames("gene")->GSE147507_count



## create  phenotypic data
GSE147507_pdat %>% 
  filter(!str_detect(title,"IAV|RSV|HPIV3|IFNB|Rux")) %>%   # remove sample not related to COVID19
  dplyr::select(1,2) %>% 
  mutate(group=str_sub(title,1,-3)) %>% 
  mutate(across(group,~str_remove_all(.x,"Series[0-9]{1,2}_"))) %>% 
  arrange(group)->
  GSE147507_pheno

## create raw data based on  phenotypic data
GSE147507_count %>% 
  dplyr::select(GSE147507_pheno$title)->GSE147507_include_exp

# differential analysis by Deseq2

## group index 
GSE147507_pheno$group %>% 
  unique()->diffind 

## for loop for preforming differential analysis in diffrent groups


mylist <- list()
for (i in seq(1,10,by=2)) {
  
  ## creat phenotypic data for specific groups
  GSE147507_pheno %>% 
    filter(group%in%c(diffind[i],diffind[i+1])) %>% 
    select(-2)-> metadata 
  
  ## creat expression data for specific groups
  GSE147507_count %>% 
    select(metadata$title)->exprdata
 
 # standard procedure  of Deseq2
  dds <-DESeqDataSetFromMatrix(countData=exprdata, 
                               colData=metadata, 
                               design=~group,
                               tidy=F)
 
  dds <- dds[rowSums(counts(dds))>1,]
  dds <- DESeq(dds,parallel = T)

  contrast <- c("group",diffind[i+1],diffind[i])
  
  res <- results(dds, contrast=contrast, alpha = 0.05)%>% 
    data.frame() %>% 
    rownames_to_column("gene_id") %>% 
    mutate(group=paste0(diffind[i+1],"vs.",diffind[i]))
  
  mylist[[(i+1)/2]] <- res
  
}

## combine output data
GSE147507_diffgene_new <- as.data.frame(do.call(rbind,mylist))

## save data

cleandata_path <- paste0(path_initial,"/2. clean data")

setwd(cleandata_path)

save(GSE147507_diffgene_new ,
     file = "GSE147507_diffgene_new.Rdata")



# # GSE166190 -------------------------------------------------------------
setwd(GEOdataset_path)
## expression data
GSE166190_count <- fread("GSE166190_Raw_counts.txt.gz") %>% 
  column_to_rownames("gene_id")



# phenotypic data
GSE166190_pheno <-colnames(GSE166190_count)

# create clean index of phenotypic data
GSE166190_pheno  %>% 
  str_remove("Cluster[1-5]{1}_") %>%
  str_remove("_Interval[1-5]{1}|_Visit[1-5]{1}") %>% 
  str_remove_all("[0-9]{1,2}") %>% 
  trimws() %>% 
  cbind(GSE166190_pheno) %>% 
  as.data.frame() %>% 
  dplyr::rename("group"=1) %>% 
  arrange(group)->GSE166190_pheno_clean

# unique comparison index
diffind <- unique(GSE166190_pheno_clean$group)





## for loop for preforming differential analysis in diffrent groups
mylist <- list()
for (i in seq(1,4,by=2)) {
  
  GSE166190_pheno_clean %>% 
    filter(group%in%c(diffind[i],diffind[i+1]))-> 
    metadata 
  
  GSE166190_count %>% 
    dplyr::select(metadata$GSE166190_pheno)->exprdata

  
  # standard procedure  of Deseq2
  dds <-DESeqDataSetFromMatrix(countData=exprdata, 
                               colData=metadata, 
                               design=~group,
                               tidy=F)
  dds <- dds[rowSums(counts(dds))>1,]
  dds <- DESeq(dds,parallel = T)
  contrast <- c("group",diffind[i+1],diffind[i])
  
  res <- results(dds, contrast=contrast, alpha = 0.05)%>% 
    data.frame() %>% 
    rownames_to_column("gene_id") %>% 
    mutate(group=paste0(diffind[i+1],"vs.",diffind[i]))
  
  mylist[[(i+1)/2]] <- res
  
}

## combine output data
GSE166190_diffgene <- as.data.frame(do.call(rbind,mylist))


###GSE166190 result needs gene name annotation 
if (T) {
  
  library(org.Hs.eg.db)
  

  gene_symbol <-clusterProfiler::bitr(GSE166190_diffgene$gene_id, 
                                      fromType="ENSEMBL",
                                      toType=c("SYMBOL", "ENTREZID"),
                                      OrgDb="org.Hs.eg.db")
  ## Unification of the column
  names(GSE166190_diffgene)[1] <- names(gene_symbol)[1]
  names(gene_symbol)[2] <- "gene_id"
  
  GSE166190_diffgene %>% 
    inner_join(dplyr::select(gene_symbol,1,2)) %>% 
    dplyr::select(-ENSEMBL)->GSE166190_diffgene
}


## save output data
cleandata_path <- paste0(path_initial,"/2. clean data")

setwd(cleandata_path)

save(GSE166190_diffgene,file = "GSE166190_diffgene.Rdata")



# # GSE157103  ------------------------------------------------------------
### expression data
setwd(GEOdataset_path)
GSE157103_count <- fread("GSE157103_genes.ec.tsv.gz")

### phenotypic data
GSE157103_med <- openxlsx::read.xlsx("GSE157103 Metadata.xlsx")

### create a regular expression  to match the eligible data
mypat <- "(COVID|NONCOVID)_([0-9]{1,3})_[0-9]{1,3}y_(male|female|unknown)_(ICU|NonICU)"

### create phenotypic data
str_match(GSE157103_med$ID,pattern = mypat) %>% 
  as.data.frame() %>% 
  `colnames<-`(c("ID","dis","num","gender","hos"))->
  metadat

## Adjust the order of phenotypic data based on expression matrix column names
metadat$num[1:9] <-as.character(1:9)
metadat$num[1:100] <- paste0("C",metadat$num[1:100]) 
metadat$num[101:109] <- as.character(1:9)
metadat$num[101:126] <- paste0("NC",metadat$num[101:126]) 

## Delete samples whose gender is unknown
metadat_clean <- filter(metadat,gender!="unknown") %>% 
  mutate(group=paste(gender,dis,hos,sep = "-"))

### rename the column name of the expression matrix
GSE157103_include_exp <- GSE157103_count
names(GSE157103_include_exp)[-1] <- paste0(metadat$dis,"-",metadat$num,"-",metadat$gender,"-",metadat$hos)


## for loop for preforming differential analysis in diffrent groups
## check the group index
metadat_clean$group %>% sort %>% unique()
## create the contrast group neme
versgroup <- list(list("female-COVID-ICU","female-NONCOVID-ICU"),
                  list("female-COVID-NonICU","female-NONCOVID-NonICU"),
                  list("male-COVID-ICU","male-NONCOVID-ICU"),
                  list("male-COVID-NonICU","male-NONCOVID-NonICU"))

mylist <- list()

for (i in seq_along(versgroup)) {

  GSE157103_count %>% 
    rename("symbol"=1) %>% 
    column_to_rownames("symbol") %>% 
    select(filter(metadat_clean,group%in%unlist(versgroup[[i]])) %>% 
             pull(num)) %>% 
    mutate_all(as.integer)->exprdata

  metadat_clean %>%
    filter(group%in%unlist(versgroup[[i]])) %>% 
    select(num,group)->med
  # standard procedure  of Deseq2
  dds <-DESeqDataSetFromMatrix(countData=exprdata, 
                               colData=med, 
                               design=~group,
                               tidy=F)

  dds <- dds[rowSums(counts(dds))>1,]
  dds <- DESeq(dds,parallel = T)

  contrast <- c("group",unlist(versgroup[[i]])[1],unlist(versgroup[[i]])[2])
  
  res <- results(dds, contrast=contrast, alpha = 0.05)%>% 
    data.frame() %>% 
    rownames_to_column("gene_id") %>% 
    mutate(group=paste0(unlist(versgroup[[i]])[1],"vs.",unlist(versgroup[[i]])[2]))
  
  mylist[[i]] <- res
  
}

cleandata_path <- paste0(path_initial,"/2. clean data")
setwd(cleandata_path)

GSE157103_diffgene <- as.data.frame(do.call(rbind,mylist))
save(GSE157103_diffgene,file = "GSE157103_diffgene.Rdata")

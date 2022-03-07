rm(list = ls())
library(tidyverse)
library(data.table)
library(DESeq2)

## get the work path

path_initial<- getwd()
TCGAdata_path <- paste0(path_initial,"/1.Study raw data/TCGA-NSCLC data")


## Counts matrix data combination

### Read data 
LUAD_count <- fread("TCGA-LUAD.htseq_counts.tsv.gz",data.table = F)
LUSC_count <- fread("TCGA-LUSC.htseq_counts.tsv.gz",data.table = F)
annotation_dat <- fread("gencode.v22.annotation.gene.probeMap")



## create annotation function for gene annotation
Annotate <- function(dat){
  dat %>% 
    dplyr::rename("id"=1) %>% 
    inner_join(select(annotation_dat,id,gene)) %>% 
    select(-id) %>% 
    mutate(ind=rowSums(across(where(is.numeric)))) %>% 
    arrange(desc(ind)) %>% 
    distinct(gene,.keep_all = T) %>%  #Delete duplicate genes
    select(-ind) %>% 
    arrange(gene) %>% 
    column_to_rownames("gene") %>%  
    mutate_all(~.x^2-1)->out_dat   # convert back from log2(count+1)
  out_dat  
  
}

LUAD_countanno <- Annotate(LUAD_count)
LUSC_countanno <- Annotate(LUSC_count)

## create function for phenotypic data

Creatmd <- function(dat){
  ind <- as.numeric(substr(colnames(dat),14,15))
  group <- ifelse(ind%in%seq(1,9),"cancer","normal")
  metadat <- data.frame(id=colnames(dat),group)
  metadat$group <- factor(group,levels = c("normal","cancer"))
  metadat 
}

LUAD_countanno %>% 
  Creatmd() %>% 
  mutate(type="LUAD")->LUAD_md
LUAD_md$group%>% table()


LUAD_md$id %>% substr(1,12) %>% n_distinct()


LUSC_countanno %>% 
  Creatmd() %>% 
  mutate(type="LUSC")->LUSC_md
LUSC_md$group %>% table()
LUSC_md$id %>% substr(1,12) %>% n_distinct()

#rbind metadata/counts

NSCLC_med <- rbind(LUAD_md,LUSC_md)

NSCLC_med$group %>% table() 
NSCLC_counts <- cbind(LUAD_countanno,LUSC_countanno)

## aviod the invalid expression value
NSCLC_counts[NSCLC_counts<0] <- 0 

# keep the genes that at least have one expression value not equal to 0
keepGene=rowSums(edgeR::cpm(NSCLC_counts)>0) >=2

table(keepGene)
mycounts <-NSCLC_counts[keepGene,]
mycounts <- mutate_all(mycounts,as.integer) 

# ## standard Deseq2 process -------------------------------------------------------

dds <-DESeqDataSetFromMatrix(countData=mycounts, 
                             colData=NSCLC_med, 
                             design=~group,
                             tidy=F)
### gene counts less than or equal to 1 in all samples, then delete it  
dds <- dds[rowSums(counts(dds))>1,]
dds <- DESeq(dds,parallel = T)

# Extract the difference analysis results

contrast <- c("group","cancer","normal")

results(dds, contrast=contrast, alpha = 0.05)%>% 
  data.frame() %>% 
  rownames_to_column("gene_id") %>% 
  mutate(group=paste0("Cancer","vs.","normal"))->
  TCGA_NSCLC_diffgene

## save data
cleandata_path <- paste0(path_initial,"/2. clean data")

setwd(cleandata_path)

save(TCGA_NSCLC_diffgene,file="TCGA-NSCLC_diffgene.Rdata")

setwd(path_initial)

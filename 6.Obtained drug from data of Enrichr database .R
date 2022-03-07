## identify the drug related to the interaction hub 
## genes between NSCLC and COV19-19
library(tidyverse)
library(data.table)
library(openxlsx)

## get the work path
path_initial<- getwd()
cleandata_path <- paste0(path_initial,"/2. clean data")

## set work place

setwd(cleandata_path)

### read the data of interactional hub genes
hubgene <- read.xlsx("interactional hub genes.xlsx")
names(hubgene) <- "name"

##Create a function to extract the drug corresponding 
## to the core interacting gene

if (T) {
  Creatdrugtab <- function(dat){
    dat %>% 
      rowwise() %>% 
      mutate(genecount=sum(!is.na(across(2:ncol(.))))) %>%  #Extract data with no missing values
      select(V1,genecount,everything()) %>% 
      pivot_longer(cols = -c(1:2),
                   names_to = "name",
                   values_to = "gene") %>% 
      drop_na() %>% 
      mutate(drug=str_split_fixed(V1,"human",2)[,1],
             source=str_split_fixed(V1,"human",2)[,2]) %>% 
      select(-V1,-name) %>% 
      select(drug,source,gene,everything()) %>% 
      filter(gene%in%hubgene$name)->ts1
    
    ts1 %>% 
      group_by(drug) %>% 
      summarise(hubgene_count=length(unique(gene))) %>% 
      inner_join(ts1)->ts2
    ts2
  }
}

##Extract drugs, and identify number of the corresponding genes 
## and the number of interactional hub genes



#Drug_Perturbations_up data and the corresponding genes  --------------


Drug_up <- Creatdrugtab(Drug_Perturbations_up_hsa)

## show the result in the table
Drug_up %>% 
  filter(hubgene_count>=11) %>% #Extract drugs corresponding to more than 11 interactional hub genes
  mutate(Regulation="Up-regulated") %>% 
  arrange(desc(hubgene_count)) %>% 
  distinct(drug,.keep_all = T) %>% 
  gt::gt()

### Present the result in another form
Drug_up %>% 
  filter(hubgene_count>=11) %>% 
  mutate(Regulation="Up-regulated") %>% 
  arrange(desc(hubgene_count))->dat

drug <- unique(dat$drug)

Drug_updf <- data.frame()

## need a for loop
for (i in seq_along(drug)) {
  
  Drug_updf[i,1] <- drug[i]
  Drug_updf[i,2] <- filter(dat,drug==drug[i]) %>% 
    pull(hubgene_count) %>% 
    unique()
  
  Drug_updf[i,3] <- filter(dat,drug==drug[i]) %>% 
    pull(genecount) %>% 
    .[1]
  
  Drug_updf[i,4] <- filter(dat,drug==drug[i]) %>% 
    pull(Regulation) %>% 
    unique()
  
  dat %>% 
    filter(drug==drug[i]) %>% 
    pull(gene) %>% 
    unique() %>% 
    sort() %>% 
    paste(collapse = "/")->Drug_updf[i,5]
  
  
}


## according to the reference from Pubmed 
## estradio and tamoxifen can not treat NSCLC or COVID-19
## thus these two drugs were removed from the result

Drug_updf %>% 
  rename("Drug"=1,
         "Hub gene count"=2,
         "Total gene count"=3,
         "Regulation"=4,
         "Hub gene list"=5) %>% 
  mutate(Data_source="Drug_Perturbations_from_GEO") %>% 
  filter(!str_detect(Drug,"estradiol|tamoxifen"))->
  Drug_updf_clean



# Drug_Perturbations_down data and the corresponding genes -------------------------------------------------------------------


Drug_down <- Creatdrugtab(Drug_Perturbations_down_hsa)

## show the result in the table
Drug_down %>% 
  mutate(Regulation="Down-regulated") %>% 
  arrange(desc(hubgene_count)) %>% 
  distinct(drug,.keep_all = T) %>% 
  filter(hubgene_count>=11) %>% 
  gt::gt(rownames_to_stub  = 1)



### Present the result in another form
Drug_down %>% 
  mutate(Regulation="Down-regulated") %>% 
  arrange(desc(hubgene_count)) %>% 
  filter(hubgene_count>=11) ->dat


drug <- unique(dat$drug)

Drug_downdf <- data.frame()

## also need a for loop
for (i in seq_along(drug)) {
  
  Drug_downdf[i,1] <- drug[i]
  Drug_downdf[i,2] <- filter(dat,drug==drug[i]) %>% 
    pull(hubgene_count) %>% 
    unique()
  Drug_downdf[i,3] <- filter(dat,drug==drug[i]) %>% 
    pull(genecount) %>% 
    .[1]
  Drug_downdf[i,4] <- filter(dat,drug==drug[i]) %>% 
    pull(Regulation) %>% 
    .[1]
  
  dat %>% 
    filter(drug==drug[i]) %>% 
    pull(gene) %>% 
    unique() %>% 
    sort() %>% 
    paste(collapse = "/")->Drug_downdf[i,5]
  
  
}



### check the drug result
Drug_downdf$V1 %>% 
  print()

## we carefully search these drugs one by one in the pubmed database
## ## according to the reference from Pubmed 
## only the following drugs are eligable  
# "ascorbic acid","carboplatin",
# "resveratrol",  "vemurafenib",
# "PD173074","azacitidine",
# "decitabine","Y15",
# "neocarzinostatin",
# "phorbol 12-myristate 13-acetate"
## thus these drugs were included


c("ascorbic acid","carboplatin",
  "resveratrol",  "vemurafenib",
  "PD173074","azacitidine",
  "decitabine","Y15",
  "neocarzinostatin",
  "phorbol 12-myristate 13-acetate") %>% 
  sort() %>% 
  paste(collapse = "|")->drug_ind1




Drug_downdf  %>% 
  rename("Drug"=1,
         "Hub gene count"=2,
         "Total gene count"=3,
         "Regulation"=4,
         "Hub gene list"=5) %>% 
  mutate(Data_source="Drug_Perturbations_from_GEO") %>% 
  filter(str_detect(Drug,drug_ind1))->
  Drug_down_df_clean





# DSigDB data -------------------------------------------------------------

#Because DSigDB data is too large, this step requires 
#higher computer performance. We provide the analyzed 
#cleaning data, please read and use it directly.  
if (F) {
  DSigDB %>% 
    select(-V2) %>% 
    rowwise() %>% 
    mutate(gene_count=sum(!is.na(across(2:ncol(.))))) %>% 
    select(V1,gene_count,everything()) %>% 
    pivot_longer(cols =-c(1:2),
                 values_to = "gene",
                 names_to = "name") %>% 
    select(-name)->DSigDB_tab1
}
## load the prepared data
load("DSigDB_tab1.Rda")

## Extract the clean data like the previous step
DSigDB_tab1 %>% 
  drop_na() %>% 
  filter(gene%in%hubgene$name) %>% 
  group_by(V1) %>% 
  summarise(hubgene_count=length(unique(gene))) %>% 
  inner_join(DSigDB_tab1)->DSigDB_tab_clean

## Present the result in the table 
DSigDB_tab_clean %>% 
  filter(hubgene_count>=11) %>% 
  arrange(desc(hubgene_count)) %>% 
  distinct(V1,.keep_all = T) %>%
  arrange(desc(hubgene_count)) %>% 
  gt::gt(rownames_to_stub = 1)

## then we also carefully serached the drug in the table one by one
## in the pubmed database.
##ultimately:
## the following drugs are included
# "cyclosporin A", "LUCANTHONE", "troglitazone",
# "COUMESTROL","resveratrol" , "irinotecan",
# "Dasatinib","quercetin","genistein",
# "menadione","trichostatin A",
# "5-Fluorouracil", "etoposide", 
# "resveratrol"

c("cyclosporin A", "LUCANTHONE", "troglitazone",
  "COUMESTROL","resveratrol" , "irinotecan",
  "Dasatinib","quercetin","genistein",
  "menadione","trichostatin A",
  "5-Fluorouracil", "etoposide", 
  "resveratrol")%>% paste(collapse  = "|")->drug_ind2

## extract the eligable drug data
DSigDB_tab_clean %>% 
  filter(hubgene_count>=11) %>% 
  arrange(desc(hubgene_count)) %>% 
  rename("Drug"=1)->dat

## the for loop to present the data in another form
drug <- unique(dat$Drug)

DSigDB_df <- data.frame()

for (i in seq_along(drug)) {
  
  
  DSigDB_df[i,1] <- drug[i]
  DSigDB_df[i,2] <- filter(dat,Drug==drug[i]) %>% 
    pull(hubgene_count) %>% 
    .[1]
  DSigDB_df[i,3] <-filter(dat,Drug==drug[i]) %>% 
    pull(gene_count) %>% 
    .[1]
  
  
  dat %>% 
    filter(Drug==drug[i]) %>% 
    pull(gene) %>% 
    unique() %>% 
    sort() %>% 
    paste(collapse = "/")->DSigDB_df[i,4]
  
  
}

DSigDB_df %>% 
  rename("Drug"=1,
         "Hub gene count"=2,
         "Total gene count"=3,
         "Hub gene list"=4) %>% 
  filter(str_detect(Drug,drug_ind2)) %>% 
  mutate(Drug=Hmisc::capitalize(Drug)) %>%
  mutate(Data_souce="DSigDB") ->DSigDB_fin



# ##write out the result --------------------------------------------------

## combine the data from same source

Drug_ud_fin <- rbind(Drug_updf_clean,Drug_down_df_clean)




library(openxlsx)


OUT <- createWorkbook()

for (i in c(
  "DSigDB_fin",
  "Drug_ud_fin")) {
  
  addWorksheet(OUT,i)
  writeData(OUT,sheet = i,x = get(i))  
  
}
saveWorkbook(OUT,"all_drug_data_from_Enrichr.xlsx",overwrite=T)

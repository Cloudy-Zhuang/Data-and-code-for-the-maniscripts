library(tidyverse)
library(data.table)

## get the work path
path_initial<- getwd()
Enrichrdata_path <- paste0(path_initial,"/1.Study raw data/Enrichr database raw data")
cleandata_path <- paste0(path_initial,"/2. clean data")


## set work path
setwd(Enrichrdata_path)

## Raw data is not tidy, build cleaning function
if (T) {
  Creatdata <- function(mymatrix) {
    
    
    mymatrix[[1]] %>% 
      strsplit("\t")->my_list
    
    
    matrix_df <- do.call(rbind,
                         lapply(my_list,`length<-`,
                                max(sapply(my_list,length))))
    matrix_df <-matrix_df %>% 
      as.data.frame( )
    matrix_df 
  }
}

# Clean the data based on Drug_Perturbations -------------------------------------------------------------------------


Drug_Perturbations_down <- readLines("Drug_Perturbations_from_GEO_down.txt") %>% 
  as.data.frame()%>% 
  Creatdata() #applied the function

Drug_Perturbations_down %>% 
  select(-2) %>% 
  filter(str_detect(V1,"human|hsa"))->Drug_Perturbations_down_hsa



Drug_Perturbations_up <- readLines("Drug_Perturbations_from_GEO_up.txt") %>% 
  as.data.frame()%>% 
  Creatdata()

###
Drug_Perturbations_up %>% 
  select(-2) %>% # Delete redundant column
  filter(str_detect(V1,"human|hsa"))-> #Extract  data with species as human
  Drug_Perturbations_up_hsa



#Clean the data based on  DSigDB ------------------------------------------------------------------

DSigDB <- readLines("DSigDB.txt") %>% 
  as.data.frame() %>% 
  Creatdata()



## save output data
setwd(cleandata_path)

save(Drug_Perturbations_up,
     Drug_Perturbations_down,
     DSigDB,
     file = "Enrichr database_clean_data.Rdata")




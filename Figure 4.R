# Figure 4
# Figure 4a is the network created by inputing  122 interactional 
# genes between NSCLC and COVID-19 into STRING database(https://string-db.org/)
# The input data was in the raw data  file;
# it named: interact_genes.xlsx
#to assure the highest confidence of the network, the minimum confidence 
#score was set at 0.90 for network construction

# The ouput data from STRING database was then analyzed by Cytoscape 3.9
#(https://cytoscape.org/download.html) to create the data for Figure 4b
# The prepared data was provided in the file of clean data

# Then the code for plotting Figure 4b are as follows
library(tidyverse)



path_initial<- getwd()
cleandata_path <- paste0(path_initial,"/2. clean data")

## set work place

setwd(cleandata_path)

## open a canvas
x11()

read.csv("data for Figure 4b.csv")%>% 
  filter(Degree>2*median(Degree)) %>%  #filter the genes with 2 fold median Degree in the network
  ggplot()+
  geom_point(aes(x=reorder(name,Degree),
                 y=Degree,
                 color=BetweennessCentrality,
                 size=ClosenessCentrality))+
  scale_color_viridis_b()+
  theme_bw(base_size =20)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_size(range = c(5,8))+
  labs(x="Gene name",y="Degree")+
  guides(color=guide_colorbar(order=1))





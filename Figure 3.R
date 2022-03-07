### draw figure 3.
library(clusterProfiler)
library(tidyverse)
library(openxlsx)

## set work path and load data
path_initial<- getwd()
cleandata_path <- paste0(path_initial,"/2. clean data")

setwd(cleandata_path)

load("GO_KEGG_result.Rdata")



### Figure 3 B

diseases_EGG@result  %>% 
  filter(pvalue<0.05) %>% 
  arrange(desc(Count)) %>% 
  mutate(richFactor =Count / as.numeric(sub("/\\d+", "", BgRatio))) %>% 
  arrange(pvalue,desc(richFactor)) %>% 
  ggplot(aes(Description,sort(richFactor),
             color=pvalue,
             size=Count))+
  geom_point()+
  coord_flip()+
  labs(
    #设定图片标题参数 
    x="KEGG pathway",
    y="Rich Factor" ,
    size="Gene Count",
    color="p value")+
  theme_bw()+ 
  theme( 
    axis.text.y = element_text(size = rel(1.8)), 
    axis.text.x = element_text(size = rel(1.8)),
    axis.title.x = element_text(size=rel(1.5)),
    axis.title.y = element_text(size = rel(1.5)),
    title = element_text(size = rel(1.5)),
    plot.title = element_text(hjust = 0.5),
    strip.text.y = element_text(size = rel(1.8)))+
  scale_size(range=c(3,10))+
  theme(legend.position = "right")+
  scale_color_viridis_c(begin = 0.6, end = 1,option="magma")+
  guides(color=guide_colorbar(order=2))->goplot
x11()
print(goplot)





# plot Figure 3 A--------------------------------------------------------------------

x = diseases_GO
## 计算富集分数
x@result$richFactor =x@result$Count / as.numeric(sub("/\\d+", "", x@result$BgRatio))
y =x@result
showCategory = 20


## By referring to the relevant published studies, cell aging
## is the significant pathways shared between NSCLC and COVID-19
## Thus, the GO terms related to cell aging were detected and presented 
## in our study.

y$Description[grepl("aging",y$Description)]

x11()
y %>% 
  filter(pvalue<0.05) %>% 
  group_by(ONTOLOGY) %>% 
  dplyr::slice(1:10) %>%
  bind_rows(filter(y,Description=="regulation of cell aging")) %>% 
  bind_rows(filter(y,Description=="cell aging")) %>% 
  filter(Description!="chromosome segregation") %>% # reduce  duplicate terms 
  filter(Description!="mitotic sister chromatid segregation") %>%  # reduce  duplicate terms 
  arrange(desc(Count),p.adjust) %>% 
  ggplot(aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_bar(aes(fill=pvalue),stat = "identity") +
  scale_fill_viridis_c(begin = 0.3, end = 1,option="magma") +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0,0))+
  facet_grid(ONTOLOGY~.,scales = "free")+
  theme_minimal(base_size = 20) + 
  xlab("Rich factor") +
  ylab("GO Term") + 
  labs(fill="p value")+
  ggtitle("")+
  guides(color=guide_colorbar(order=2))

## Note: these plots could be combined by adobe Adobe Illustrator


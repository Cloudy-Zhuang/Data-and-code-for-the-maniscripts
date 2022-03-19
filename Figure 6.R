## Figure 6
library(enrichplot)
require(clusterProfiler)



## get the work path
path_initial<- getwd()
cleandata_path <- paste0(path_initial,"/2. clean data")

## set work place

setwd(cleandata_path)

## load data
load("sub-network_KEGG_GO.Rdata")


## check the result
cluster_go_res <- cluster_go@result

##check the results based on different p value 
cluster_go_res  %>% 
  filter(pvalue<0.05)->x



# Figure 6: A-C -----------------------------------------------------------

## create a function to draw the plots
myempplot <- function(dat,term){
  
  dat %>% 
    filter(ONTOLOGY==term) %>% 
    arrange(desc(Count)) %>% 
    dplyr::slice(1:10) %>% 
    new("enrichResult",
        result=.) %>% 
    pairwise_termsim() ->ego
  
  p1 <- emapplot(ego , cex_label_category=.8, cex_line=.5) + 
    scale_fill_continuous(low = "#e06663", high = "#327eba", name = "pvalue",
                          guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')+
    labs(title = term,color="p value")+
    theme_void(base_size = 18)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  p1
}

x11()
#Figure 6: A
myempplot(x,"BP") 
#Figure 6: B
myempplot(x,"CC") 
#Figure 6: C
myempplot(x,"MF") 


# Figure 6D ---------------------------------------------------------------

library(ggraph)
library(tidygraph)

## load the NSCLC-DEGs to add the log2Foldchange
## information for the genes

load("TCGA-NSCLC_diffgene.Rdata")

##create two helpful functions to creat the network plot of
## the KEGG pathway term and enrichment genes

if (T) {
  gather_graph_edge <- function(df,index=NULL,root=NULL){
    require(dplyr)
    if (length(index) < 2){
      stop("please specify at least two index column(s)")
    } else if (length(index)==2){
      data <- df %>% mutate(from=.data[[index[[1]]]]) %>%
        tidyr::unite(to,index,sep="/") %>%
        select(from,to) %>%
        mutate_at(c("from","to"),as.character)
    } else {
      list <- lapply(seq(2,length(index)), function(i){
        dots <- index[1:i]
        df %>% tidyr::unite(from,dots[-length(dots)],sep = "/",remove = F)  %>%
          tidyr::unite(to,dots,sep="/") %>%
          select(from,to) %>%
          mutate_at(c("from","to"),as.character)
      })
      data <- do.call("rbind",list)
    }
    data <- as_tibble(data)
    if (is.null(root)){
      return(data)
    } else {
      root_data <- df %>% group_by(.dots=index[[1]]) %>%
        summarise(count=n()) %>%
        mutate(from=root,to=as.character(.data[[index[[1]]]] )) %>%
        select(from,to)
      rbind(root_data,data)
    }
    
  }
  
  
  gather_graph_node <- function(df,index=NULL,value=tail(colnames(df),1),root=NULL){
    require(dplyr)
    if (length(index) < 2){
      stop("please specify at least two index column(s)")
    } else {
      list <- lapply(seq_along(index), function(i){
        dots <- index[1:i]
        df %>%
          group_by(.dots=dots) %>%
          summarise(node.size=sum(.data[[value]]),
                    node.level=index[[i]],
                    node.count=n()) %>%
          mutate(node.short_name=as.character(.data[[ dots[[length(dots)]] ]]),
                 node.branch = as.character(.data[[ dots[[1]]]])) %>%
          tidyr::unite(node.name,dots,sep = "/")
      })
      data <- do.call("rbind",list) %>% as_tibble()
      data$node.level <- factor(data$node.level,levels = index)
      
      if (is.null(root)){
        return(data)
      } else {
        root_data <- data.frame(node.name=root,
                                node.size=sum(df[[value]]),
                                node.level=root,
                                node.count=1,
                                node.short_name=root,
                                node.branch=root,
                                stringsAsFactors = F)
        data <- rbind(root_data,data)
        data$node.level <- factor(data$node.level, levels = c(root,index))
        return(data)
      }
    }
  }  
  
}


cluster_EGG@result %>% 
  filter(p.adjust<0.05) %>%  # filter the terms with p.djust<0.05
  dplyr::select(Description, geneID) %>% 
  dplyr::rename("term"=Description,
                "genes"=geneID) %>% 
  separate_rows(genes,sep="/") %>% 
  inner_join(TCGA_NSCLC_diffgene %>% 
               dplyr::rename("genes"=gene_id,
                             "logFC"=log2FoldChange) %>% 
             select(genes,logFC))->
  df

df$genes %>% unique()
geneSpecial <- data.frame(Gene=unique(df$genes),
                          Type=paste0("type",seq(1,n_distinct(df$genes))))
geneCol <- geneSpecial$Type
names(geneCol) <- geneSpecial$Gene
geneCol

nodes <- gather_graph_node(df, 
                            index = c("term", "genes"),
                            value = "logFC", 
                            root="all")
edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
nodes <- nodes %>% mutate_at(c("node.level","node.branch"),as.character)
head(nodes, 10)




nodes$color <- "normal"
nodes[nodes$node.short_name %in% geneSpecial$Gene,]$color <- geneCol[nodes[nodes$node.short_name %in% geneSpecial$Gene,]$node.short_name]
nodes[nodes$node.short_name %in% geneSpecial$Gene,]

nodes$color <- factor(nodes$color, levels = unique(nodes$color))

# With node and edge data, 'tbl_graph()' produces a graph.  
graph <- tbl_graph(nodes, edges)


gc1 <-  ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  
  geom_edge_diagonal(aes(color = node2.color,
                         filter=node1.node.level!="all"), 
                     alpha = 0.5, 
                     edge_width=2.5) + 
  geom_node_point(aes(size = node.size,
                      filter=node.level!="all"), 
                  #alpha = 1/3,
                  color = "#61C3ED") +
  scale_size(range = c(0.5,30)) +
  theme(legend.position = "none") + 
  
  # Add the surrounding comment text, in this case they are gene  name
  geom_node_text(
    aes(
      x = 1.05 * x, 
      y = 1.05 * y, 
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf
    ),
    color="black", 
    size = 6, hjust = 'outward') +
  
 # add inner ring text, here is the pathway name 
  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all")
    ),
    color="black", 
    fontface="bold",
    size=6,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) + 
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) 
x11()
gc1



## Note: these vene plots could be combined by adobe Adobe Illustrator







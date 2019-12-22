data_biogrid<-read.table("/home/wuzefeng/MyResearch/networks/PPI/Arabidopsis/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-3.4.134.tab2.txt",header = TRUE,sep = "\t", check.names=FALSE)
data_biogrid<-data_biogrid[c(6,7)]
data<-as.data.frame(t(apply(data_biogrid,1,sort)))
data1<-unique(data)
data2<-subset(data1,as.character(data1$V1)!=as.character(data1$V2))
library(igraph)
g<-graph_from_data_frame(data2,directed = FALSE)
V(g)$size = 0.2
plot.igraph(g,vertex.label=NA,layout=layout.fruchterman.reingold)
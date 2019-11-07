#! /usr/bin/R

set.seed(1000)
options(digits = 5)

df2network<-function(df){
  data<-read.table(df,header = TRUE,stringsAsFactors = FALSE,sep="\t")
  library(igraph)
  g<-graph.data.frame(d = data, directed = FALSE)
  g<-simplify(g,remove.multiple = TRUE,remove.loops = TRUE)
  return(g)
}

#### read network data with different cufoff 

g<-df2network("myresearch/network/1network_cutoff10.txt")  

### degree kept network
gg<-rewire(g,with = keeping_degseq(niter = vcount(g) * 100))


## network attributions calculation
#0.0 network diameter

#farthest.nodes(g,weights = NA) # AT4G04200 AT5G67580; 807.44

# network betweeness
system.time(AGFN_betweenness<-betweenness(g,directed = FALSE,weights = NA,normalized = TRUE))
write.table(AGFN_betweenness)

### Fit power-law
data <- degree(g)
data.dist <- data.frame(k=0:max(data),p_k=degree_distribution(g))
data.dist <- data.dist[data.dist$p_k>0,]
library(ggplot2)
ggplot(data.dist) + geom_point(aes(x=k, y=p_k)) +  # log(k)
  theme_bw()+
  theme(text=element_text(size = 20)) + ylab("Probability of degree")+xlab("Degree")

library(poweRlaw)
m_pl <- displ$new(data)
est_pl <- estimate_xmin(m_pl)
m_pl$setXmin(est_pl)
plot.data <- plot(m_pl, draw = F)
fit.data <- lines(m_pl, draw = F)

ggplot(plot.data) + geom_point(aes(x=log(x), y=log(y))) + 
                  labs(x="log(k)", y="log(CDF)") + theme_bw() + 
                  geom_line(data=fit.data, aes(x=log(x), y=log(y)), colour="red")+
                  theme(text=element_text(size=20))
### network closeness
system.time(AGFN_closeness<-closeness(g,directed = FALSE,weights = NA,normalized = TRUE))

###1.1 network validation using tair interaction data
tair_inter<-read.table("/home/wuzefeng/MyResearch/networks/PPI/Arabidopsis/TairProteinInteraction.20090527.txt",header = TRUE,stringsAsFactors = FALSE,sep="\t")
focused_genes_list <-unique(c(tair_inter$Locus_name,tair_inter$InteractorLocus_name))
for (x in focused_genes_list){
  focused_genes<-x
  if (x %in% names(V(g))){
  sub_g<-induced_subgraph(g, ego(g, 1, focused_genes)[[1]])
  
  focus_links<-subset(tair_inter,tair_inter$Locus_name==focused_genes | tair_inter$InteractorLocus_name==focused_genes )
  focus_links<-focus_links[,c(1,3,5,6)]
  focus_links<-unique(as.data.frame(t(apply(focus_links[,c(1,2)],1,sort)),stringsAsFactors = FALSE))
  colnames(focus_links)<-c("V1","V2")
  focus_links<-subset(focus_links,focus_links$V1!=focus_links$V2)
  message(c(x,";","experiment links number is: ",dim(focus_links)[1],";","of these the predicted links is: ",length(intersect(names(V(sub_g)),unique(as.vector(t(as.matrix(focus_links))))))-1))
  write.table(file="4network_analysis_result/1network_validation/tair_interaction_validation",data.frame(gene=x,experiment=dim(focus_links)[1],predicted=length(intersect(names(V(sub_g)),unique(as.vector(t(as.matrix(focus_links))))))-1),sep="\t",append=TRUE,col.names = !file.exists("out"),row.names = FALSE,quote = FALSE)
  focus_links<-subset(focus_links,focus_links$V2%in%names(V(sub_g))&focus_links$V1%in%names(V(sub_g))) # common from exper and predicte 
  }
} 
## result plot
tair<-read.table("4network_analysis_result/1network_validation/tair_interaction_validation",header = TRUE)
tair<-subset(tair,tair$experiment>0)
tair$accuaracy<-tair$predicted/tair$experiment

### compare with permute g  for 20 times
for(m in 1:20){
  message(m)
  permute_g<-rewire(g,with = keeping_degseq(niter = vcount(g) * 100))
  tair_inter<-read.table("/home/wuzefeng/MyResearch/networks/PPI/Arabidopsis/TairProteinInteraction.20090527.txt",header = TRUE,stringsAsFactors = FALSE,sep="\t")
  focused_genes_list <-unique(c(tair_inter$Locus_name,tair_inter$InteractorLocus_name))
  for (x in focused_genes_list){
    focused_genes<-x
    if (x %in% names(V(permute_g))){
      sub_g<-induced_subgraph(permute_g, ego(permute_g, 1, focused_genes)[[1]])
      
      focus_links<-subset(tair_inter,tair_inter$Locus_name==focused_genes | tair_inter$InteractorLocus_name==focused_genes )
      focus_links<-focus_links[,c(1,3,5,6)]
      focus_links<-unique(as.data.frame(t(apply(focus_links[,c(1,2)],1,sort)),stringsAsFactors = FALSE))
      colnames(focus_links)<-c("V1","V2")
      focus_links<-subset(focus_links,focus_links$V1!=focus_links$V2)
      message(c(x,";","experiment links number is: ",dim(focus_links)[1],";","of these the predicted links is: ",length(intersect(names(V(sub_g)),unique(as.vector(t(as.matrix(focus_links))))))-1))
      write.table(file=paste("4network_analysis_result/1network_validation/TIP_validation_20_random/tair_interaction_validation.random",m,sep=""),
                  data.frame(gene=x,experiment=dim(focus_links)[1],
                             predicted=length(intersect(names(V(sub_g)),unique(as.vector(t(as.matrix(focus_links))))))-1),
                  sep="\t",append=TRUE,
                  col.names = !file.exists(paste("4network_analysis_result/1network_validation/TIP_validation_20_random/tair_interaction_validation.random",m,sep="")),
                  row.names = FALSE,quote = FALSE)
      
      #focus_links<-subset(focus_links,focus_links$V2%in%names(V(sub_g))&focus_links$V1%in%names(V(sub_g))) 
    }
  } 
}
## plot
random <- c()
# create directory names
random_dir <- dir("/home/wuzefeng/MyResearch/networks/2network_prediction/4network_analysis_result/1network_validation/TIP_validation_20_random",full.names = TRUE)
# loop through all directories and grab fpkm columns
for( i in 1:length(random_dir) ){
  fname <- random_dir[i]
  x <- read.table(file=fname, sep="\t", header=T, as.is=T)
  random <- cbind(random, x[,"predicted"]/x[,"experiment"])
}
# name the columns
colnames(random) <- stringr::str_split_fixed(basename(random_dir),pattern = "\\.",n = 2)[,2]
# name the rows, they're all in the same order
rownames(random) <- x[,1]

tair_random<-data.frame(apply(random,1,mean))
colnames(tair_random)<-"random_accuaracy"
tair_random<-subset(tair_random,tair_random$random_accuaracy>=0)

dd1<-data.frame(accuaracy=tair$accuaracy,class="Original")
dd2<-data.frame(accuaracy=tair_random$random_accuaracy,class="Permutation")
dd<-rbind(dd1,dd2)
ggplot(dd,aes(x=accuaracy,fill=class))+geom_histogram(position = "dodge")+
                                       theme(text = element_text(size = 20),legend.position=c(0.85,0.9))+
                                        ylab("Count")+xlab("Predictive accuracy")+
                                        scale_fill_manual("Network type",values=c('#9ED2F0','#E6A429'))

#### 1.2 individule target gene validation and plot
tair_inter<-read.table("/home/wuzefeng/MyResearch/networks/PPI/Arabidopsis/TairProteinInteraction.20090527.txt",header = TRUE,stringsAsFactors = FALSE,sep="\t")
focused_genes<-"AT2G18790" # AT1G02340 #AT2G18790
sub_g<-induced_subgraph(g, ego(g, 1, focused_genes)[[1]])

focus_links<-subset(tair_inter,tair_inter$Locus_name==focused_genes | tair_inter$InteractorLocus_name==focused_genes)
focus_links<-focus_links[,c(1,3,5,6)]
focus_links<-unique(as.data.frame(t(apply(focus_links[,c(1,2)],1,sort)),stringsAsFactors = FALSE))
focus_links<-subset(focus_links,focus_links$V1!=focus_links$V2)
message(c("experiment links number is: ",dim(focus_links)[1],";","of these the predicted links is: ",length(intersect(names(V(sub_g)),unique(as.vector(t(as.matrix(focus_links))))))-1))
focus_links<-subset(focus_links,focus_links$V2%in%names(V(sub_g))&focus_links$V1%in%names(V(sub_g))) 

### plot attribution
plot_personal_graph<-function(graph_object,focused_list){
  #V(sub_g)$color[names(V(sub_g))%in%focus_links$V2]<-"purple"
  V(sub_g)$color<-ifelse(names(V(sub_g))%in%focused_list,"red","steelblue") # vertex color 
  V(sub_g)$size <- 8              # vertex size
  V(sub_g)$label.cex <- 1          # vertex label size
  V(sub_g)$label.color<-"brown"
  E(sub_g)$color <- "gray"        # edge color 
  E(sub_g)$width=2#E(sub_g)$weight/max(E(sub_g)$weight)
  E(sub_g,P=as.vector(t(as.matrix(focus_links))))$color<-"purple"
  E(sub_g,P=as.vector(t(as.matrix(focus_links))))$width<-5
  #E(sub_g)$width=edge.betweenness(sub_g)
  plot.igraph(sub_g,layout=layout.fruchterman.reingold)
  legend('topleft',
         legend=c("Predicted","Validation"),
         pch=95, #shape
         box.lty=2, # 
         lty=1:2,
         pt.cex= 3, #lines size 
         cex=1, #box size
         col=c("gray","purple"))
}
plot_personal_graph(graph_object = sub_g,focused_list = focused_genes)



## 2.kegg validation
###2.1  network validation using kegg pathway
require(org.At.tair.db)
xx<- as.list(org.At.tairPATH2TAIR)

combins<-combn(seq(1:length(xx)),2)
for (col in 1:ncol(combins)){
  gene_list1<-xx[[combins[,col][1]]]
  gene_list2<-xx[[combins[,col][2]]]
  
  group_in_possible<-choose(length(gene_list1),2)+choose(length(gene_list2),2)-choose(length(intersect(gene_list1,gene_list2)),2)
  group_in_prediction <- ecount(induced_subgraph(g,intersect(gene_list1,names(V(g)))) %u% induced_subgraph(g,intersect(gene_list2,names(V(g)))))
  
  between_group_possible <-length(gene_list1)*length(gene_list2)-length(intersect(gene_list1,gene_list2))^2+choose(length(intersect(gene_list1,gene_list2)),2)
  between_group_prediction <-ecount(induced_subgraph(g,intersect(unique(c(gene_list1,gene_list2)),names(V(g)))))-group_in_prediction
  
  print(c(col, group_in_possible,group_in_prediction,between_group_possible,between_group_prediction))
  message(fisher.test(matrix(c(group_in_prediction,group_in_possible,between_group_prediction,between_group_possible),nrow=2),alternative="greater")$p.value)
  write.table(file="KEGG_validation",data.frame(pathway=paste(names(xx)[combins[,col][1]],names(xx)[combins[,col][2]],sep = "_"),
                                                predicted_in_pathway = group_in_prediction,
                                                in_pathway = group_in_possible,
                                                prediction_between_pathway = between_group_prediction,
                                                between_pathway = between_group_possible,
                                                p_value = fisher.test(matrix(c(group_in_prediction,group_in_possible,between_group_prediction,between_group_possible),nrow=2),alternative="greater")$p.value),
                                                sep="\t",append=TRUE,col.names = !file.exists("KEGG_validation"),
                                                row.names = FALSE,quote = FALSE)
  }


####2.2 reanalysis kegg data with permutations network

for (m in 1:20){
message(m)
permute_g<-rewire(g,with = keeping_degseq(niter = vcount(g) * 100))

for (col in 1:ncol(combins)){
  gene_list1<-xx[[combins[,col][1]]]
  gene_list2<-xx[[combins[,col][2]]]
  
  group_in_possible<-choose(length(gene_list1),2)+choose(length(gene_list2),2)-choose(length(intersect(gene_list1,gene_list2)),2)
  group_in_prediction <- ecount(induced_subgraph(permute_g,intersect(gene_list1,names(V(permute_g)))) %u% induced_subgraph(permute_g,intersect(gene_list2,names(V(permute_g)))))
  
  between_group_possible <-length(gene_list1)*length(gene_list2)-length(intersect(gene_list1,gene_list2))^2+choose(length(intersect(gene_list1,gene_list2)),2)
  between_group_prediction <-ecount(induced_subgraph(permute_g,intersect(unique(c(gene_list1,gene_list2)),names(V(permute_g)))))-group_in_prediction
  
  print(c(col, group_in_possible,group_in_prediction,between_group_possible,between_group_prediction))
  message(fisher.test(matrix(c(group_in_prediction,group_in_possible,between_group_prediction,between_group_possible),nrow=2),alternative="greater")$p.value)
  write.table(file=paste("4network_analysis_result/1network_validation/KEGG_validation_20_random/KEGG_validation.random",m,sep=""),data.frame(pathway=paste(names(xx)[combins[,col][1]],names(xx)[combins[,col][2]],sep = "_"),
                                                predicted_in_pathway = group_in_prediction,
                                                in_pathway = group_in_possible,
                                                prediction_between_pathway = between_group_prediction,
                                                between_pathway = between_group_possible,
                                                p_value = fisher.test(matrix(c(group_in_prediction,group_in_possible,between_group_prediction,between_group_possible),nrow=2),alternative="greater")$p.value),
                                                sep="\t",append=TRUE,col.names = !file.exists(paste("4network_analysis_result/1network_validation/KEGG_validation_20_random/KEGG_validation.random",m,sep="")),
                                                row.names = FALSE,quote = FALSE)
}
}


##ggplot
data1<-read.table("4network_analysis_result/1network_validation/KEGG_validation",stringsAsFactors = FALSE,header = TRUE)
data1$qvalue<-p.adjust(data1$p_value)
### 
random <- c()
# create directory names
random_dir <- dir("/home/wuzefeng/MyResearch/networks/2network_prediction/4network_analysis_result/1network_validation/KEGG_validation_20_random",full.names = TRUE)
# loop through all directories and grab fpkm columns
for( i in 1:length(random_dir) ){
  fname <- random_dir[i]
  x <- read.table(file=fname, sep="\t", header=T, as.is=T)
  random <- cbind(random, x[,"p_value"])
}
# name the columns
colnames(random) <- stringr::str_split_fixed(basename(random_dir),pattern = "\\.",n = 2)[,2]
# name the rows, they're all in the same order
rownames(random) <- x[,1]

tair_random<-data.frame(apply(random,1,mean))
colnames(tair_random)<-"p_value"


p_values<-rbind(data.frame(Pvalue=data1$p_value,class="original",stringsAsFactors = FALSE),data.frame(Pvalue=tair_random$p_value,class="permutation",stringsAsFactors = FALSE))
library(ggplot2)
p<-ggplot(p_values,aes(x=Pvalue,fill=class))+geom_histogram(position = "dodge")+
      theme(text = element_text(size=20),legend.position=c(0.85,0.9))+
      ylab("Frequency")+xlab("P-value")+
      scale_fill_manual("Network type",values=c('#9ED2F0','#E6A429')) #legend modifications


#####
####3 compare wit aranet2
ara_net2 <-read.table("other_sourse/AraNet.txt",header = FALSE,stringsAsFactors = FALSE)
colnames(ara_net2)<-c("X1","X2","weight")
ara_net2<-ara_net2[!duplicated(ara_net2[c(1,2)]),]
library(igraph)
ara_net_g<-graph.data.frame(d = ara_net2, directed = FALSE)
ara_net_g<-simplify(ara_net_g,remove.multiple = TRUE,remove.loops = TRUE) #22894 *895000 

intersect_g<-ara_net_g%s%g #14112 common edges
intersect_permute_g<-permute_g%s%ara_net_g #199

### 3.5 compare with AtPIN

AtPIN <-read.table("other_sourse/AtPIN_PPI.txt",header = FALSE,stringsAsFactors = FALSE,sep="\t")
AtPIN<-AtPIN[,-3]
colnames(AtPIN)<-c("X1","X2")
AtPIN<-AtPIN[!duplicated(AtPIN[c(1,2)]),]
library(igraph)
AtPIN_g<-graph.data.frame(d = AtPIN, directed = FALSE)
AtPIN_g<-simplify(AtPIN_g,remove.multiple = TRUE,remove.loops = TRUE) #15163 * 95043

intersect_g<-ara_net_g%s%g # 14112
intersect_permute_g<-permute_g%s%ara_net_g #1824


###4 analysis flowering geneset
##nalysis flowering geneset

flowering_genes<-read.csv("other_resources/Flowering Interactive Database FLOR-ID - Flowering time.csv",header = TRUE,stringsAsFactors = FALSE)
commen_flower_genes<- intersect(flowering_genes$Gene.details,names(V(g)))

#short distence among these flowering genes
d_flower<-shortest.paths(g,commen_flower_genes, commen_flower_genes,weights = NA)
hist(d_flower)

#random sampled genes

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),285,replace = FALSE)
  d_random<-shortest.paths(g,rand_genes, rand_genes,weights = NA)
  return(mean(as.numeric(d_random)[is.finite(as.numeric(d_random))]))
}                              
library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 7))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue")+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean shortest path")+
  annotate("segment", x=2.95, xend=2.95, y=25, yend=0, color="black", size=2, arrow=arrow())
  

###4.5 flc genes
####4.5 flc genes

sub_g<-induced_subgraph(g, ego(g, 1, "AT5G10140")[[1]])
tair_no_annotated <-c("AT1G02840","AT1G32320","AT1G54440","AT2G20050","AT2G26330","AT2G35110","AT3G26790","AT3G46520","AT5G60410","AT5G45830")

V(sub_g)$size <- 8              # vertex size
V(sub_g)$label.cex <- 1          # vertex label size
V(sub_g)$label.color<-"brown"
E(sub_g)$color <- "gray"        # edge color 
V(sub_g)$color<-"steelblue"
V(sub_g)$color[names(V(sub_g))=="AT5G10140"]<-"orange"
V(sub_g)$color[names(V(sub_g))%in%tair_no_annotated]<-"lightgreen"

plot(sub_g,vertex.frame.color= "white")
legend('topleft',
       legend=c("FLC","Known","Unknown"),
       pch=19, #shape
       box.lty=2, # 
       pt.cex= 3, #lines size 
       cex=1, #box size
       col=c("orange","steelblue","lightgreen"),
       y.intersp=1.5
       )

## gene degree and dn/ds
genes_degree<-degree(g)
dnds<-read.table("myresearch/network/data/TAIR2lyrata.ds.mart_export.txt",sep="\t",header = TRUE,stringsAsFactors = FALSE)
dnds<-dnds[complete.cases(dnds),]
dnds<-unique(dnds)
dnds<-dnds %>% group_by(Gene.stable.ID) %>% slice(which.min(dS.with.Arabidopsis.lyrata)) # select minumum ds for same genes 
dnds$dnds<-dnds$dN.with.Arabidopsis.lyrata/dnds$dS.with.Arabidopsis.lyrata
rownames(dnds)<-dnds$Gene.stable.ID

genes_dnds<-dnds[names(degree(g)),]$dnds
df<-data.frame(degree=degree(g),dnds=genes_dnds)
df<-na.omit(df)
df<-df[is.finite(df$dnds),]
df<-df[df$dnds<5,]
cor.test(df$degree,df$dnds) # pcc = -0.12



##7 imprinted genes analysis
###7.0 import imprinted genes 
Arabidospis_imprinted_genes_parse<-function(ara_imprinted_genes_file){
  imprinted_data <-read.table(ara_imprinted_genes_file,stringsAsFactors = FALSE)
  imprinted_genes<-unique(imprinted_data$V1)
  imprinted_genes_in_network<- intersect(imprinted_genes,names(V(g)))
  
  if (ncol(imprinted_data)==2){
    paternal_imprint <- unique(imprinted_data$V1[imprinted_data$V2=="f"])
    maternal_imprint <- unique(imprinted_data$V1[imprinted_data$V2=="m"]) 
  }
  if (ncol(imprinted_data)==3){
    paternal_imprint <- unique(imprinted_data$V1[imprinted_data$V3=="f"])
    maternal_imprint <- unique(imprinted_data$V1[imprinted_data$V3=="m"]) 
  }
  
  IGs<-list()
  IGs$IG <- imprinted_genes
  IGs$PEG <- paternal_imprint
  IGs$MEG <- maternal_imprint
  
  IGs$IG_in_networks <- imprinted_genes_in_network
  IGs$PEG_in_networks <- intersect(paternal_imprint,names(V(g)))
  IGs$MEG_in_networks <- intersect(maternal_imprint,names(V(g)))
  message("Imprinted genes number in Arabidopsis is: ",length(IGs$IG))
  return(IGs)
}

imprinted_data <-Arabidospis_imprinted_genes_parse("myresearch/network/data/imp2+.list")
#imprinted_data <-Arabidospis_imprinted_genes_parse("/home/wuzefeng/MyResearch/Imprinting_prediction/imprint_gene_list/TAIR10/6imprinted.list") #526

imprinted_genes<-imprinted_data$IG  
imprinted_genes_in_network<-imprinted_data$IG_in_networks ### 82|470
maternal_imprint<-imprinted_data$MEG_in_networks  ###  #376
paternal_imprint<-imprinted_data$PEG_in_networks  ###  #94

sub_g<-induced_subgraph(g, imprinted_genes_in_network)

V(sub_g)$color[names(V(sub_g))%in%maternal_imprint]<-"tomato" # vertex color 
V(sub_g)$color[names(V(sub_g))%in%paternal_imprint]<-"steelblue"
V(sub_g)$distance <- ifelse(names(V(sub_g))%in%componet_IGs,3,iso2comp_min_dis$Freq[match(names(V(sub_g)),iso2comp_min_dis$Var1)])              # vertex size
V(sub_g)$label.cex <- 0.8          # vertex label size
V(sub_g)$label.color<-"black"
E(sub_g)$color <- "gray"        # edge color 
E(sub_g)$width=3
V(sub_g)$class<-ifelse(names(V(sub_g))%in%maternal_imprint,"MEG","PEG")
 
## plot
library(intergraph)
library(ggnetwork)

dat <- ggnetwork(sub_g, layout="fruchtermanreingold", arrow.gap=0, cell.jitter=0)

ggplot(dat,aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_edges(color="grey50", curvature=0.1, size=0.5) +
  geom_nodes(aes(size=distance,color=class)) +
  geom_nodetext(aes(label=vertex.names),size=3, color="#8856a7") +
  theme_blank()+scale_color_brewer(palette = "Set2")+
  scale_size_area(max_size = 9)



## 7.0.1imprinted and neighbors netwoerk

sub_g<-induced_subgraph(g, unique(names(unlist(ego(g,order = 1,imprinted_genes_in_network)))))
plot_personal_graph<-function(graph_object){
  V(sub_g)$color[names(V(sub_g))%in%maternal_imprint]<-"tomato" # vertex color 
  V(sub_g)$color[names(V(sub_g))%in%paternal_imprint]<-"steelblue"
  V(sub_g)$color[!names(V(sub_g))%in%imprinted_genes_in_network]<-"black"
  V(sub_g)$size <- 2              # vertex size
  V(sub_g)$label.cex <- 0.8          # vertex label size
  V(sub_g)$label.color<-"black"
  E(sub_g)$color <- "gray"        # edge color 
  E(sub_g)$width=1
  
  plot.igraph(sub_g,layout=layout.fruchterman.reingold,vertex.frame.color= "white",vertex.label=NA)
  legend('topleft',
         legend=c("Maternal","Paternal","Partners"),
         pch=19, #shape
         box.lty=2, # 
         pt.cex= 3, #lines size 
         cex=1, #box size
         col=c("tomato","steelblue","black"),
         y.intersp=1.5)
}
plot_personal_graph(graph_object = sub_g)


### perform network module analysis
nd <- as.data.frame(get.edgelist(sub_g),stringsAsFactors = FALSE)
colnames(nd)<-c("X1","X2")
gene_list<-unique(c(nd$X1,nd$X2))
## make gene2num table 
gene_list_numbers <-seq(0,length(gene_list)-1)
gene_name2number<-data.frame(gene_list,gene_list_numbers,stringsAsFactors = FALSE)

## assign network vetex to number
nd$g1<-gene_name2number$gene_list_numbers[match(nd$X1,gene_name2number$gene_list)]
nd$g2<-gene_name2number$gene_list_numbers[match(nd$X2,gene_name2number$gene_list)]

## import imprinted genes
imprinted_genes<-read.table("~/MyResearch/Imprinting_prediction/imprint_gene_list/TAIR10/ara_imp_by_paper_num/imp2+.list",stringsAsFactors = FALSE)
imprinted_genes$color<-ifelse(imprinted_genes$V2=="m",1,2)

###
nd$c1<- imprinted_genes$color[match(nd$X1,imprinted_genes$V1)]
nd$c2<- imprinted_genes$color[match(nd$X2,imprinted_genes$V1)]
nd$c1[is.na(nd$c1)]=0
nd$c2[is.na(nd$c2)]=0

fanmond_input<-nd[,c(3,4,5,6)]
write.table(fanmond_input,"fanmod_inout.txt",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")


imprinted_genes<-read.table("~/MyResearch/Imprinting_prediction/imprint_gene_list/TAIR10/ara_imp_by_paper_num/imp2+.list",stringsAsFactors = FALSE)
nd$c1<- imprinted_genes$color[match(nd$X1,imprinted_genes$V1)]



## 7.1 imp-imp distence
imp_distance <-shortest.paths(g,imprinted_genes_in_network, imprinted_genes_in_network,weights = NA)
imp_distance <- imp_distance[upper.tri(imp_distance)]

## 7.1.5 random simulation

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(imprinted_genes_in_network),replace = FALSE)
  d_random<-shortest.paths(g,rand_genes, rand_genes,weights = NA)
  d_random <-d_random[upper.tri(d_random)]
  d_random <-d_random[is.finite(d_random)]
  return(mean(d_random))
}                              
library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 60))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 200)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean shortest path")+
  annotate("segment", x=mean(imp_distance[is.finite(imp_distance)]), xend=mean(imp_distance[is.finite(imp_distance)]), y=5, yend=0, color="black", size=2, arrow=arrow())+
  annotate("segment", x=median(imp_distance[is.finite(imp_distance)]), xend=median(imp_distance[is.finite(imp_distance)]), y=5, yend=0, color="red", size=2, arrow=arrow())+
  annotate("text", x = median(imp_distance[is.finite(imp_distance)]), y = 6, label = "Median")+
  annotate("text", x = mean(imp_distance[is.finite(imp_distance)]), y = 6, label = "Mean")


### paternal and maternal distance

maternal_imp_distance <-shortest.paths(g,maternal_imprint, maternal_imprint,weights = NA)
maternal_imp_distance<-maternal_imp_distance[upper.tri(maternal_imp_distance)]
maternal_imp_distance<-maternal_imp_distance[is.finite(maternal_imp_distance)]

paternal_imp_distance <-shortest.paths(g,paternal_imprint, paternal_imprint,weights = NA)
paternal_imp_distance<-paternal_imp_distance[upper.tri(paternal_imp_distance)]
paternal_imp_distance<-paternal_imp_distance[is.finite(paternal_imp_distance)]

maternal_paternal_distance <- shortest.paths(g,paternal_imprint, maternal_imprint,weights = NA)
maternal_paternal_distance<-maternal_paternal_distance[upper.tri(maternal_paternal_distance)]
maternal_paternal_distance<-maternal_paternal_distance[is.finite(maternal_paternal_distance)]

dis_maternal<-data.frame(shortest.paths=maternal_imp_distance,class="MEG-MEG")
dis_paternal<-data.frame(shortest.paths=paternal_imp_distance,class="PEG-PEG")
dis_maternal2paternal<-data.frame(shortest.paths=maternal_paternal_distance,class="MEG-PEG")
df<-rbind(dis_maternal,dis_paternal,dis_maternal2paternal)
df$shortest.paths<-ifelse(df$shortest.paths>5,5,df$shortest.paths)


ggplot(df,aes(y=shortest.paths,x=class,fill=class))+
      geom_violin()+
      theme_bw(base_size = 20)+
      geom_signif(comparisons = list(c("MEG", "PEG"),c("MEG","MEG->PEG"),c("MEG->PEG","PEG")),
      test = "wilcox.test",map_signif_level = "FALSE",
      test.args = list(alternative = "greater"),
      step_increase = 0.05,
      tip_length = 0.01)+
      scale_fill_manual(values = c("tomato","steelblue", "orange")) +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
      ylab("Shortest path") +
      stat_summary(fun.y = mean, geom = "point", size = 2, color = "red")

### optional plot
df<-df%>%group_by(class,shortest.paths)%>%summarise(n=n())%>%mutate(y=n/sum(n))
ggplot(df, aes(fill=as.factor(shortest.paths), y=y, x=class))+
  theme_bw()+
  theme(text=element_text(size=20))+
  geom_bar( stat="identity", position="fill")+
  scale_fill_brewer(name="Shortest path",palette = "Set3",labels=c(c(1,2,3,4),expression(phantom(x)>=5)))+
  geom_text(aes(label=paste(round(y,2)*100,"%",sep="")),color="black", size=3.5,position=position_fill(0.5))+
  ylab("Percentage")+xlab("")
  #scale_fill_manual(name="Shortest path", values=c("#458B74", "#CDAA7D", "#8968CD", "#CD5555", "#1874CD", "#EE7600"))

### keeping degree unchanged to count the frequency of the largestest component by sampling same number of genes in random networks
IGON<-induced_subgraph(g, imprinted_genes_in_network)
max_component_size_of_IGON<-max(components(IGON)$csize)
max_component_size_of_random<-c()  
for (m in (seq(1,1000))){
  sub_gg<-induced_subgraph(g,vids = sample(names(V(g)),length(imprinted_genes_in_network)))
  max_component_size_of_sub_gg<-max(components(sub_gg)$csize)
  max_component_size_of_random<-c(max_component_size_of_random,max_component_size_of_sub_gg)
}
ggplot(data.frame(Max_component_size = max_component_size_of_random),aes(x=Max_component_size))+
        geom_histogram(bins = 200,fill="steelblue")+
        theme_bw(base_size = 20)+
        xlab("Maximal component size")+ylab("Count")+
        #xlim(0,45)+
        annotate("segment",x=max_component_size_of_IGON,xend = max_component_size_of_IGON,yend=0,y=5,size=2, arrow=arrow(),color="red")+
        annotate("text", x = max_component_size_of_IGON, y = 6, label = "Imprinted genes")


### calulate isolated genes to component genes distance
iso_IGs <-names(components(IGON)$membership)[components(IGON)$membership!=1] #48 IGs
componet_IGs <- names(components(IGON)$membership)[components(IGON)$membership==1] #34
iso2comp_dis<- shortest.paths(g,v = iso_IGs,to = componet_IGs,weights = NA)
iso2comp_dis_df<-as.data.frame(as.table(iso2comp_dis))
iso2comp_min_dis<-data.frame(iso2comp_dis_df%>%group_by(Var1)%>%summarise_at("Freq",min))
distance1_ratio<-sum(iso2comp_min_dis$Freq==2)/nrow(iso2comp_min_dis) #60.4
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.00    2.00    2.00    2.42    3.00    4.00

### random sampled genes to component distance (nort significant)

fun <- function(x){
  require(igraph)
  random_genes<-sample(names(V(g))[!names(V(g))%in%componet_IGs],length(iso_IGs))
  iso2comp_dis<- shortest.paths(g,v = random_genes,to = componet_IGs,weights = NA)
  iso2comp_dis_df<-as.data.frame(as.table(iso2comp_dis))
  a<-data.frame(iso2comp_dis_df%>%group_by(Var1)%>%summarise_at("Freq",min))
  return(mean(a$Freq))
}                              
library(parallel)
system.time(random2comp_dis<-mclapply(1:1000,fun,mc.cores = 60))
random2comp_dis<-unlist(random2comp_dis)
ggplot(data.frame(shortest.path = random2comp_dis[is.finite(random2comp_dis)]),aes(x=shortest.path))+
  geom_histogram(fill="steelblue")+
  theme_bw(base_size = 20)+
  xlab("Maximal component size")+ylab("Count")+
  #xlim(0,45)+
  annotate("segment",x=mean(iso2comp_min_dis$Freq),xend = mean(iso2comp_min_dis$Freq),yend=0,y=25,size=2, arrow=arrow(),color="red")+
  annotate("text", x = mean(iso2comp_min_dis$Freq), y = 30, label = "Imprinted genes")


##### the ratio of non-component genes to IG component genes distance is 1 in random network 
 
fun<-function(x){
  random_genes<-sample(names(V(g))[!names(V(g))%in%componet_IGs],length(iso_IGs))
  iso2comp_dis<- shortest.paths(g,v = random_genes,to = componet_IGs,weights = NA)
  iso2comp_dis_df<-as.data.frame(as.table(iso2comp_dis))
  a<-data.frame(iso2comp_dis_df%>%group_by(Var1)%>%summarise_at("Freq",min))
  return(sum(a$Freq==2)/nrow(a))
}
system.time(random2comp_dis1_ratio<-mclapply(1:1000,fun,mc.cores = 60))
random2comp_dis1_ratio<-unlist(random2comp_dis1_ratio)

ggplot(data.frame(Max_component_size = random2comp_dis1_ratio),aes(x=Max_component_size))+
  geom_histogram(fill="steelblue",bins=200)+
  theme_bw(base_size = 14)+
  xlab("Pecentage of genes seperated by one genes to component")+ylab("Count")+
  #xlim(0,45)+
  annotate("segment",x=distance1_ratio,xend = distance1_ratio,yend=0,y=5,size=2, arrow=arrow(),color="red")+
  annotate("text", x = distance1_ratio, y = 6, label = "IGON")


### non-component genes to component gene distance1 in random networks

fun<-function(x){
  sub_gg<-induced_subgraph(g,vids = sample(names(V(g)),length(imprinted_genes_in_network)))
  componet_Gs <-names(V(sub_gg)[which.max(components(sub_gg)$csize)==components(sub_gg)$membership]) #48 IGs
  iso_Gs <- names(V(sub_gg)[which.max(components(sub_gg)$csize)!=components(sub_gg)$membership])
  iso2comp_dis<- shortest.paths(g,v = iso_Gs,to = componet_Gs,weights = NA)
  iso2comp_dis_df<-as.data.frame(as.table(iso2comp_dis))
  b<-data.frame(iso2comp_dis_df%>%group_by(Var1)%>%summarise_at("Freq",min))
  ratio_distance1<-sum(b$Freq==2)/nrow(b)
  return(ratio_distance1)
}
system.time(noncomp2comp_distance1_ratio<-mclapply(1:1000,fun,mc.cores = 60))
noncomp2comp_distance1_ratio<-unlist(noncomp2comp_distance1_ratio)


ggplot(data.frame(Max_component_size = noncomp2comp_distance1_ratio),aes(x=Max_component_size))+
  geom_histogram(fill="steelblue",bins=200)+
  theme_bw(base_size = 20)+
  xlab("Pecentage of isolated genes in subnetwork")+ylab("Count")+
  #xlim(0,45)+
  annotate("segment",x=distance1_ratio,xend = distance1_ratio,yend=0,y=6,size=2, arrow=arrow(),color="red")+
  annotate("text", x = distance1_ratio, y = 5, label = "IGON")


### test whether non-connected or iso_IGs can link to the IG componnet by permissive imprinted genes
intermediate_nodes <-c()
iso_IG_can_reach_component_via_1_step<-c()
for (node in iso_IGs){
  temp_shortest_path <- get.shortest.paths(graph = g,from = node,to = componet_IGs,weights = NA)
  print(listLen(temp_shortest_path$vpath))
  if (min(listLen(temp_shortest_path$vpath))>3){message("Shortest path with more than 3 nodes, not consider！")}
  if (min(listLen(temp_shortest_path$vpath))==3){
    for (path in temp_shortest_path$vpath){
      if (length(path)==3){
      #print (path)
        iso_IG_can_reach_component_via_1_step<-unique(c(iso_IG_can_reach_component_via_1_step,node))
      intermediate_nodes<-c(intermediate_nodes,names(path[2]))
      }
    }
  }
}
intermediate_nodes <-unique(intermediate_nodes) # 90 permissive

## intersect with permissive IGs
permissive_imprinted_data <-Arabidospis_imprinted_genes_parse("myresearch/network/data/6imprinted.list") #528
permissive_imp_not_in_strict <- permissive_imprinted_data$IG_in_networks[!permissive_imprinted_data$IG_in_networks%in% imprinted_genes_in_network] # 470

PINIS_intermediate<-length(intersect(intermediate_nodes,permissive_imp_not_in_strict))

random_ratio <-c()
for (m in (1:1000)){
  random_genes <-sample(names(V(g))[!names(V(g))%in%imprinted_genes_in_network],length(permissive_imp_not_in_strict))
  random_ratio<-c(random_ratio,length(intersect(intermediate_nodes,random_genes)))
}

ggplot(data.frame(Number = random_ratio),aes(x=random_ratio))+
  geom_histogram(fill="steelblue",bins=30,stat = "count")+
  theme_bw(base_size = 20)+
  xlab("Number of intermediated genes")+ylab("Frequency")+
  #xlim(0,45)+
  annotate("segment",x=PINIS_intermediate,xend=PINIS_intermediate,yend=0,y=30,size=2, arrow=arrow(),color="red")+
  annotate("text", x =PINIS_intermediate, y = 35, label = "Obseved")+
  scale_x_discrete(limits=seq(0,10,2))





### 7.1.6 connectivity
connectivity<-function(igraph_object){
  return(ecount(igraph_object)*2/(vcount(igraph_object)*vcount(igraph_object)-1))
}

AGFN_connectivity<-connectivity(g) #0.0021

IGPN<-induced_subgraph(g, unique(names(unlist(ego(g,order = 1,imprinted_genes_in_network))))) 
IGPN_connnectivity<-connectivity(IGPN) # 0.018

IGON<-induced_subgraph(g, imprinted_genes_in_network)
IGON_connectivity<-connectivity(IGON)  #0.02 / 0.0082

IGPN_sub_imp<-induced_subgraph(g,unique(names(unlist(ego(g,order = 1,imprinted_genes_in_network))))[!unique(names(unlist(ego(g,order = 1,imprinted_genes_in_network))))%in%imprinted_genes_in_network])
IGPN_sub_imp_conectivity <-connectivity(IGPN_sub_imp) #0.018 / 0.0085


random_connectivity <-c()
for (m in 1:1000){
  message(m)
  sample_g<-induced_subgraph(graph = g,vids = sample(names(V(g)),length(imprinted_genes_in_network)))
  random_connectivity<-c(random_connectivity,connectivity(sample_g))
}

library(ggplot2)
ggplot(data.frame(Connectivity=random_connectivity),
       aes(x=Connectivity))+
       geom_histogram(bins=200,fill="steelblue")+
       theme_bw()+
       theme(text = element_text(size=20))+ylab("Frequency")+
       annotate("segment",x=IGPN_connnectivity,xend = IGPN_connnectivity,yend=0,y=10,size=2, arrow=arrow(),color="steelblue")+
       annotate("segment",x=IGON_connectivity,xend = IGON_connectivity,  yend=0,y=10,size=2, arrow=arrow(),color="orange")+
       annotate("text", x = IGPN_connnectivity, y = 12, label = "IGPN")+
       annotate("text", x = IGON_connectivity, y = 12, label = "IGON")

### 7.1.6 cluster coefficent (transitivity)
AGFN_transitivity<-transitivity(g)    # 0.24
IGPN_transitivity<-transitivity(IGPN) # 0.38  
IGON_transitivity<-transitivity(IGON) # 0.31  #
IGPN_sub_imp_transitivity<- transitivity(IGPN_sub_imp) #0.38
AGFN_permute_transitivity <-transitivity(gg)

AGFN_transitivity_local<-mean(transitivity(g,type = "local"),na.rm = TRUE)    # 0.25
IGPN_transitivity_local<-mean(transitivity(IGPN,type = "local"),na.rm = TRUE) # 0.40
IGON_transitivity_local<-mean(transitivity(IGON,type = "local"),na.rm = TRUE) # 0.35
IGPN_sub_imp_transitivity_local <-mean(transitivity(IGPN_sub_imp,type = "local"),na.rm = TRUE)
AGFN_permute_transitivity_local <-mean(transitivity(gg,type = "local"),na.rm = TRUE)

permute_g_transitivity<- transitivity(permute_g) #0.24

ggtexttable(rbind(data.frame(Data="AGFN",Clust.Coef=AGFN_transitivity, mode = "global"),
                  data.frame(Data="IGPN",Clust.Coef=IGPN_transitivity, mode = "global") ,
                  data.frame(Data="IGON",Clust.Coef=IGON_transitivity, mode = "global"),
                  data.frame(Data="IGPN-imprinted",Clust.Coef = IGPN_sub_imp_transitivity, mode = "global"),
                  data.frame(Data="AGFN_permute",Clust.Coef = AGFN_permute_transitivity, mode = "global"),
                  
                  data.frame(Data="IGFN",Clust.Coef=AGFN_transitivity_local, mode = "local"),
                  data.frame(Data="IGPN",Clust.Coef=IGPN_transitivity_local, mode = "local"),
                  data.frame(Data="IGON",Clust.Coef=IGON_transitivity_local, mode = "local"),
                  data.frame(Data="IGPN-imprinted",Clust.Coef = IGPN_sub_imp_transitivity_local, mode = "local"),
                  data.frame(Data="AGFN_permute",Clust.Coef = AGFN_permute_transitivity_local, mode = "local")),
            theme = ttheme("mOrange",base_size = 15))


#### partner gene ontolgy enrichemnts (topgo)
partners<- names(V(IGPN))[!names(V(IGPN))%in% imprinted_genes]
source("/home/wuzefeng/R/topGO.R")
aa<-GO_enrichemnts(gene2go_file = "/home/wuzefeng/MyResearch/Imprinting_prediction/imprint_gene_list/TAIR10/imp_go_enrich/1gene2go.out",interested_genes = partners)
aa<-aa[apply(aa,1,function(x) x[7]<0.001),]  
write.csv(aa,file = "partners.go.csv",quote = FALSE,col.names = TRUE,row.names = FALSE)

### plant GSEA using enricher by clusterProfiler for non-imprited partener and all partner
library(clusterProfiler)
tf_target <- read.table("myresearch/network/data/1reformat/Ara_TFT.txt_out.txt",sep="\t",stringsAsFactors = FALSE)
Gfam <- read.table("myresearch/network/data/1reformat/Ara_GFam.txt_out.txt",sep="\t",stringsAsFactors = FALSE)
Kegg <- read.table("myresearch/network/data/1reformat/Ara_KEGG.txt_out.txt",sep="\t",stringsAsFactors = FALSE)
Lit <-read.table("myresearch/network/data/1reformat/Ara_LIT.txt_out.txt",sep="\t",stringsAsFactors = FALSE,quote = "")
MiR <- read.table("myresearch/network/data/1reformat/Ara_MIR.txt_out.txt",sep="\t",stringsAsFactors = FALSE,quote = "'")
PO <- read.table("myresearch/network/data/1reformat/Ara_PO.txt_out.txt",sep="\t",stringsAsFactors = FALSE)
Cyc <- read.table("myresearch/network/data/1reformat/Ara_Cyc.txt_out.txt",sep="\t",stringsAsFactors = FALSE)

length(unique(tf_target$V1))+length(unique(Gfam$V1))+length(unique(Kegg$V1))+length(unique(MiR$V1))+length(unique(PO$V1))+length(unique(Cyc$V1))

Term_size_filter<-function(df,size=5){
  df<-df%>%group_by(V1)%>%filter(n()>size)
}

Term_list = list(TF_target = tf_target,
                 Gene_family = Gfam,
                 KEGG_pathway = Kegg,
                 #Liture = Lit,
                 Micro_RNA = MiR,
                 Plant_ontology = PO,
                 Mata_pathway=Cyc)

for (name in names(Term_list)){
    message(name, "-Orignal term number: ",length(unique(Term_list[[name]]$V1)))
    temp_type <- Term_list[[name]]
    for (term_size in c(5,10,15)){
      type <- Term_size_filter(df = temp_type,size = term_size)
      message("...After filtering the term number: ",length(unique(type$V1)))
      enrich<-enricher(gene = names(V(IGPN_sub_imp)),
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    universe = names(V(g)),
                    qvalueCutoff = 0.2,
                    TERM2GENE = type[,c(1,3)],
                    TERM2NAME = unique(type[,c(1,2)]))
      enrich_df<-enrich@result[enrich@result$qvalue<0.2&enrich@result$p.adjust<0.05,] # filter
      enrich_df<-enrich_df%>%arrange(p.adjust) # sort by adjusted p-value
      
      message("...Enrich term number:",nrow(enrich_df))
      
      ### filter by term number (optional)
      if (nrow(enrich_df)>=10){
        enrich_df<-enrich_df[1:10,]
        }
      if (nrow(enrich_df)>=5&nrow(enrich_df)<10){
        enrich_df<- enrich_df[1:5,]     
        }
      
      ### visulization 
      
      if (nrow(enrich_df)>0){
        out_df<- data.frame(term_class = name, min_term_size = term_size, Description = enrich_df$Description, p.adjust = enrich_df$p.adjust)
        #write.table(out_df,file = "myresearch/network/result_data/PlantGSEA_enrich.txt",sep="\t", quote = FALSE,append = TRUE,col.names = !file.exists("myresearch/network/result_data/PlantGSEA_enrich.txt"),row.names = FALSE)
        write.table(out_df,file = "myresearch/network/result_data/PlantGSEA_permissive_enrich.txt",sep="\t", quote = FALSE,append = TRUE,col.names = !file.exists("myresearch/network/result_data/PlantGSEA_permissive_enrich.txt"),row.names = FALSE)
        
      #print(dotplot(enrich,showCategory=Enrichmed_term_number))
      }
      else{message("....NO enriched terms!")}
}
}
#dd<-read.table("myresearch/network/result_data/PlantGSEA_enrich.txt",sep="\t",header = TRUE)
dd<-read.table("myresearch/network/result_data/PlantGSEA_permissive_enrich.txt",sep="\t",header = TRUE)
ggplot(data = dd,aes(as.factor(min_term_size),Description,fill=-log(p.adjust),10))+
  geom_tile(color="black")+
  facet_grid(term_class~.,scales="free_y",space="free")+
  scale_fill_gradient(low = "red",high = "green")+
  theme_bw()+
  xlab("Minimum term size")

##### plant GSEA using enricher by clusterProfiler for all imprited partener

for (name in names(Term_list)){
  message(name, "-Orignal term number: ",length(unique(Term_list[[name]]$V1)))
  temp_type <- Term_list[[name]]
  for (term_size in c(5,10,15)){
    type <- Term_size_filter(df = temp_type,size = term_size)
    message("...After filtering the term number: ",length(unique(type$V1)))
    enrich<-enricher(gene = names(V(IGPN)),
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     universe = names(V(g)),
                     qvalueCutoff = 0.2,
                     TERM2GENE = type[,c(1,3)],
                     TERM2NAME = unique(type[,c(1,2)]))
    enrich_df<-enrich@result[enrich@result$qvalue<0.2&enrich@result$p.adjust<0.05,] # filter
    enrich_df<-enrich_df%>%arrange(p.adjust) # sort by adjusted p-value
    
    ### filter by term number (optional)
    if (nrow(enrich_df)>=10){
      enrich_df<-enrich_df[1:10,]
    }
    if (nrow(enrich_df)>=5&nrow(enrich_df)<10){
      enrich_df<- enrich_df[1:5,]     
    }
    
    ### visulization 
    
    if (nrow(enrich_df)>0){
      out_df<- data.frame(term_class = name, min_term_size = term_size, Description = enrich_df$Description, p.adjust = enrich_df$p.adjust)
      write.table(out_df,file = "myresearch/network/result_data/PlantGSEA_enrich_add_imp.txt",sep="\t", quote = FALSE,append = TRUE,
                  col.names = !file.exists("myresearch/network/result_data/PlantGSEA_enrich_add_imp.txt"),row.names = FALSE)
      message("...Enrich term number:",Enrichmed_term_number)
      #print(dotplot(enrich,showCategory=Enrichmed_term_number))
    }
    else{message("....NO enriched terms!")}
  }
}
dd<-read.table("myresearch/network/result_data/PlantGSEA_enrich_add_imp.txt",sep="\t",header = TRUE)
ggplot(data = dd,aes(as.factor(min_term_size),Description,fill=-log(p.adjust),10))+
  geom_tile(color="black")+
  facet_grid(term_class~.,scales="free_y",space="free")+
  scale_fill_gradient(low = "red",high = "green")+
  theme_bw()+
  xlab("Minimum term size")

################ merge enrich out with or without imprinted genes

dd1<-read.table("myresearch/network/result_data/PlantGSEA_enrich.txt",sep="\t",header = TRUE)
dd1$type<-"Partners Without IGs"

dd2<-read.table("myresearch/network/result_data/PlantGSEA_enrich_add_imp.txt",sep="\t",header = TRUE)
dd2$type<-"Partners With IGs"

dd<-rbind(dd1,dd2)
ggplot(data = dd,aes(as.factor(min_term_size),Description,fill=-log(p.adjust),10))+
  geom_tile(color="black")+
  facet_grid(term_class~type,scales="free_y",space="free")+
  scale_fill_gradient(low = "red",high = "green")+
  theme_bw()+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"))+
  xlab("Minimum term size")


## GO enrichment for non imprinted partners
## (GO slim)
library(clusterProfiler)
GO_slim_BP <- read.table("myresearch/network/data/1reform_GO/GO_slim/biological.go",sep="\t",stringsAsFactors = FALSE)
GO_slim_MF <- read.table("myresearch/network/data/1reform_GO/GO_slim/molecular_function.go",sep="\t",stringsAsFactors = FALSE)
GO_slim_CC <- read.table("myresearch/network/data/1reform_GO/GO_slim/cell_component.GO",sep="\t",stringsAsFactors = FALSE)


Term_list = list(BP = GO_slim_BP,MF = GO_slim_MF,CC = GO_slim_CC)

for (name in names(Term_list)){
  message(name, "-Orignal term number: ",length(unique(Term_list[[name]]$V1)))
  temp_type <- Term_list[[name]]
  for (term_size in c(5,10,15)){
    type <- Term_size_filter(df = temp_type,size = term_size)
    message("...After filtering the term number: ",length(unique(type$V1)))
    enrich<-enricher(gene = names(V(IGPN_sub_imp)),
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     universe = names(V(g)),
                     qvalueCutoff = 0.2,
                     TERM2GENE = type[,c(1,3)],
                     TERM2NAME = unique(type[,c(1,2)]))
    enrich_df<-enrich@result[enrich@result$qvalue<0.2&enrich@result$p.adjust<0.05,] # filter
    enrich_df<-enrich_df%>%arrange(p.adjust) # sort by adjusted p-value
    
    message("...Enrich term number:",nrow(enrich_df))
    
    ### filter by term number (optional)
    if (nrow(enrich_df)>=10){
      enrich_df<-enrich_df[1:10,]
    }
    if (nrow(enrich_df)>=5&nrow(enrich_df)<10){
      enrich_df<- enrich_df[1:5,]     
    }
    
    ### visulization 
    
    if (nrow(enrich_df)>0){
      out_df<- data.frame(term_class = name, min_term_size = term_size, Description = enrich_df$Description, p.adjust = enrich_df$p.adjust)
      #write.table(out_df,file = "myresearch/network/result_data/PlantGSEA_enrich.txt",sep="\t", quote = FALSE,append = TRUE,col.names = !file.exists("myresearch/network/result_data/PlantGSEA_enrich.txt"),row.names = FALSE)
      write.table(out_df,file = "myresearch/network/result_data/GO/GO_slim_enrich.txt",sep="\t", quote = FALSE,append = TRUE,
                  col.names = !file.exists("myresearch/network/result_data/GO/GO_slim_enrich.txt"),row.names = FALSE)
      
      #print(dotplot(enrich,showCategory=Enrichmed_term_number))
    }
    else{message("....NO enriched terms!")}
  }
} #non-imprinted partners
for (name in names(Term_list)){
  message(name, "-Orignal term number: ",length(unique(Term_list[[name]]$V1)))
  temp_type <- Term_list[[name]]
  for (term_size in c(5,10,15)){
    type <- Term_size_filter(df = temp_type,size = term_size)
    message("...After filtering the term number: ",length(unique(type$V1)))
    enrich<-enricher(gene = names(V(IGPN)),
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     universe = names(V(g)),
                     qvalueCutoff = 0.2,
                     TERM2GENE = type[,c(1,3)],
                     TERM2NAME = unique(type[,c(1,2)]))
    enrich_df<-enrich@result[enrich@result$qvalue<0.2&enrich@result$p.adjust<0.05,] # filter
    enrich_df<-enrich_df%>%arrange(p.adjust) # sort by adjusted p-value
    
    message("...Enrich term number:",nrow(enrich_df))
    
    ### filter by term number (optional)
    if (nrow(enrich_df)>=10){
      enrich_df<-enrich_df[1:10,]
    }
    if (nrow(enrich_df)>=5&nrow(enrich_df)<10){
      enrich_df<- enrich_df[1:5,]     
    }
    
    ### visulization 
    
    if (nrow(enrich_df)>0){
      out_df<- data.frame(term_class = name, min_term_size = term_size, Description = enrich_df$Description, p.adjust = enrich_df$p.adjust)
      #write.table(out_df,file = "myresearch/network/result_data/PlantGSEA_enrich.txt",sep="\t", quote = FALSE,append = TRUE,col.names = !file.exists("myresearch/network/result_data/PlantGSEA_enrich.txt"),row.names = FALSE)
      write.table(out_df,file = "myresearch/network/result_data/GO/GO_slim_all_enrich.txt",sep="\t", quote = FALSE,append = TRUE,
                  col.names = !file.exists("myresearch/network/result_data/GO/GO_slim_all_enrich.txt"),row.names = FALSE)
      
      #print(dotplot(enrich,showCategory=Enrichmed_term_number))
    }
    else{message("....NO enriched terms!")}
  }
} # all patners

## (plantgsea)

library(clusterProfiler)
GO_plantgsea_BP <- read.table("myresearch/network/data/1reform_GO/PlantGSEA/BP.GO",sep="\t",stringsAsFactors = FALSE)
GO_plantgsea_MF <- read.table("myresearch/network/data/1reform_GO/PlantGSEA/MF.GO",sep="\t",stringsAsFactors = FALSE)
GO_plantgsea_CC <- read.table("myresearch/network/data/1reform_GO/PlantGSEA/CC.GO",sep="\t",stringsAsFactors = FALSE)


Term_list = list(BP = GO_plantgsea_BP,MF = GO_plantgsea_MF,CC = GO_plantgsea_CC)

for (name in names(Term_list)){
  message(name, "-Orignal term number: ",length(unique(Term_list[[name]]$V1)))
  temp_type <- Term_list[[name]]
  for (term_size in c(5,10,15)){
    type <- Term_size_filter(df = temp_type,size = term_size)
    message("...After filtering the term number: ",length(unique(type$V1)))
    enrich<-enricher(gene = names(V(IGPN_sub_imp)),
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     universe = names(V(g)),
                     qvalueCutoff = 0.2,
                     TERM2GENE = type[,c(1,3)],
                     TERM2NAME = unique(type[,c(1,2)]))
    enrich_df<-enrich@result[enrich@result$qvalue<0.2&enrich@result$p.adjust<0.05,] # filter
    enrich_df<-enrich_df%>%arrange(p.adjust) # sort by adjusted p-value
    
    message("...Enrich term number:",nrow(enrich_df))
    
    ### filter by term number (optional)
    if (nrow(enrich_df)>=10){
      enrich_df<-enrich_df[1:10,]
    }
    if (nrow(enrich_df)>=5&nrow(enrich_df)<10){
      enrich_df<- enrich_df[1:5,]     
    }
    
    ### visulization 
    
    if (nrow(enrich_df)>0){
      out_df<- data.frame(term_class = name, min_term_size = term_size, Description = enrich_df$Description, p.adjust = enrich_df$p.adjust)
      #write.table(out_df,file = "myresearch/network/result_data/PlantGSEA_enrich.txt",sep="\t", quote = FALSE,append = TRUE,col.names = !file.exists("myresearch/network/result_data/PlantGSEA_enrich.txt"),row.names = FALSE)
      write.table(out_df,file = "myresearch/network/result_data/GO/GO_plantgsea_enrich.txt",sep="\t", quote = FALSE,append = TRUE,
                  col.names = !file.exists("myresearch/network/result_data/GO/GO_plantgsea_enrich.txt"),row.names = FALSE)
      
      #print(dotplot(enrich,showCategory=Enrichmed_term_number))
    }
    else{message("....NO enriched terms!")}
  }
}  # for mon-impritned partners
for (name in names(Term_list)){
  message(name, "-Orignal term number: ",length(unique(Term_list[[name]]$V1)))
  temp_type <- Term_list[[name]]
  for (term_size in c(5,10,15)){
    type <- Term_size_filter(df = temp_type,size = term_size)
    message("...After filtering the term number: ",length(unique(type$V1)))
    enrich<-enricher(gene = names(V(IGPN)),
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     universe = names(V(g)),
                     qvalueCutoff = 0.2,
                     TERM2GENE = type[,c(1,3)],
                     TERM2NAME = unique(type[,c(1,2)]))
    enrich_df<-enrich@result[enrich@result$qvalue<0.2&enrich@result$p.adjust<0.05,] # filter
    enrich_df<-enrich_df%>%arrange(p.adjust) # sort by adjusted p-value
    
    message("...Enrich term number:",nrow(enrich_df))
    
    ### filter by term number (optional)
    if (nrow(enrich_df)>=10){
      enrich_df<-enrich_df[1:10,]
    }
    if (nrow(enrich_df)>=5&nrow(enrich_df)<10){
      enrich_df<- enrich_df[1:5,]     
    }
    
    ### visulization 
    
    if (nrow(enrich_df)>0){
      out_df<- data.frame(term_class = name, min_term_size = term_size, Description = enrich_df$Description, p.adjust = enrich_df$p.adjust)
      #write.table(out_df,file = "myresearch/network/result_data/PlantGSEA_enrich.txt",sep="\t", quote = FALSE,append = TRUE,col.names = !file.exists("myresearch/network/result_data/PlantGSEA_enrich.txt"),row.names = FALSE)
      write.table(out_df,file = "myresearch/network/result_data/GO/GO_plantgsea_all_enrich.txt",sep="\t", quote = FALSE,append = TRUE,
                  col.names = !file.exists("myresearch/network/result_data/GO/GO_plantgsea_all_enrich.txt"),row.names = FALSE)
      
      #print(dotplot(enrich,showCategory=Enrichmed_term_number))
    }
    else{message("....NO enriched terms!")}
  }
} # for all partners


### plot
dd1<-read.table("myresearch/network/result_data/GO/GO_plantgsea_enrich.txt",sep="\t",header = TRUE)
dd1$type<-"Partners Without IGs"

dd2<-read.table("myresearch/network/result_data/GO/GO_plantgsea_all_enrich.txt",sep="\t",header = TRUE)
dd2$type<-"Partners With IGs"

dd<-rbind(dd1,dd2)
ggplot(data = dd,aes(as.factor(min_term_size),Description,fill=-log(p.adjust),10))+
  geom_tile(color="black")+
  facet_grid(term_class~type,scales="free_y",space="free")+
  scale_fill_gradient(low = "red",high = "green")+
  theme_bw()+
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"))+
  xlab("Minimum term size")














###7.3  imprinted genes degree 
maternal_degree <-degree(g,maternal_imprint)
paternal_degree <-degree(g,paternal_imprint)
imprint_degree <-degree(g,imprinted_genes_in_network)
partner_degree<-degree(g,partners)
boxplot(list(as.numeric(degree(g,imprinted_genes_in_network)),as.numeric(maternal_degree),as.numeric(paternal_degree),as.numeric(degree(g))[as.numeric(degree(g))<=200]),
            names = c("Imprinted","Maternal","Paternal","All"),
            col=c("purple","red","blue","green"),ylab='Degree',
            cex.lab=2,cex.axis=1.5)

library(ggplot2)

df1<-data.frame(Degree=maternal_degree,class="MEG")
df2<-data.frame(Degree=paternal_degree,class="PEG")
df3<-data.frame(Degree=imprint_degree,class='IGs')
df4<-data.frame(Degree=partner_degree,class="Partners")
df5<-data.frame(Degree=as.numeric(degree(g,v = names(V(g))[!names(V(g))%in%imprinted_genes_in_network])),class="Non-IGs")

df<-rbind(df1,df2,df3,df4,df5)

ggplot(df,aes(y=log(Degree,2),x=class,col=class))+
  geom_violin(trim = FALSE,aes(fill=class),alpha=0.5,col="white")+
  theme_bw(base_size = 20)+
  geom_boxplot(width=0.2)+
  geom_signif(comparisons = list(c("MEG","Non-IGs"),c("PEG","Non-IGs"),c("IGs","Non-IGs"),c("Partners","Non-IGs")),
              test = "wilcox.test",map_signif_level = "TRUE",
              test.args = list(alternative = "greater"),
              step_increase = 0.06,
              tip_length = 0.01)+
  #scale_fill_manual(values = c("tomato","steelblue", "orange")) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
  ylab("Degree (log2)")


#### 7.4.1 imprinted genes betweeness 
IGFN_betweenness<-read.table("myresearch/network/data/0.1AFGN_betweenness_normalized.txt",header = TRUE,stringsAsFactors = FALSE)
imprinted_betweenness <-IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%imprinted_genes_in_network]
maternal_betweenness <-IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%maternal_imprint]
paternal_betweenness <-IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%paternal_imprint]
IGPN_betweenness <- IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%names(V(IGPN))]
IGPN_betweenness_remov_imp <-IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%names(V(IGPN))& !IGFN_betweenness$genes%in% imprinted_genes_in_network]

df1<-data.frame(Betweenness = log(IGFN_betweenness$betweenness[!IGFN_betweenness$genes%in%imprinted_genes_in_network],2), class="Non_IGs")
df2<-data.frame(Betweenness = log(imprinted_betweenness,2),class="IGs")
df3<-data.frame(Betweenness = log(maternal_betweenness,2), class="MEG")
df4<-data.frame(Betweenness = log(paternal_betweenness,2), class="PEG")
df5<-data.frame(Betweenness = log(IGPN_betweenness,2), class = "IGPN")
df6<-data.frame(Betweenness = log(IGPN_betweenness_remov_imp,2),class="Non_IG_IGPN")

df<-rbind(df1,df2,df3,df4,df5,df6)
df$class<-factor(x = df$class,levels = c("IGs","MEG","PEG","Non_IGs","IGPN","Non_IG_IGPN"))
ggplot(df,aes(y=Betweenness,col=class,x=class))+
        geom_violin(trim = FALSE,aes(fill=class),alpha=0.5,col="white")+
        geom_boxplot(width=0.2)+
        theme_bw()+
        theme(text = element_text(size = 20),legend.position = "none")+
        ylab("Betweenness (log2)")+
        xlab("")+
        geom_signif(comparisons = list(c("IGs","Non_IGs"),c("MEG","Non_IGs"),c("PEG","Non_IGs"),c("IGPN","IGs")),
                    map_signif_level = TRUE,
                    step_increase = 0.08,
                    tip_length = 0.01,
                    test.args = list(alternative = "greater"))

### 7.4.5 imprinted genes closenness

IGFN_closeness<-read.table("myresearch/network/data/0.2closeness_normlized.txt",header = TRUE,stringsAsFactors = FALSE)
imprinted_closeness <-IGFN_closeness$closesness[IGFN_closeness$gene_name%in%imprinted_genes_in_network]
maternal_closeness <-IGFN_closeness$closesness[IGFN_closeness$gene_name%in%maternal_imprint]
paternal_closeness <-IGFN_closeness$closesness[IGFN_closeness$gene_name%in%paternal_imprint]
IGPN_closeness <- IGFN_closeness$closesness[IGFN_closeness$gene_name%in%names(V(IGPN))]
IGPN_closeness_remov_imp <-IGFN_closeness$closesness[IGFN_closeness$gene_name%in%names(V(IGPN))& !IGFN_closeness$gene_name%in% imprinted_genes_in_network]


df1<-data.frame(Closeness=IGFN_closeness$closesness[!IGFN_closeness$gene_name%in%imprinted_genes_in_network], class="Non_IGs")
df2<-data.frame(Closeness=imprinted_closeness,class="IGs")
df3<-data.frame(Closeness=maternal_closeness, class="MEG")
df4<-data.frame(Closeness=paternal_closeness, class="PEG")
df5<-data.frame(Closeness=IGPN_closeness, class="IGPN")
df6<-data.frame(Closeness=IGPN_closeness_remov_imp,class="Non_IGs_IGPN")

df<-rbind(df1,df2,df3,df4,df5,df6)
df<-df[df$Closeness>0.001,]
df$class<-factor(x = df$class,levels = c("IGs","MEG","PEG","Non_IGs","IGPN","Non_IGs_IGPN"))
ggplot(df,aes(y=Closeness,col=class,x=class))+
  geom_violin(trim = FALSE,aes(fill=class),alpha=0.5,col="white")+
  geom_boxplot(width=0.2)+
  theme_bw()+
  theme(text = element_text(size = 20),legend.position = "none")+
  ylab("Closeness")+
  xlab("")+
  geom_signif(comparisons = list(c("IGs","Non_IGs"),c("MEG","Non_IGs"),c("PEG","Non_IGs"),c("IGPN","IGs")),
                map_signif_level = FALSE,
                step_increase = 0.08,
                tip_length = 0.01,
                test.args = list(alternative = "greater"))
              
### add AGPN
dff<-rbind(df5,df2,df1)
dff<-dff[dff$Closeness>0.001,]
dff$class<-factor(x = dff$class,levels = c("IGPN","IGs","All"))
ggplot(dff,aes(y=Closeness,fill=class,x=class))+geom_boxplot()+
  theme_bw()+
  theme(text = element_text(size = 20),legend.position = "none")+
  ylab("Closeness")+
  xlab("")+
  geom_signif(comparisons = list(c("IGPN","All"),c("IGPN","IGs"),c("IGs","All")),
              map_signif_level = FALSE,
              step_increase = 0.08,tip_length = 0.01)

### remove imprinted genes from IGPN
dff<-rbind(df6,df2,df1)
dff<-dff[dff$Closeness>0.001,]
dff$class<-factor(x = dff$class,levels = c("IGPN","IGs","All"))
ggplot(dff,aes(y=Closeness,fill=class,x=class))+geom_boxplot()+
  theme_bw()+
  theme(text = element_text(size = 20),legend.position = "none")+
  ylab("Closeness")+
  xlab("")+
  geom_signif(comparisons = list(c("IGPN","All"),c("IGPN","IGs")),
              map_signif_level = FALSE,
              step_increase = 0.08,tip_length = 0.01)


### 7.5 imporinted genes eigenvector
eigns<-eigen_centrality(graph = g,directed = FALSE,weights = NA)
maternal_eigns<-eigns$vector[maternal_imprint_in_network]
paternal_eigns<-eigns$vector[paternal_imprint_in_network]
boxplot(list(log(maternal_eigns,2),log(paternal_eigns,2)),col=c("red","steelblue"),
        names = c("Maternal","Paternal"),ylab='Eigenvector (log2)',
        cex.lab=2,cex.axis=2)

######### overlap% impritned genes with hub genes (seting degree) (optional )
IGPN_sub_imp_ratio<-c()
imprinted_genes_imp_ratio <-c()
random_genes_ratio<-c()

random_genes<-sample(x = names(V(g)),size = length(V(IGON)),replace = FALSE)

for (r in seq(0,800,10)){
  hub_genes<- names(V(g))[degree(g)>=r]
  IGPN_sub_imp_ratio<-c(IGPN_sub_imp_ratio,(length(intersect(names(V(IGPN_sub_imp)),hub_genes))/vcount(IGPN_sub_imp)))
  imprinted_genes_imp_ratio<-c(imprinted_genes_imp_ratio,length(intersect(imprinted_genes_in_network,hub_genes))/length(imprinted_genes_in_network))
  random_genes_ratio<-c(random_genes_ratio,length(intersect(random_genes,hub_genes))/length(random_genes))
}


overlap_hub_ratio<-data.frame(threshold=rep(seq(0,800,10),3),ratio=c(IGPN_sub_imp_ratio,imprinted_genes_imp_ratio,random_genes_ratio),class=c(rep("IGPs",81),rep("IGs",81),rep("Random",81)))
ggplot(overlap_hub_ratio,aes(x=threshold,y=ratio,color=class))+
      geom_line(size=1.5)+
      theme_bw(base_size = 20)+
      theme(legend.position = c(0.85,0.85))+
      ylab("Percentage of genes")+xlab("Hub threshold")+
      annotate("segment", x=20, xend=150, y=0.75, yend=0.80, color="black", size=1, arrow=arrow())+
      annotate("text", x =250, y = 0.80, label = "P-value = 2.8e-10")+#8e-15
      annotate("segment", x=20, xend=150, y=0.5, yend=0.55, color="black", size=1, arrow=arrow())+
      annotate("text", x =250, y = 0.55, label = "P-value = 0.094")#8e-15


#####  hub top 5%, top 10% and 15% (Figure 2B)
threshold = c(0.90,0.95,0.99)

hub_genes<- names(V(g))[degree(g)>=quantile(degree(g),probs = threshold[3])]
IG_belong_hubs<-intersect(imprinted_genes_in_network,hub_genes) #2
IG_non_hubs<-imprinted_genes_in_network[!imprinted_genes_in_network%in%IG_belong_hubs] # 80

IG_non_hubs2hub_genes<-shortest.paths(graph = g,v = IG_non_hubs,to = hub_genes,weights = NA)
IG_non_hubs2hub_min_distance<-apply(IG_non_hubs2hub_genes,1,min)
IG_non_hubs2hub_direct_linked<-sum(IG_non_hubs2hub_min_distance==1) # 55  (55/80)

### random sampling imprinted geens and count the number of genes directing interact with  hub genes   
fun<-function(x){
  random_genes <- sample(names(V(g)),82)
  random_belong_hubs <- intersect(random_genes,hub_genes)
  random_non_hubs <-random_genes[!random_genes%in%hub_genes]
  random_non_hubs2hub_genes<-shortest.paths(graph = g,v = random_non_hubs,to = hub_genes,weights = NA)
  random_non_hubs2hub_min_distance<-apply(random_non_hubs2hub_genes,1,min)
  random_non_hubs2hub_direct_linked<-sum(random_non_hubs2hub_min_distance==1) # 55  (55/80)
  return(random_non_hubs2hub_direct_linked/length(random_non_hubs))
}
system.time(random_results<-mclapply(1:1000,fun,mc.cores = 60))
random_results<-unlist(random_results)

ggplot(as.data.frame(random_results),aes(x=random_results))+geom_histogram(bins=200,fill="steelblue")+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("% genes interacting with hub genes")+
  annotate("segment", x=IG_non_hubs2hub_direct_linked/length(IG_non_hubs), xend=IG_non_hubs2hub_direct_linked/length(IG_non_hubs), y=5, yend=0, color="red", size=2, arrow=arrow())+
  annotate("text", x =IG_non_hubs2hub_direct_linked/length(IG_non_hubs), y = 6, label = "IGs",col="black")
  

### ROC or prediction of partner genes
score_list<-list()  # record each cor
label_list<-list()
modenames<- c("Partners","PEGs","MEGs")

score_list[[1]]<-c(degree(g,names(V(IGPN_sub_imp))),
                   degree(g,names(V(g))[!names(V(g))%in%names(V(IGPN_sub_imp))]))
label_list[[1]]<-c(rep(1,length(names(V(IGPN_sub_imp)))),
                   rep(0,length(names(V(g))[!names(V(g))%in%names(V(IGPN_sub_imp))])))

score_list[[2]]<-c(degree(g,paternal_imprint),
                   degree(g,names(V(g))[!names(V(g))%in%paternal_imprint]))
label_list[[2]]<-c(rep(1,length(paternal_imprint)),
                   rep(0,length(names(V(g))[!names(V(g))%in%paternal_imprint])))

score_list[[3]]<-c(degree(g,maternal_imprint),
                   degree(g,names(V(g))[!names(V(g))%in%maternal_imprint]))
label_list[[3]]<-c(rep(1,length(maternal_imprint)),
                   rep(0,length(names(V(g))[!names(V(g))%in%maternal_imprint])))
dsids <-seq(0,2)

mdata <- mmdata(score_list, labels = label_list,modnames = modenames, dsids = dsids)
mmcurves <- evalmod(mdata)
autoplot(mmcurves,curvetype = "ROC")

### 7.6 imprinted genes years
library(readxl)
gene_olds <-as.data.frame(read_excel(path = "other_resources/nature11394-s2_gene_time.xls",sheet = "first dataset",col_names = TRUE))
gene_olds$Gene<-toupper(gene_olds$Gene)
boxplot(gene_olds$Phylostratum,gene_olds$Phylostratum[gene_olds$Gene%in%imprinted_genes],
        gene_olds$Phylostratum[gene_olds$Gene%in%maternal_imprint],
        gene_olds$Phylostratum[gene_olds$Gene%in%paternal_imprint],
        col=c('#9ED2F0','#E6A429',"tomato","steelblue"),
        names=c("All genes","Imprinted genes","Maternal","Paternal"),
        ylab="Phylostratum",cex.lab=2,cex.axis=1.2) # no difference in gene year between imprinted and non-imprinted 

## class imprinted genes into two classed based on gene olds
#7.6.1 degree by old
old_imprinted_genes<-gene_olds$Gene[gene_olds$Phylostratum<=4&gene_olds$Gene%in%imprinted_genes_in_network]
younger_imprinted_genes<-gene_olds$Gene[gene_olds$Phylostratum>4&gene_olds$Gene%in%imprinted_genes_in_network]

old_paternal_genes <-gene_olds$Gene[gene_olds$Phylostratum<=4&gene_olds$Gene%in%paternal_imprint]
younger_paternal_genes <-gene_olds$Gene[gene_olds$Phylostratum>4&gene_olds$Gene%in%paternal_imprint]

old_maternal_genes <-gene_olds$Gene[gene_olds$Phylostratum<=4&gene_olds$Gene%in%maternal_imprint]
younger_maternal_genes <-gene_olds$Gene[gene_olds$Phylostratum>4&gene_olds$Gene%in%maternal_imprint]

par(mar=c(15,5,5,5))
boxplot(degree(g,old_imprinted_genes),
        degree(g,younger_imprinted_genes),
        degree(g,old_maternal_genes),
        degree(g,younger_maternal_genes),
        degree(g,old_paternal_genes),
        degree(g,younger_paternal_genes),
        
        col=c("#CDC1C5", "#8B8970", "#FF7F50", "#CD3333", "#87CEFA", "#009ACD"),
        names=c("Older imprinted genes","Younger imprinted genes","Older maternal genes","Younger maternal genes","Older paternal genes","Younger paternal genes"),
        ylab="Degree",cex.lab=2,cex.axis=1.2,las=2)

older_genes<-gene_olds$Gene[gene_olds$Phylostratum<=3&gene_olds$Gene%in%names(V(g))]
younger_genes<-gene_olds$Gene[gene_olds$Phylostratum>3&gene_olds$Gene%in%names(V(g))]
boxplot(degree(g,older_genes),degree(g,younger_genes),names=c("old","younger"))

## 7.6.2 shortest path by gene old

par(mar=c(15,5,5,5))
boxplot(as.numeric(distances(g,old_imprinted_genes,old_imprinted_genes,weights = NA)),
        as.numeric(distances(g,younger_imprinted_genes,younger_imprinted_genes,weights = NA)),
        as.numeric(distances(g,old_maternal_genes,old_maternal_genes,weights = NA)),
        as.numeric(distances(g,younger_maternal_genes,younger_maternal_genes,weights = NA)),
        as.numeric(distances(g,old_paternal_genes,old_paternal_genes,weights = NA)),
        as.numeric(distances(g,younger_paternal_genes,younger_paternal_genes,weights = NA)),
        
        col=c("#CDC1C5", "#8B8970", "#FF7F50", "#CD3333", "#87CEFA", "#009ACD"),
        names=c("Older imprinted genes","Younger imprinted genes","Older maternal genes","Younger maternal genes","Older paternal genes","Younger paternal genes"),
        ylab="Shortest path",cex.lab=2,cex.axis=1.2,las=2)

older_genes<-gene_olds$Gene[gene_olds$Phylostratum<=3&gene_olds$Gene%in%names(V(g))]
younger_genes<-gene_olds$Gene[gene_olds$Phylostratum>3&gene_olds$Gene%in%names(V(g))]
boxplot(degree(g,older_genes),degree(g,younger_genes),names=c("old","younger"))


## 7.7 imprinted genes pcc classification
edges_list2df<-as.data.frame(as_edgelist(g)) # some genes not in networks
gene_expression<-read.table("~/MyResearch/microarray/AT40_RMA_MAS/7all_gene_expresion.txt",header = TRUE)
edges_list2df1<-subset(edges_list2df,edges_list2df$V1%in%rownames(gene_expression)&edges_list2df$V2%in%rownames(gene_expression))

gg_pcc<-cor(t(gene_expression))
edges_list2df1$PCC<-apply(edges_list2df1,1,function(x)gg_pcc[x[1],x[2]])
rm(gg_pcc)&gc()
#edges_list2df<-merge(edges_list2df,edges_list2df1,by.x = c("V1","V2"),all.x=T)
#edges_list2df$PCC[is.na(edges_list2df$PCC)]<-0
E(g)$weight<-edges_list2df$PCC

## plot
sub_g<-induced_subgraph(graph = g,vids = imprinted_genes_in_network)
sub_g<-induced_subgraph(graph = sub_g,names(components(sub_g)$membership[components(sub_g)$membership==1]))
par(bg="white",mar=c(2,2,2,2))

V(sub_g)$color[names(V(sub_g))%in%maternal_imprint_in_network]<-"tomato" # vertex color 
V(sub_g)$color[names(V(sub_g))%in%paternal_imprint_in_network]<-"steelblue"

V(sub_g)$size <- 8              # vertex size
V(sub_g)$label.cex <- 0.8         # vertex label size
V(sub_g)$label.color<-"black"
E(sub_g)$color<-ifelse(E(sub_g)$weight>0.3,"orange","gray")        # edge color 

E(sub_g)$width=3
plot(sub_g,layout=layout.fruchterman.reingold,
     vertex.frame.color= NA,
     vertex.color=V(sub_g)$color,
     vetex.shape=V(sub_g)$shape)

legend('topleft',
       #bg="white",
       text.col="black",
       legend=c("Paternal","Maternal","Correlated","Anti-correlated"),
       pch=c(19,19,NA,NA), #shape
       lty = c(0.5,0.5,1,1),
       box.lty = 2, # 
       pt.cex = 3, #lines size 
       cex=1, #box size
       col=c("tomato","steelblue","orange","gray"),
       y.intersp=1.5)

#text(-1, -1,"bbbbbb",col = "white")


#########################################
wc<-cluster_edge_betweenness(sub_g,weights = NULL)
plot(wc,sub_g)

### imp distance to exp corrs
IGs_distance<-shortest.paths(graph = g,v = intersect(colnames(hh),imprinted_genes_in_network),intersect(colnames(hh),imprinted_genes_in_network),weights = NA)
IGs_distance<-data.frame(row=rownames(IGs_distance)[row(IGs_distance)[upper.tri(IGs_distance)]],
                         col=colnames(IGs_distance)[col(IGs_distance)[upper.tri(IGs_distance)]],
                         d=IGs_distance[upper.tri(IGs_distance)])
IGs_distance<-subset(IGs_distance,is.finite(d))
IGs_distance$d<-ifelse(IGs_distance$d>=5,5,IGs_distance$d)

IGs_distance$exp_cor<-apply(IGs_distance,1,function(x)cor(hh[,c(x[1],x[2])])[2,1])
ggplot(IGs_distance,aes(y=exp_cor,x=as.factor(d),fill=as.factor(d)))+
  geom_boxplot()+
  theme_bw()+
  theme(text=element_text(size = 20),legend.position = "none")+
  ylab("Gene expression correlation")+
  xlab("Shortest path")
  #geom_smooth(method = "lm", se=TRUE, color="orange", aes(group=1),alpha=0.1,linetype = "dashed")

## optional visuliazation  
ggplot(IGs_distance,aes(y=exp_cor,x=as.factor(d),col=as.factor(d)))+
  geom_violin(trim = FALSE,aes(fill=as.factor(d)),alpha=0.5,col="white")+
  theme_bw(base_size = 20)+
  geom_boxplot(width=0.2)+
  ylab("Gene expression correlation (PCC)")+
  xlab("Shortest path")+
  theme(legend.position = "none")
  
  
  
  
## paternal 2 maternal distancee
PEG2MEGs<-shortest.paths(g,v = paternal_imprint[paternal_imprint%in%rownames(gene_expression)],
                           maternal_imprint[maternal_imprint%in%rownames(gene_expression)],weights = NA)
PEG2MEGs<-as.data.frame(as.table(PEG2MEGs))
PEG2MEGs$pcc<-apply(PEG2MEGs,1,function(x)(gg_pcc[x[1],x[2]]))
ggplot(data = PEG2MEGs) +
  aes(x = Freq, fill = Freq,y=pcc) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

## paternal 2 maternal distence and comapared with to parters
paternal_in_exp <-paternal_imprint[paternal_imprint%in%rownames(gene_expression)]
maternal_in_exp <-maternal_imprint[maternal_imprint%in%rownames(gene_expression)]

paternal_partners<-unique(unlist(neighborhood(graph = g,nodes = paternal_in_exp))) # paternal neighours
paternal_partners<-paternal_partners[!paternal_partners%in%paternal_imprint]

maternal_partners<-unique(unlist(neighborhood(graph = g,nodes = maternal_in_exp)))
maternal_partners<-maternal_partners[!maternal_partners%in%maternal_imprint]

PEG2partners<-shortest_paths(graph = g,from = paternal_imprint[paternal_imprint%in%rownames(gene_expression)])


### add ds value from lyrata  
ds <- read.table("other_resources/DS_to_lyrata/TAIR2lyrata.ds.mart_export.txt",header = TRUE,sep="\t")
ds <- na.omit(ds)
ds <- as.data.frame(ds%>%group_by(Gene.stable.ID)%>%summarise_at('dS.with.Arabidopsis.lyrata',max))

Delta_ds<-function(g1,g2){
  if(g1%in%ds$Gene.stable.ID&g2%in%ds$Gene.stable.ID){
    delta_ds<-abs(ds$dS.with.Arabidopsis.lyrata[ds$Gene.stable.ID==g1]-ds$dS.with.Arabidopsis.lyrata[ds$Gene.stable.ID==g2])
  }
  else{
    Delta_ds<-NA
  }
}

PEG2MEGs$ds_diff<-unlist(apply(PEG2MEGs,1,function(x)Delta_ds(x[1],x[2])))
PEG2MEGs<-na.omit(PEG2MEGs)
PEG2MEGs<-PEG2MEGs[PEG2MEGs$ds_diff<=1,]
ggplot(PEG2MEGs,aes(x=as.factor(Freq),y=ds_diff))+geom_boxplot(notch = FALSE)


##7.8 imprinted and transcription factors distance
##7.8 imprinted and transcription factors distance

ara_TFs<-read.table("/home/wuzefeng/MyResearch/motif_dbs/2plantTFDB/Arabidopsis_TF_list/Ath_TF_list",header = FALSE,stringsAsFactors = FALSE)
ara_TFs<-unique(ara_TFs$V2)
TFs_in_networks <- intersect(ara_TFs,names(V(g))) #1526
sub_g<-induced_subgraph(g,vids = c(imprinted_genes_in_network,TFs_in_networks))

plot_personal_graph<-function(graph_object){
  V(sub_g)$color[names(V(sub_g))%in%maternal_imprint]<-"tomato" # vertex color 
  V(sub_g)$color[names(V(sub_g))%in%paternal_imprint]<-"steelblue"
  V(sub_g)$size <- 8              # vertex size
  V(sub_g)$label.cex <- 0.8          # vertex label size
  V(sub_g)$label.color<-"black"
  E(sub_g)$color <- "gray"        # edge color 
  E(sub_g)$width=3
  
  plot.igraph(sub_g,layout=layout.fruchterman.reingold,vertex.frame.color= "white")
  legend('topleft',
         legend=c("Maternal","Paternal"),
         pch=19, #shape
         box.lty=2, # 
         pt.cex= 3, #lines size 
         cex=1, #box size
         col=c("tomato","steelblue"),
         y.intersp=1.5)
}
plot_personal_graph(graph_object = sub_g)

summary(as.numeric(distances(g,paternal_imprint,TFs_in_networks,weights = NA))) #no difference
summary(as.numeric(distances(g,maternal_imprint,TFs_in_networks,weights = NA)))
summary(as.numeric(distances(g,sample(names(V(g)),length(TFs_in_networks)),TFs_in_networks,weights = NA)))

##7.7.1 TF-imprinted martix heatmaps

imprinted_paternal_TF <-intersect(TFs_in_networks,paternal_imprint)
imprinted_maternal_TF <-intersect(TFs_in_networks,maternal_imprint)

m<-matrix(data = 0,nrow = length(TFs_in_networks),ncol = length(imprinted_genes_in_network)) # 1526 * 82/468
rownames(m)<-TFs_in_networks
colnames(m)<-imprinted_genes_in_network

for (imp in imprinted_genes_in_network){
  nbs<-names(unlist(neighborhood(g,order = 1,nodes = imp)))
  tf_nbs<-intersect(nbs,TFs_in_networks)
  m[tf_nbs,imp]<-1
}

mm<-m[,colSums(m)>0]
mmm<-mm[rowSums(mm)>0,] #226 * 42 | 501 * 229
## pheatmap
annotation <- data.frame(Pattern = factor(c(rep("MEG",length(maternal_imprint)),
                                          rep("PEG",length(paternal_imprint))), 
                                          labels = c("Maternal", "Paternal")))
rownames(annotation) <- c(maternal_imprint,paternal_imprint) # check out the row names of annotation
annotation<-subset(annotation,rownames(annotation)%in%colnames(mmm))
library(pheatmap)
pheatmap(t(mmm),legend_breaks = 0:1, legend_labels = c("0","1"),fontsize_col =5,annotation_row = annotation)


## TF-imprinted genes turn into networks module
tf_g<-as.data.frame(as.table(mmm),stringsAsFactors = FALSE)
tf_g<-tf_g[tf_g$Freq>0,]
tf_g<-tf_g[tf_g$Var1!=tf_g$Var2,]

tf_imp_net<-graph.data.frame(d = tf_g, directed = FALSE)
tf_imp_net<-simplify(tf_imp_net,remove.multiple = TRUE,remove.loops = TRUE)
sub_g<-tf_imp_net



### 7.9.2 integrating  motif enrichment data analysis (no use)

motif_gene_mappping<-read.table("/home/wuzefeng/MyResearch/motif_dbs/1jaspar2018/5gene_mapping_from2database")
mapping_function<-function(x){
  if (sum(is.na(x[c(3,4,5)]))==length(x[c(3,4,5)])) return("NA")
  else return(na.omit(unique(x[c(3,4,5)])))
}
motif_gene_mappping$gene_id<-apply(motif_gene_mappping,1,FUN = mapping_function)

motif_enrich_tf<-read.table("/home/wuzefeng/MyResearch/motif_dbs/1jaspar2018/1fimo_out/Tair/imp_vs_all.txt",header = TRUE,sep = "\t")
motif_enrich_tf$gene_id <-unlist(motif_gene_mappping$gene_id[match(rownames(motif_enrich_tf),motif_gene_mappping$V2)])
intersect(rownames(mmm),motif_enrich_tf$gene_id) # "AT1G32640" "AT2G03710" "AT5G03150" "AT5G44160" "AT4G32730"

### 7.9.3 integrating ppi data analysis for visulization

ara_reg<-read.table("../PPI/Arabidopsis/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-3.4.134.tab2.txt",sep="\t",header = TRUE,stringsAsFactors = FALSE)
ppi_imprinted<-ara_reg[ara_reg$Systematic.Name.Interactor.A%in%imprinted_genes_in_network | ara_reg$Systematic.Name.Interactor.A%in%imprinted_genes_in_network,]
ppi_imprinted<-data.frame(a=ppi_imprinted$Systematic.Name.Interactor.A,b=ppi_imprinted$Systematic.Name.Interactor.B,stringsAsFactors = FALSE)
ppi_imprinted<-ppi_imprinted[ppi_imprinted$a!=ppi_imprinted$b&ppi_imprinted$b!="-",]
ppi_imprinted$a_info<-ifelse(ppi_imprinted$a%in%imprinted_genes_in_network,"Y","N")
ppi_imprinted$b_info<-ifelse(ppi_imprinted$b%in%imprinted_genes_in_network,"Y","N")
ppi_imprinted<-ppi_imprinted[ppi_imprinted$a%in%names(V(tf_imp_net))&ppi_imprinted$b%in%names(V(tf_imp_net)),]
ppi_imprinted<-ppi_imprinted[apply(ppi_imprinted,1,function(x) are.connected(graph = tf_imp_net,v1 = x[1],v2 = x[2])),]

par(bg="white",mar=c(2,2,2,2))
plot_personal_graph<-function(graph_object){
  #V(sub_g)$color[names(V(sub_g))%in%focus_links$V2]<-"purple"
  V(sub_g)$color<-"green"
  V(sub_g)$frame.color<-"white"
  V(sub_g)$color[names(V(sub_g))%in%maternal_imprint]<-"tomato" # vertex color 
  V(sub_g)$color[names(V(sub_g))%in%paternal_imprint]<-"steelblue"
  V(sub_g)$color[names(V(sub_g))%in%imprinted_paternal_TF]<-"pink"
  V(sub_g)$color[names(V(sub_g))%in%imprinted_maternal_TF]<-"yellow"
  
  V(sub_g)$size <- 4              # vertex size
  V(sub_g)$label.cex <- 0.5         # vertex label size
  V(sub_g)$label.color<-"black"
  
  E(sub_g)$color <- "gray"        # edge color 
  E(sub_g)$width=2#E(sub_g)$weight/max(E(sub_g)$weight)
  E(sub_g,P=as.vector(t(as.matrix(ppi_imprinted[,c(1,2)]))))$color<-"purple"
  E(sub_g,P=as.vector(t(as.matrix(ppi_imprinted[,c(1,2)]))))$width<-3
  #E(sub_g)$width=edge.betweenness(sub_g)
  plot.igraph(sub_g,layout=layout.fruchterman.reingold)
  legend('topleft',
         #bg="white",
         text.col="tomato",
         legend=c("MEG","PEG","TF","MEG TFs","PEG TFs"),
         pch=19, #shape
         box.lty=2, # 
         pt.cex= 3, #lines size 
         cex=1, #box size
         col=c("tomato","steelblue","green","pink","yellow"),
         y.intersp=1.5)
}
plot_personal_graph(tf_imp_net)

#### integrating micro data to calcualte co-regulation of TF-TFs or TF2IGs correlation
micro_data<-read.table("myresearch/network/data/7all_gene_expresion.txt",sep = "\t")
hh<-t(micro_data)
co_regulated_IGs<-colSums(mmm)[colSums(mmm)>=2]
co_regulated_tfs_corrs<-c()

for (name in names(co_regulated_IGs)){
  co_regulated_tfs<-rownames(mmm)[mmm[,name]!=0]
  message(length(co_regulated_tfs))
  if (length(intersect(colnames(hh),co_regulated_tfs))>=2){
    temp_corrs<-cor(hh[,intersect(colnames(hh),co_regulated_tfs)])
    diag(temp_corrs)<-0
    #message(mean(temp_corrs[upper.tri(temp_corrs)]))
    message(max(temp_corrs))
    print(which(temp_corrs == max(temp_corrs), arr.ind = TRUE))
    co_regulated_tfs_corrs<-c(co_regulated_tfs_corrs,temp_corrs[upper.tri(temp_corrs)])
  }
} 

tfs_random_corrs<-c()
for (m in as.numeric(co_regulated_IGs)){
  message(m)
  random_tfs<-sample(rownames(mmm),m)
  if (length(intersect(colnames(hh),random_tfs))>=2){
    temp_corrs<-cor(hh[,intersect(colnames(hh),random_tfs)])
    message(mean(temp_corrs[upper.tri(temp_corrs)]))
    tfs_random_corrs<-c(tfs_random_corrs,temp_corrs[upper.tri(temp_corrs)])
  }
}

tfs_corrs<-data.frame(corr=c(co_regulated_tfs_corrs,tfs_random_corrs),
                      class=c(rep("coregulator",length(co_regulated_tfs_corrs)),
                              rep("non-coregulator",length(tfs_random_corrs)))
                      )

library(dplyr)
mu<-tfs_corrs%>%group_by(class)%>%summarise_all(mean)

p<-ggplot(tfs_corrs,aes(x=corr,fill=class))+
  geom_density(alpha = 0.5)+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = c(0.8,0.9),legend.text = element_text(size = 12))+
  xlab("Expression correlation")+
  ylab("Density")+
  geom_vline(data=mu, aes(xintercept=corr, color=class),
             linetype="dashed",size=1)+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"))

#### co-regulated imprinted genes expression correlation analysis

co_regulated_TFs<-rowSums(mmm)[rowSums(mmm)>=2]  # these tfs regulated more than two genes
co_regulated_IGs_corrs<-c()

for (name in names(co_regulated_TFs)){
  co_regulated_IGs<-colnames(mmm)[mmm[name,]!=0]
  message(length(co_regulated_IGs))
  if (length(intersect(colnames(hh),co_regulated_IGs))>=2){  ## some genes not in mircro expresison
    temp_corrs<-cor(hh[,intersect(colnames(hh),co_regulated_IGs)])
    message(mean(temp_corrs[upper.tri(temp_corrs)]))
    co_regulated_IGs_corrs<-c(co_regulated_IGs_corrs,temp_corrs[upper.tri(temp_corrs)])
  }
} 

IGs_random_corrs<-c()
for (m in as.numeric(co_regulated_TFs)){
  message(m)
  random_IGs<-sample(colnames(mmm),m)
  if (length(intersect(colnames(hh),random_IGs))>=2){
    temp_corrs<-cor(hh[,intersect(colnames(hh),random_IGs)])
    message(mean(temp_corrs[upper.tri(temp_corrs)]))
    IGs_random_corrs<-c(IGs_random_corrs,temp_corrs[upper.tri(temp_corrs)])
  }
}
### plot
IGs_corrs<-data.frame(corr=c(co_regulated_IGs_corrs,IGs_random_corrs),
                      class=c(rep("coregulated",length(co_regulated_IGs_corrs)),
                              rep("non-coregulated",length(IGs_random_corrs)))
)
p<-ggplot(IGs_corrs,aes(y=corr,x=class,fill=class))+
  geom_boxplot()+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  ylab("Expression correlation")+
  xlab("")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"))+
  geom_signif(comparisons = list(c("coregulated","non-coregulated")))

## 7.9 transcription tissue specificty
gene_specificity<-function(expression_matrix){
  y<-apply(expression_matrix,1,function(x)sum(1-x/max(x))/(ncol(expression_matrix)-1))
  return(y)
}

endosperm_ratio<-function(expression_matrix){
  y<-apply(expression_matrix,1,function(x)(x[1]+x[2]+x[3])/(3*sum(x)))
  return(y)
}

rna_expression<-read.table("/data1/SRA/Arabodopsis_rna/3all_tissues_integration/reads_cpm.txt",header = TRUE)
#micro_pp<-read.table("other_resources/pp_micro_arrary/GSE69995_re-analyzed_data_matrix.txt",header = TRUE,sep="\t",row.names = 1)
gs<-gene_specificity(rna_expression)
er<-endosperm_ratio(rna_expression)

df_IGs<-data.frame(Specificity=gs[colnames(mmm)],Relative_Abandunce=er[colnames(mmm)],class="IGs")
df_TFs<-data.frame(Specificity=gs[rownames(mmm)],Relative_Abandunce=er[rownames(mmm)],class="TFs")
rownames(df_TFs)<-rownames(mmm)
df_TFs<-na.omit(df_TFs)
#df_All<-data.frame(Specificity=gs,Relative_Abandunce=er,class="All")
df<-rbind(df_IGs,df_TFs)#,df_All)

df$color[df$class=="IGs"]<-"tomato"
df$color[df$class=="TFs"]<-"steelblue"
df$id<-rownames(df)
#df$color[is.na(df$color)]<-"gray"

library(ggrepel)
plot(df$Specificity,df$Relative_Abandunce,col=adjustcolor(df$color,alpha=0.2),pch=16)
ggplot(df,aes(x=Specificity,y=Relative_Abandunce,color=class))+geom_point(alpha=0.5)+
  geom_smooth(method = "lm")+
  theme_bw()+
  theme(text=element_text(size=20),legend.position = c(0.15,0.9))+
  xlab("Expression specificity")+
  ylab("Expression in endopem")+
  geom_label_repel(data = subset(df,df$Specificity>0.8&df$Relative_Abandunce>0.3),
                   aes(label = id,
                       fill = class), 
                       color = 'white',
                       size = 3.5,
                       segment.color="black",
                   point.padding=0.5)+
  ## 0.66 vs 0.21
  
# tf and imp correlation
tf_g_gs<-data.frame(x=gs[tf_g$Var1],y=gs[tf_g$Var2])
ggplot(tf_g_gs,aes(x=x,y=y))+geom_point(alpha=0.5)+
  geom_smooth(method = "lm")+
  geom_density_2d()+
  theme_bw()+
  theme(text=element_text(size=20),legend.position = c(0.15,0.9))+
  xlab("TF expression specificity")+
  ylab("IG expression specificity")+
  annotate(x=0.3,y=0.95,
           label=paste("R = ", round(cor(tf_g_gs$x,tf_g_gs$y,use = "complete"),2)),
           geom="text", size=8, col="darkred")





#7.12 ######### orthologous imprinted genes annalsis
######### rice ortholgous imprinted genes

Compare_orthologous_genes<-function(imprinted_genes_file_in_other_species, orthologous_file,calculation="degree",species="Rice"){
  library(ggsignif)
  library(ggplot2)
  message(c("Species:",species))
  ## get IG, PEG, MEG from other species
  other_imprinted_data<-read.table(imprinted_genes_file_in_other_species,stringsAsFactors = FALSE) 
  other_imprinted_genes<-unique(other_imprinted_data$V1) # imprinted genes list in other species
  maternal_genes<-unique(other_imprinted_data$V1[other_imprinted_data$V2=="m"])
  paternal_genes<-unique(other_imprinted_data$V1[other_imprinted_data$V2=="f"])
  
  
  ##arabidopsis orthologs of impritned genes of other species 
  ortholougs<-read.table(orthologous_file,stringsAsFactors = FALSE) #39072 from biomart
  imp_orthlogous_genes<-unique(ortholougs$V1[ortholougs$V2%in% other_imprinted_genes])
  
  message(c("other imprinted genes:",length(other_imprinted_genes)))
  message(c("ortholgous to TAIR:",length(imp_orthlogous_genes)))
  message(c("ortholgous is imprinted in TAIR:",sum(imp_orthlogous_genes%in%imprinted_genes))) # 
  message(c("ortholgous is imprinted and in network in TAIR:",sum(imp_orthlogous_genes%in%imprinted_genes_in_network))) # 
 
  #imp_orthlogous_genes<-imp_orthlogous_genes[!imp_orthlogous_genes%in%imprinted_genes_in_network] # whether drop the othologous genes also imprinted in arabidopsis
  message(c("ortholgous is not imprinted and in network in TAIR:",length(imp_orthlogous_genes)))
  imp_orthlogous_genes_in_ara_network <- intersect(imp_orthlogous_genes,names(V(g)))
  
  if (calculation=="degree"){ ### calulated degree for different set of genes
                            message("degree calculation!")
                            df_ggplot1<-data.frame(degree=degree(g,names(V(g))[!names(V(g))%in%imprinted_genes_in_network]),class="All genes")
                            df_ggplot2<-data.frame(degree=degree(g,imp_orthlogous_genes_in_ara_network),class="orthologous")
                            df_ggplot3<-data.frame(degree=degree(g,imprinted_genes_in_network),class="Imprinted_genes")
  
                            ##plot for ggplot
                            df<-rbind(df_ggplot1,df_ggplot2,df_ggplot3)
                            df$species<-species
                            message(c(nrow(df),"-",ncol(df)))
                            p<-ggplot(df,aes(x=class,y=degree,fill=class))+geom_boxplot()+
                                              geom_signif(comparisons = list(c("All genes","orthologous"),c("Imprinted_genes","orthologous")),
                                              test="wilcox.test", test.args=list(alternative="greater"),
                                              step_increase = 0.05,tip_length = 0.01)+
                                              theme_bw(base_size = 20)+
                                              scale_x_discrete(labels=c("AG","OG","IG"))+
                                              scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
                                              theme(legend.position="none",axis.title.x=element_blank(),
                                                    #axis.text.x=element_text(angle = 45,vjust=0.7),
                                                    plot.title = element_text(hjust = 0.5,size=20,face = "bold"))+
                                              ylab("Degree")+ggtitle(species)
  
                            
  } ## d
  
  if(calculation=="shortest.path"){ 
                                  message("shortest.paths calculating!")
                                  #df_ggplot1<-data.frame(calculation = as.numeric(shortest.paths(g,names(V(g))[!names(V(g))%in%imprinted_genes_in_network],weights = NA),class="All genes"))
                                  
                                  d2<-shortest.paths(g,v = imp_orthlogous_genes_in_ara_network,to =imp_orthlogous_genes_in_ara_network, weights = NA)
                                  d2<-d2[upper.tri(d2)]
                                  d2<-d2[is.finite(d2)]
                                  df_ggplot2<-data.frame(calculation = d2,class="orthologous")
                                  
                                  d3<-shortest.paths(g,v = imprinted_genes_in_network,to =imprinted_genes_in_network, weights = NA)
                                  d3<-d3[upper.tri(d3)]
                                  d3<-d3[is.finite(d3)]
                                  df_ggplot3<-data.frame(calculation = d3,class="Imprinted_genes")
   
                                  ##plot for ggplot
                                  df<-rbind(df_ggplot2,df_ggplot3)
                                  df$species<-species
                                  #p<-ggplot(df,aes(x=class,y=calculation,fill=class))+geom_boxplot()+
                                  p<-ggplot(df,aes(x=class,y=calculation,fill=class))+geom_violin()+
                                                 geom_signif(comparisons = list(c("Imprinted_genes","orthologous")),
                                                 test="wilcox.test", test.args=list(alternative="less"),
                                                 step_increase = 0.05,tip_length = 0.01)+
                                                 theme_bw(base_size = 20)+
                                                 scale_x_discrete(labels=c("OGs","IGs"))+
                                                 scale_fill_manual(values = c("#999999", "#E69F00"))+
                                                 theme(legend.position="none",axis.title.x=element_blank(),
                                                       plot.title = element_text(hjust = 0.5,size=20,face = "bold"))+ # title posistion
                                                 ylab("Shortest path")+ggtitle(species)
                                                 
                                  }  ## compare shortest path of IGs and OGs (have difference)
  
  if(calculation=="betweenness"){
                                message("betweenness calculating!")
                                betweenness_all_genes<-read.table("4network_analysis_result/2high_confidence_imp_result/0.1AFGN_betweenness_normalized.txt",stringsAsFactors = FALSE,header = TRUE)
                                rownames(betweenness_all_genes)<-betweenness_all_genes$genes
                                
                                df_ggplot1<-data.frame(calculation = betweenness_all_genes[names(V(g))[!names(V(g))%in%imprinted_genes_in_network],]$betweenness,class="All genes")
                                df_ggplot2<-data.frame(calculation = betweenness_all_genes[imp_orthlogous_genes_in_ara_network,]$betweenness,class="orthologous")
                                df_ggplot3<-data.frame(calculation = betweenness_all_genes[imprinted_genes_in_network,]$betweenness,class="Imprinted_genes")
                                
                                df<-rbind(df_ggplot1,df_ggplot2,df_ggplot3)
                                df$species<-species
                                #p<-ggplot(df,aes(x=class,y=calculation,fill=class))+geom_boxplot()+
                                p<-ggplot(df,aes(x=class,y=log(calculation,2),fill=class))+geom_violin()+
                                  geom_signif(comparisons = list(c("All genes","orthologous"),c("Imprinted_genes","orthologous")),
                                              test="wilcox.test", test.args=list(alternative="greater"),
                                              step_increase = 0.05,tip_length = 0.01)+
                                  theme_bw(base_size = 20)+
                                  scale_x_discrete(labels=c("AG","OG","IG"))+
                                  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
                                  theme(legend.position="none",axis.title.x=element_blank(),
                                        plot.title = element_text(hjust = 0.5,size=12,face = "bold"))+ # title posistion
                                  ylab("Betweenness")+ggtitle(species)+
                                  stat_summary(fun.y=median, geom="point", size=2, color="red")
                                
                                  } ## differentce of betweenness of IGs and OGs
  
  if(calculation=="shortest.path.compared.random"){
    message("shortest.path.compared.random")
    d2<-shortest.paths(g,v = imp_orthlogous_genes_in_ara_network,to =imp_orthlogous_genes_in_ara_network, weights = NA)
    d2<-d2[upper.tri(d2)]
    d2<-d2[is.finite(d2)]
    
    # simulate 1000 times
    fun <- function(x){
      require(igraph)
      rand_genes<-sample(names(V(g)),length(imp_orthlogous_genes_in_ara_network),replace = FALSE)
      d_random<-shortest.paths(g,rand_genes, rand_genes,weights = NA)
      d_random<-d_random[upper.tri(d_random)]
      d_random<-d_random[is.finite(d_random)]
      return(mean(d_random))
    }                              
    library(parallel)
    system.time(simulations<-mclapply(1:1000,fun,mc.cores = 7))
    simulations<-unlist(simulations)
    dff<-data.frame(simulations=simulations,species=species)
    p<-ggplot(dff,aes(x=simulations))+geom_histogram(fill="steelblue")+
      theme_bw()+
      theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5,size=20,face = "bold"))+
      ylab("Frequency")+
      xlab("Mean shortest path")+
      annotate("segment", x=mean(d2), xend=mean(d2), y=25, yend=0, color="black", size=1, arrow=arrow())+
      annotate("segment", x=median(d2), xend=median(d2), y=25, yend=0, color="red", size=1, arrow=arrow())+
      annotate("text", x = median(d2), y = 28, label = "Median")+
      annotate("text", x = mean(d2), y = 28, label = "Mean")+ggtitle(species)
  } ## mean shortest path of OG , with random gene mean distance
  
  if(calculation=="distance_to_imp"){
    message("distance_to_imp")
    d2<-shortest.paths(g,v = imp_orthlogous_genes_in_ara_network,to =imprinted_genes_in_network, weights = NA)
    d2<-d2[upper.tri(d2)]
    d2<-d2[is.finite(d2)]
   ## distence of OG to IGs
    fun <- function(x){
      require(igraph)
      rand_genes<-sample(names(V(g)),length(imp_orthlogous_genes_in_ara_network),replace = FALSE)
      d_random<-shortest.paths(g,rand_genes, imprinted_genes_in_network,weights = NA)
      d_random<-d_random[upper.tri(d_random)]
      d_random<-d_random[is.finite(d_random)]
      return(mean(d_random))
    }
    library(parallel)
    system.time(simulations<-mclapply(1:1000,fun,mc.cores = 7))
    simulations<-unlist(simulations)
    dff<-data.frame(simulations=simulations,species=species)
    ##plot
    p<-ggplot(dff,aes(x=simulations))+geom_histogram(fill="steelblue")+
      theme_bw()+
      theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5,size=20,face = "bold"))+
      ylab("Frequency")+
      xlab("Mean shortest path of OG to IGs")+
      annotate("segment", x=mean(d2), xend=mean(d2), y=25, yend=0, color="black", size=1, arrow=arrow())+
      annotate("segment", x=median(d2), xend=median(d2), y=25, yend=0, color="red", size=1, arrow=arrow())+
      annotate("text", x = median(d2), y = 28, label = "Median")+
      annotate("text", x = mean(d2), y = 28, label = "Mean")+ggtitle(species)
 } # shortest path of OG to IG
  
  if(calculation=="cluster_coefficient"){
    message("cluster coefficient calculating!")
    cluster_coefficient<-transitivity(g,weights = NA,type = "local")
    names(cluster_coefficient)<-names(V(g))
    
    df_ggplot1<-data.frame(calculation = cluster_coefficient[names(V(g))[!names(V(g))%in%imprinted_genes_in_network]],class="All genes")
    df_ggplot2<-data.frame(calculation = cluster_coefficient[imp_orthlogous_genes_in_ara_network],class="orthologous")
    df_ggplot3<-data.frame(calculation = cluster_coefficient[imprinted_genes_in_network],class="Imprinted_genes")
    
    df<-rbind(df_ggplot1,df_ggplot2,df_ggplot3)
    df$species<-species
    #p<-ggplot(df,aes(x=class,y=calculation,fill=class))+geom_boxplot()+
    p<-ggplot(df,aes(x=class,y=calculation,fill=class))+geom_boxplot(notch = TRUE)+
      geom_signif(comparisons = list(c("All genes","orthologous"),c("Imprinted_genes","orthologous")),
                  test="wilcox.test", test.args=list(alternative="greater"),
                  step_increase = 0.05,tip_length = 0.01)+
      theme_bw(base_size = 20)+
      scale_x_discrete(labels=c("AG","OG","IG"))+
      scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
      theme(legend.position="none",axis.title.x=element_blank(),
            plot.title = element_text(hjust = 0.5,size=12,face = "bold"))+ # title posistion
      ylab("Cluster coefficient")+ggtitle(species)+
      stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  }
  
  if(calculation=="OG_distane_to_imprinting_center"){
    IG2other<-shortest.paths(g,v = imprinted_genes_in_network,names(V(g)),weights = NA)
    other2IG_mean_distance<-colMeans(IG2other)
    IG_center<-unique(names(sort(other2IG_mean_distance)))[50] # 50 most closed genes to imprinted genes as imprinting center
    OG2center<-shortest.paths(graph = g,v = imp_orthlogous_genes_in_ara_network,to = IG_center,weights = NA)
    OG2center_mean_distance<-mean(OG2center[is.finite(OG2center)])
    return(OG2center_mean_distance)
  }
  return(p)
  
}

#degree calcualtion
p1<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/RICE/imp_by_papers/imp_2.txt",
                          orthologous_file = "/data1/SRA/rice_chip/3scripts/1pesudo_gene_prediction/2rna_seq_20173D/3orthlogs/3ara2rice",calculation = "degree",species = "O. sativa")
p2<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/MAIZE/maize_imp_by_papers/maize2+.imp",
                          orthologous_file = "~/MyResearch/genome.db/Maize/gtf/gene_id_convert/TAIR2maizev3_ortholgous",calculation = "degree",species = "Z. mays")
p3<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/lyrata/lyrata_maped.imp",
                          orthologous_file = "~/MyResearch/genome.db/Arabidopsis_lyrata/gtf/TAIR2lyrata",calculation = "degree",species = "A. lyrata")
p4<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sy/strong.sy.imp",
                          orthologous_file = "~/MyResearch/genome.db/Solanum_lycopersicum/gtf/TAIR2sly",calculation = "degree",species = "S. lycopersicum")
p5<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sbicolor/sb.impV2",
                          orthologous_file = "~/MyResearch/genome.db/Sbicolor/gtf/sb2tair.txt",calculation = "degree",species = "S. bicolor")
p6<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/castor.imp",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/2TAIR2R.comm.mRNA.txt",calculation = "degree",species = "R. communis")
p7<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/2imprinted_list",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/tair2c.rubella",calculation = "degree",species = "C. rubella")

#grid.arrange(plot1, plot2, nrow=1, ncol=2)
#library(ggpubr) 
#ggarrange(p1,p2,p3,p4,p5,p6,p7)
pd<-rbind(p1$data,p2$data,p3$data,p4$data,p5$data,p6$data,p7$data)
pd$species<-factor(pd$species,levels=c("A. lyrata","C. rubella","R. communis", "S. lycopersicum", "O. sativa", "S. bicolor", "Z. mays"))
p<-ggplot(pd,aes(x=class,y=log(degree,2),fill=class))+geom_boxplot(notch=TRUE)+
  facet_wrap(~ species,ncol = 7)+
  geom_signif(comparisons = list(c("All genes","orthologous"),
                                 c("Imprinted_genes","orthologous")),
              test="wilcox.test", 
              #test.args=list(alternative="greater"),
              #map_signif_level = TRUE,
              step_increase = 0.05,
              tip_length = 0.01)+
  theme_bw(base_size = 20)+
  scale_x_discrete(labels=c("AG","OG","IG"))+
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        #axis.text.x=element_text(angle = 45,vjust=0.7),
        plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        panel.spacing=unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12,face = "italic"))+
  ylab("Degree (log2)")

## shortest path calulation
p1<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/RICE/imp_by_papers/imp_2.txt",
                          orthologous_file = "/data1/SRA/rice_chip/3scripts/1pesudo_gene_prediction/2rna_seq_20173D/3orthlogs/3ara2rice",calculation = "shortest.path",species = "O. sativa")
p2<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/MAIZE/maize_imp_by_papers/maize2+.imp",
                          orthologous_file = "~/MyResearch/genome.db/Maize/gtf/gene_id_convert/TAIR2maizev3_ortholgous",calculation = "shortest.path",species = "Z. mays")
p3<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/lyrata/lyrata_maped.imp",
                          orthologous_file = "~/MyResearch/genome.db/Arabidopsis_lyrata/gtf/TAIR2lyrata",calculation = "shortest.path",species = "A. lyrata")
p4<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sy/strong.sy.imp",
                          orthologous_file = "~/MyResearch/genome.db/Solanum_lycopersicum/gtf/TAIR2sly",calculation = "shortest.path",species = "S. lycopersicum")
p5<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sbicolor/sb.impV2",
                          orthologous_file = "~/MyResearch/genome.db/Sbicolor/gtf/sb2tair.txt",calculation = "shortest.path",species = "S. bicolor")
p6<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/castor.imp",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/2TAIR2R.comm.mRNA.txt",calculation = "shortest.path",species = "R. communis")
p7<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/2imprinted_list",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/tair2c.rubella",calculation = "shortest.path",species = "C. rubella")

#grid.arrange(plot1, plot2, nrow=1, ncol=2)
#library(ggpubr) 
#ggarrange(p1,p2,p3,p4,p5)
pd<-rbind(p1$data,p2$data,p3$data,p4$data,p5$data,p6$data,p7$data)
pd$species<-factor(pd$species,levels=c("A. lyrata","C. rubella","R. communis", "S. lycopersicum", "O. sativa", "S. bicolor", "Z. mays"))

p<-ggplot(pd,aes(x=class,y=calculation,fill=class))+geom_violin()+
  facet_wrap(~ species,ncol = 7)+
  geom_signif(comparisons = list(c("Imprinted_genes","orthologous")),
              test="wilcox.test", 
              #test.args=list(alternative="greater"),
              #map_signif_level = TRUE,
              step_increase = 0.05,
              tip_length = 0.01)+
  theme_bw(base_size = 20)+
  scale_x_discrete(labels=c("OG","IG"))+
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        #axis.text.x=element_text(angle = 45,vjust=0.7),
        plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        panel.spacing=unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12,face = "italic"))+
  ylab("Shortest path")+
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "red")
##betweennesss

p1<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/RICE/imp_by_papers/imp_2.txt",
                              orthologous_file = "/data1/SRA/rice_chip/3scripts/1pesudo_gene_prediction/2rna_seq_20173D/3orthlogs/3ara2rice",calculation = "betweenness",species = "O. sativa")
p2<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/MAIZE/maize_imp_by_papers/maize2+.imp",
                              orthologous_file = "~/MyResearch/genome.db/Maize/gtf/gene_id_convert/TAIR2maizev3_ortholgous",calculation = "betweenness",species = "Z. mays")
p3<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/lyrata/lyrata_maped.imp",
                              orthologous_file = "~/MyResearch/genome.db/Arabidopsis_lyrata/gtf/TAIR2lyrata",calculation = "betweenness",species = "A. lyrata")
p4<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sy/strong.sy.imp",
                              orthologous_file = "~/MyResearch/genome.db/Solanum_lycopersicum/gtf/TAIR2sly",calculation = "betweenness",species = "S. lycopersicum")
p5<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sbicolor/sb.impV2",
                              orthologous_file = "~/MyResearch/genome.db/Sbicolor/gtf/sb2tair.txt",calculation = "betweenness",species = "S. bicolor")
p6<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/castor.imp",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/2TAIR2R.comm.mRNA.txt",calculation = "betweenness",species = "R. communis")
p7<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/2imprinted_list",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/tair2c.rubella",calculation = "betweenness",species = "C. rubella")

##
pd<-rbind(p1$data,p2$data,p3$data,p4$data,p5$data,p6$data,p7$data)
pd$species<-factor(pd$species,levels=c("A. lyrata","C. rubella","R. communis", "S. lycopersicum", "O. sativa", "S. bicolor", "Z. mays"))

p<-ggplot(pd,aes(x=class,y=log(calculation,2),fill=class))+geom_boxplot(notch = TRUE)+
  facet_wrap(~ species,ncol = 7)+
  geom_signif(comparisons = list(c("All genes", "orthologous"),c("Imprinted_genes","orthologous")),
              test="wilcox.test", 
              #test.args=list(alternative="greater"),
              #map_signif_level = TRUE,
              step_increase = 0.05,
              tip_length = 0.01)+
  theme_bw(base_size = 20)+
  scale_x_discrete(labels=c("AG","OG","IG"))+
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        #axis.text.x=element_text(angle = 45,vjust=0.7),
        plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        panel.spacing=unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12,face = "italic"))+
  ylab("Betweeness (log2)")



### distance compared with random

p1<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/RICE/imp_by_papers/imp_2.txt",
                              orthologous_file = "/data1/SRA/rice_chip/3scripts/1pesudo_gene_prediction/2rna_seq_20173D/3orthlogs/3ara2rice",calculation = "shortest.path.compared.random",species = "O. sativa")
p2<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/MAIZE/maize_imp_by_papers/maize2+.imp",
                              orthologous_file = "~/MyResearch/genome.db/Maize/gtf/gene_id_convert/TAIR2maizev3_ortholgous",calculation = "shortest.path.compared.random",species = "Z. mays")
p3<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/lyrata/lyrata_maped.imp",
                              orthologous_file = "~/MyResearch/genome.db/Arabidopsis_lyrata/gtf/TAIR2lyrata",calculation = "shortest.path.compared.random",species = "A. lyrata")
p4<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sy/strong.sy.imp",
                              orthologous_file = "~/MyResearch/genome.db/Solanum_lycopersicum/gtf/TAIR2sly",calculation = "shortest.path.compared.random",species = "S. lycopersicum")
p5<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sbicolor/sb.impV2",
                              orthologous_file = "~/MyResearch/genome.db/Sbicolor/gtf/sb2tair.txt",calculation = "shortest.path.compared.random",species = "S. bicolor")
p6<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/castor.imp",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/2TAIR2R.comm.mRNA.txt",calculation = "shortest.path.compared.random",species = "R. communis")
p7<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/2imprinted_list",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/tair2c.rubella",calculation = "shortest.path.compared.random",species = "C. rubella")

pd<-rbind(p1$data,p2$data,p3$data,p4$data,p5$data,p6$data,p7$data)
pd$species<-factor(pd$species,levels=c("A. lyrata","C. rubella","R. communis", "S. lycopersicum", "O. sativa", "S. bicolor", "Z. mays"))

arrow_pos = data.frame(name=c(sum(p1$data$simulations<=p1$layers[[2]]$data[,1])/1000,
                              sum(p2$data$simulations<=p2$layers[[2]]$data[,1])/1000,
                              sum(p3$data$simulations<=p3$layers[[2]]$data[,1])/1000,
                              sum(p4$data$simulations<=p4$layers[[2]]$data[,1])/1000,
                              sum(p5$data$simulations<=p5$layers[[2]]$data[,1])/1000,
                              sum(p6$data$simulations<=p6$layers[[2]]$data[,1])/1000,
                              sum(p7$data$simulations<=p7$layers[[2]]$data[,1])/1000),
                      species=as.character(unique(pd$species)),
                       xvar=c(p1$layers[[2]]$data[,1],
                              p2$layers[[2]]$data[,1],
                              p3$layers[[2]]$data[,1],
                              p4$layers[[2]]$data[,1],
                              p5$layers[[2]]$data[,1],
                              p6$layers[[2]]$data[,1],
                              p7$layers[[2]]$data[,1]),
                       yvar=c(p1$layers[[2]]$data[,3],
                              p2$layers[[2]]$data[,3],
                              p3$layers[[2]]$data[,3],
                              p4$layers[[2]]$data[,3],
                              p5$layers[[2]]$data[,3],
                              p6$layers[[2]]$data[,3],
                              p7$layers[[2]]$data[,3])
                       )

arrow_pos$species = factor(arrow_pos$species,levels=levels(pd$species))

p<-ggplot(pd,aes(x=simulations))+geom_histogram(fill="#E69F00",bins = 100)+
  facet_wrap(~ species,ncol = 7)+
  theme_bw(base_size = 12)+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        #axis.text.x=element_text(angle = 45,vjust=0.7),
        plot.title = element_text(hjust = 0.5,size=12,face = "bold"),
        panel.spacing=unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12,face = "italic"))+
  ylab("Frequency")+
  xlab("Mean shortest path")+
  geom_segment(data = arrow_pos, aes(x=xvar,xend=xvar, y=yvar, yend=0), color="black", size=1, arrow=arrow())+
  geom_text(data = arrow_pos, aes(label=paste("p=",name),x=xvar,y=yvar+2))


#library(ggpubr) 
#ggarrange(p1,p2,p3,p4,p5)

### orthologous to imp
p1<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/RICE/imp_by_papers/imp_2.txt",
                              orthologous_file = "/data1/SRA/rice_chip/3scripts/1pesudo_gene_prediction/2rna_seq_20173D/3orthlogs/3ara2rice",calculation = "distance_to_imp",species = "O. sativa")
p2<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/MAIZE/maize_imp_by_papers/maize2+.imp",
                              orthologous_file = "~/MyResearch/genome.db/Maize/gtf/gene_id_convert/TAIR2maizev3_ortholgous",calculation = "distance_to_imp",species = "Z. mays")
p3<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/lyrata/lyrata_maped.imp",
                              orthologous_file = "~/MyResearch/genome.db/Arabidopsis_lyrata/gtf/TAIR2lyrata",calculation = "distance_to_imp",species = "A. lyrata")
p4<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sy/strong.sy.imp",
                              orthologous_file = "~/MyResearch/genome.db/Solanum_lycopersicum/gtf/TAIR2sly",calculation = "distance_to_imp",species = "S. lycopersicum")
p5<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sbicolor/sb.impV2",
                              orthologous_file = "~/MyResearch/genome.db/Sbicolor/gtf/sb2tair.txt",calculation = "distance_to_imp",species = "S. bicolor")
p6<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/castor.imp",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/2TAIR2R.comm.mRNA.txt",calculation = "distance_to_imp",species = "R. communis")
p7<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/2imprinted_list",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/tair2c.rubella",calculation = "distance_to_imp",species = "C. rubella")

pd<-rbind(p1$data,p2$data,p3$data,p4$data,p5$data,p6$data,p7$data)
pd$species<-factor(pd$species,levels=c("A. lyrata","C. rubella","R. communis", "S. lycopersicum", "O. sativa", "S. bicolor", "Z. mays"))

arrow_pos = data.frame(name=c(sum(p1$data$simulations<=p1$layers[[2]]$data[,1])/1000,
                              sum(p2$data$simulations<=p2$layers[[2]]$data[,1])/1000,
                              sum(p3$data$simulations<=p3$layers[[2]]$data[,1])/1000,
                              sum(p4$data$simulations<=p4$layers[[2]]$data[,1])/1000,
                              sum(p5$data$simulations<=p5$layers[[2]]$data[,1])/1000,
                              sum(p6$data$simulations<=p6$layers[[2]]$data[,1])/1000,
                              sum(p7$data$simulations<=p7$layers[[2]]$data[,1])/1000),
                       species=as.character(unique(pd$species)),
                       xvar=c(p1$layers[[2]]$data[,1],
                              p2$layers[[2]]$data[,1],
                              p3$layers[[2]]$data[,1],
                              p4$layers[[2]]$data[,1],
                              p5$layers[[2]]$data[,1],
                              p6$layers[[2]]$data[,1],
                              p7$layers[[2]]$data[,1]),
                       yvar=c(p1$layers[[2]]$data[,3],
                              p2$layers[[2]]$data[,3],
                              p3$layers[[2]]$data[,3],
                              p4$layers[[2]]$data[,3],
                              p5$layers[[2]]$data[,3],
                              p6$layers[[2]]$data[,3],
                              p7$layers[[2]]$data[,3])
)

p<-ggplot(pd,aes(x=simulations))+geom_histogram(fill="#E69F00",bins = 100)+
  facet_wrap(~ species,ncol = 7)+
  theme_bw(base_size = 12)+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        #axis.text.x=element_text(angle = 45,vjust=0.7),
        plot.title = element_text(hjust = 0.5,size=12,face = "bold"),
        panel.spacing=unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12,face = "italic"))+
  ylab("Frequency")+
  xlab("Mean shortest path to IGs")+
  geom_segment(data = arrow_pos, aes(x=xvar,xend=xvar, y=yvar, yend=0), color="black", size=1, arrow=arrow())+
  geom_text(data = arrow_pos, aes(label=paste("p=",name),x=xvar,y=yvar+2))

### cluster coefficient
##

p1<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/RICE/imp_by_papers/imp_2.txt",
                              orthologous_file = "/data1/SRA/rice_chip/3scripts/1pesudo_gene_prediction/2rna_seq_20173D/3orthlogs/3ara2rice",calculation = "cluster_coefficient",species = "O. sativa")
p2<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/MAIZE/maize_imp_by_papers/maize2+.imp",
                              orthologous_file = "~/MyResearch/genome.db/Maize/gtf/gene_id_convert/TAIR2maizev3_ortholgous",calculation = "cluster_coefficient",species = "Z. mays")
p3<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/lyrata/lyrata_maped.imp",
                              orthologous_file = "~/MyResearch/genome.db/Arabidopsis_lyrata/gtf/TAIR2lyrata",calculation = "cluster_coefficient",species = "A. lyrata")
p4<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sy/strong.sy.imp",
                              orthologous_file = "~/MyResearch/genome.db/Solanum_lycopersicum/gtf/TAIR2sly",calculation = "cluster_coefficient",species = "S. lycopersicum")
p5<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sbicolor/sb.impV2",
                              orthologous_file = "~/MyResearch/genome.db/Sbicolor/gtf/sb2tair.txt",calculation = "cluster_coefficient",species = "S. bicolor")
p6<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/castor.imp",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/2TAIR2R.comm.mRNA.txt",calculation = "cluster_coefficient",species = "R. communis")
p7<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/2imprinted_list",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/tair2c.rubella",calculation = "cluster_coefficient",species = "C. rubella")

##
pd<-rbind(p1$data,p2$data,p3$data,p4$data,p5$data,p6$data,p7$data)
pd$species<-factor(pd$species,levels=c("A. lyrata","C. rubella","R. communis", "S. lycopersicum", "O. sativa", "S. bicolor", "Z. mays"))

p<-ggplot(pd,aes(x=class,y=calculation,fill=class))+geom_boxplot(notch = TRUE)+
  facet_wrap(~ species,ncol = 7)+
  geom_signif(comparisons = list(c("All genes", "orthologous"),c("Imprinted_genes","orthologous")),
              test="wilcox.test", 
              #test.args=list(alternative="greater"),
              #map_signif_level = TRUE,
              step_increase = 0.05,
              tip_length = 0.01)+
  theme_bw(base_size = 20)+
  scale_x_discrete(labels=c("AG","OG","IG"))+
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        #axis.text.x=element_text(angle = 45,vjust=0.7),
        plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        panel.spacing=unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12,face = "italic"))+
  ylab("Cluster coefficient")

### OG distance to imprinting center 
p1<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/RICE/imp_by_papers/imp_2.txt",
                              orthologous_file = "/data1/SRA/rice_chip/3scripts/1pesudo_gene_prediction/2rna_seq_20173D/3orthlogs/3ara2rice",calculation = "OG_distane_to_imprinting_center",species = "O. sativa")
p2<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/MAIZE/maize_imp_by_papers/maize2+.imp",
                              orthologous_file = "~/MyResearch/genome.db/Maize/gtf/gene_id_convert/TAIR2maizev3_ortholgous",calculation = "OG_distane_to_imprinting_center",species = "Z. mays")
p3<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/lyrata/lyrata_maped.imp",
                              orthologous_file = "~/MyResearch/genome.db/Arabidopsis_lyrata/gtf/TAIR2lyrata",calculation = "OG_distane_to_imprinting_center",species = "A. lyrata")
p4<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sy/strong.sy.imp",
                              orthologous_file = "~/MyResearch/genome.db/Solanum_lycopersicum/gtf/TAIR2sly",calculation = "OG_distane_to_imprinting_center",species = "S. lycopersicum")
p5<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sbicolor/sb.impV2",
                              orthologous_file = "~/MyResearch/genome.db/Sbicolor/gtf/sb2tair.txt",calculation = "OG_distane_to_imprinting_center",species = "S. bicolor")
p6<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/castor.imp",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/2TAIR2R.comm.mRNA.txt",calculation = "OG_distane_to_imprinting_center",species = "R. communis")
p7<-Compare_orthologous_genes(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/2imprinted_list",
                              orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/tair2c.rubella",calculation = "OG_distane_to_imprinting_center",species = "C. rubella")



## OG distance between species based on sicence <ds> = d(ab) - (d(aa)+d(bb))/2
V(g)$ara_imp[names(V(g))%in%maternal_imprint]="MEGs"
V(g)$ara_imp[names(V(g))%in%paternal_imprint]="PEGs"

OG_class<-function(imprinted_genes_file_in_other_species,orthologous_file){
  other_imprinted_data<-read.table(imprinted_genes_file_in_other_species,stringsAsFactors = FALSE) 
  other_imprinted_genes<-unique(other_imprinted_data$V1) # imprinted genes list in other species
  maternal_genes<-unique(other_imprinted_data$V1[other_imprinted_data$V2=="m"])
  paternal_genes<-unique(other_imprinted_data$V1[other_imprinted_data$V2=="f"])
  
  
  ##arabidopsis orthologs of impritned genes of other species 
  ortholougs<-read.table(orthologous_file,stringsAsFactors = FALSE) #39072 from biomart
  imp_orthlogous_genes<-unique(ortholougs$V1[ortholougs$V2%in% other_imprinted_genes])
  
  message(c("other imprinted genes:",length(other_imprinted_genes)))
  message(c("ortholgous to TAIR:",length(imp_orthlogous_genes)))
  message(c("ortholgous is imprinted in TAIR:",sum(imp_orthlogous_genes%in%imprinted_genes))) # 
  message(c("ortholgous is imprinted and in network in TAIR:",sum(imp_orthlogous_genes%in%imprinted_genes_in_network))) # 
  
  #imp_orthlogous_genes<-imp_orthlogous_genes[!imp_orthlogous_genes%in%imprinted_genes_in_network] # whether drop the othologous genes also imprinted in arabidopsis
  message(c("ortholgous is not imprinted and in network in TAIR:",length(imp_orthlogous_genes)))
  imp_orthlogous_genes_in_ara_network <- intersect(imp_orthlogous_genes,names(V(g)))
  return(imp_orthlogous_genes_in_ara_network)
}
#rice
rice_OG<-OG_class(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/RICE/imp_by_papers/imp_2.txt",orthologous_file = "/data1/SRA/rice_chip/3scripts/1pesudo_gene_prediction/2rna_seq_20173D/3orthlogs/3ara2rice")
V(g)$rice_imp=ifelse(names(V(g))%in%rice_OG,"rice_OG","other") #75
## mazie
maize_OG<-OG_class(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/MAIZE/maize_imp_by_papers/maize2+.imp",orthologous_file = "~/MyResearch/genome.db/Maize/gtf/gene_id_convert/TAIR2maizev3_ortholgous")
V(g)$maize_imp=ifelse(names(V(g))%in%maize_OG,"maize_OG","other")
#lyrata
lyrata_OG<-OG_class(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/lyrata/lyrata_maped.imp",orthologous_file = "~/MyResearch/genome.db/Arabidopsis_lyrata/gtf/TAIR2lyrata")
V(g)$lyrata_imp=ifelse(names(V(g))%in%lyrata_OG,"lyrata_OG","other")
# soly
sly_OG<-OG_class(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sy/strong.sy.imp",orthologous_file = "~/MyResearch/genome.db/Solanum_lycopersicum/gtf/TAIR2sly")
V(g)$sly_imp=ifelse(names(V(g))%in%sly_OG,"sly_OG","other")
# sbicolor
sb_OG<-OG_class(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/sbicolor/sb.impV2",orthologous_file = "~/MyResearch/genome.db/Sbicolor/gtf/sb2tair.txt")
V(g)$sb_imp = ifelse(names(V(g))%in%sb_OG,"sb_OG","other")
## castor
castor_OG<-OG_class(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/castor.imp", orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/caster_bean/2TAIR2R.comm.mRNA.txt")
V(g)$castor_imp= ifelse(names(V(g))%in%castor_OG,"castor_OG","other")
#rubella
rubella_OG<-OG_class(imprinted_genes_file_in_other_species = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/2imprinted_list",orthologous_file = "~/MyResearch/Imprinting_prediction/imprint_gene_list/capsla_rubela/tair2c.rubella")
V(g)$rubella_imp=ifelse(names(V(g))%in%rubella_OG,"rubella_OG","other")

module_distances<-function(gene_list1,gene_list2){
  d1_dis_matrix<-shortest.paths(graph = g,v = gene_list1,to = gene_list1,weights = NA)
  d1_dis<-mean(d1_dis_matrix[upper.tri(d1_dis_matrix)&is.finite(d1_dis_matrix)],na.rm = TRUE)
  
  d2_dis_matrix<-shortest.paths(graph = g,v = gene_list2,to = gene_list2,weights = NA)
  d2_dis<-mean(d2_dis_matrix[upper.tri(d2_dis_matrix)&is.finite(d2_dis_matrix)],na.rm = TRUE)
  
  d12_dis_matrix<-shortest.paths(graph = g,v = gene_list1,gene_list2,weights = NA)
  if(identical(rownames(d12_dis_matrix),colnames(d12_dis_matrix))){
    d12_dis<-mean(d12_dis_matrix[upper.tri(d12_dis_matrix)&is.finite(d12_dis_matrix)],na.rm = TRUE)
  }
  else{
    d12_dis<-mean(d12_dis_matrix[is.finite(d12_dis_matrix)],na.rm = TRUE)
  }
  
  
  D<-d12_dis-(d1_dis+d2_dis)/2
  return(D)
}
md<-module_distances(gene_list1 = names(V(g))%in%imprinted_genes_in_network,gene_list2 = names(V(g))[V(g)$castor_imp=="castor_OG"])

###plot orthologous genes
all_OGs<-c(names(V(g))[V(g)$castor_imp=="castor_OG"],names(V(g))[V(g)$rice_imp=="rice_OG"],names(V(g))[V(g)$maize_imp=="maize_OG"],
           names(V(g))[V(g)$lyrata_imp=="lyrata_OG"],names(V(g))[V(g)$castor_imp=="castor"],names(V(g))[V(g)$sly_imp=="sly_OG"],
           names(V(g))[V(g)$sb_imp=="sb_OG"],imprinted_genes_in_network)

sub_g<-induced_subgraph(graph = g,vids = unique(all_OGs))


plot_personal_graph<-function(graph_object){
  #V(sub_g)$color[names(V(sub_g))%in%focus_links$V2]<-"purple"
  #V(sub_g)$color<-"green"
  V(sub_g)$frame.color<-"white"
  V(sub_g)$color[names(V(sub_g))%in%imprinted_genes_in_network]<-"#CD3333" # vertex color 
  V(sub_g)$color[names(V(sub_g))%in%names(V(g))[V(g)$lyrata_imp=="lyrata_OG"]]<-"#FF7F00"
  V(sub_g)$color[names(V(sub_g))%in%names(V(g))[V(g)$rice_imp=="rice_OG"]]<-"#EEC591"
  V(sub_g)$color[names(V(sub_g))%in%names(V(g))[V(g)$maize_imp=="maize_OG"]]<-"#53868B"
  V(sub_g)$color[names(V(sub_g))%in%names(V(g))[V(g)$rubella_imp=="rubella_OG"]]<-"#458B00" # vertex color 
  V(sub_g)$color[names(V(sub_g))%in%names(V(g))[V(g)$castor_imp=="castor_OG"]]<-"#AB82FF"
  V(sub_g)$color[names(V(sub_g))%in%names(V(g))[V(g)$sly_imp=="sly_OG"]]<-"#0000FF"
  V(sub_g)$color[names(V(sub_g))%in%names(V(g))[V(g)$sb_imp=="sb_OG"]]<-"#FF6EB4"
  
  V(sub_g)$size <- 4              # vertex size
  #V(sub_g)$label.cex <- 0.5         # vertex label size
  V(sub_g)$label<-NA
  V(sub_g)$label.color<-"black"
  
  E(sub_g)$color <- "black"        # edge color 
  E(sub_g)$width=2#E(sub_g)$weight/max(E(sub_g)$weight)
  
  plot.igraph(graph_object,layout=layout.fruchterman.reingold)
  legend('topleft',
         #bg="white",
         text.col="tomato",
         legend=c("Arabidopsis","Lyrata","Rice","Maize","Rubella","Castor", "Sly","Sbicolor"),
         pch=19, #shape
         box.lty=2, # 
         pt.cex= 3, #lines size 
         cex=1, #box size
         col=c("#CD3333", "#FF7F00", "#EEC591", "#53868B", "#458B00", "#AB82FF", "#0000FF", "#FF6EB4"),
         y.intersp=1.5)
}
plot_personal_graph(sub_g)

## simulation
simulations<-c()
for (m in seq(1:1000)){
  random_genes<-sample(names(V(g)),vcount(sub_g))
  random_net<-induced_subgraph(g,vids = random_genes)
  random_comp_size<-max(components(random_net)$csize)
  simulations<-c(simulations,random_comp_size)
}
p<-ggplot(data.frame(simulations=simulations),aes(x=simulations))+geom_histogram(fill="#E69F00",bins = 100)+
  theme_bw(base_size = 20)+
  ylab("Frequency")+
  xlab("Maximum component size")+
  annotate("segment", x=max(components(sub_g)$csize), xend=max(components(sub_g)$csize), y=25, yend=0, color="red", size=2, arrow=arrow())+
  annotate("text", x = max(components(sub_g)$csize), y = 28, label = "Observed")




######### populaiton data
#########7.13 population data
popu<-read.csv("other_sourse/Ath_chr_protein_coding.tsv",header = TRUE,sep = "\t")
rownames(popu)<-popu$gene_id
popu$gene_pi[popu$gene_id%in%imprinted_genes_in_network]
#boxplot(popu$promoter_pi[popu$gene_id%in%imprinted_genes],popu$promoter_pi,popu$promoter_pi[popu$gene_id%in%maternal_imprint])
g_pop_common_common_genes<-intersect(names(V(g)),popu$gene_id)

df_ggplot1<-data.frame(pi=popu$promoter_pi[popu$gene_id%in%imprinted_genes],class="IG",region="Promoter")
df_ggplot2<-data.frame(pi=popu$promoter_pi[popu$gene_id%in%maternal_imprint],class="MEG",region="Promoter")
df_ggplot3<-data.frame(pi=popu$promoter_pi[popu$gene_id%in%paternal_imprint],class="PEG",region="Promoter")
df_ggplot4<-data.frame(pi=popu$promoter_pi,class="All genes",region="Promoter")

df_ggplot5<-data.frame(pi=popu$gene_pi[popu$gene_id%in%imprinted_genes],class="IG",region="Genebody")
df_ggplot6<-data.frame(pi=popu$gene_pi[popu$gene_id%in%maternal_imprint],class="MEG",region="Genebody")
df_ggplot7<-data.frame(pi=popu$gene_pi[popu$gene_id%in%paternal_imprint],class="PEG",region="Genebody")
df_ggplot8<-data.frame(pi=popu$gene_pi,class="All genes",region="Genebody")


df<-rbind(df_ggplot1,df_ggplot2,df_ggplot3,df_ggplot4,df_ggplot5,df_ggplot6,df_ggplot7,df_ggplot8)
library(ggsignif)
ggplot(df,aes(x=class,y=pi,fill=class))+geom_boxplot()+
  geom_signif(comparisons = list(c("IG","All genes"),c("MEG","All genes"),c("PEG","All genes")),
              test="wilcox.test", test.args=list(alternative="greater"),step_increase = 0.05,tip_length = 0.01)+
  theme_bw(base_size = 20)+ facet_wrap(.~region,ncol = 2,scales = "free")+
  #scale_x_discrete(labels=c("All genes","Maize orthologs","Imprinted genes"))+
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9","#9370DB"))+
  theme(legend.position="none")+ylab("Pi")


#####7.14  imprinted genes and essential genes
#### imprinted genes and essential genes

gene_class<-read.table("myresearch/network/data/gene_classes.txt",header = FALSE,sep="\t",stringsAsFactors = FALSE)
ESN<-gene_class$V1[gene_class$V2=="1ESN"] #841
ESN<-intersect(ESN,names(V(g))) #820

MRP<-gene_class$V1[gene_class$V2=="2MRP"] #2149
MRP<-intersect(MRP,names(V(g))) #2009

CLB<-gene_class$V1[gene_class$V2=="3CLB"] #399
CLB<-intersect(CLB,names(V(g))) #385

CND<-gene_class$V1[gene_class$V2=="4CND"] #586
CND<-intersect(CND,names(V(g))) #574

NO_PHE5NOTYPE<-gene_class$V1[gene_class$V2=="5NO PHE5NOTYPE"] #369
NO_PHE5NOTYPE<-intersect(NO_PHE5NOTYPE,names(V(g))) #341

ESN2IGs <-as.numeric(shortest.paths(graph = g,v = imprinted_genes_in_network,ESN,weights = NA))
ESN2IGs_P <-as.numeric(shortest.paths(graph = g,v = paternal_imprint,ESN,weights = NA))
ESN2IGs_M <-as.numeric(shortest.paths(graph = g,v = maternal_imprint,ESN,weights = NA))



MRP2IGs <-as.numeric(shortest.paths(graph = g,v = imprinted_genes_in_network,MRP,weights = NA))
CLB2IGs <-as.numeric(shortest.paths(graph = g,v = imprinted_genes_in_network,CLB,weights = NA))
CND2IGs <-as.numeric(shortest.paths(graph = g,v = imprinted_genes_in_network,CND,weights = NA))
NO_PHE5NOTYPE2IGs <-as.numeric(shortest.paths(graph = g,v = imprinted_genes_in_network,NO_PHE5NOTYPE,weights = NA))


##distance of imprinted genes to essential genes (Figure 3A)

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(ESN),replace = FALSE)
  d_random<-shortest.paths(g,imprinted_genes_in_network, rand_genes,weights = NA)
  return(mean(as.numeric(d_random)[is.finite(as.numeric(d_random))]))
}                              
library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 60))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 500)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean shortest path")+
  annotate("segment", x=mean(ESN2IGs[is.finite(ESN2IGs)]), xend=mean(ESN2IGs[is.finite(ESN2IGs)]), y=5, yend=0, color="black", size=2, arrow=arrow())+
  annotate("text", x = mean(ESN2IGs[is.finite(ESN2IGs)]), y = 6, label = "IGs")+
  annotate("segment", x=mean(ESN2IGs_P[is.finite(ESN2IGs_P)]), xend=mean(ESN2IGs_P[is.finite(ESN2IGs_P)]), y=5, yend=0, color="blue", size=2, arrow=arrow())+
  annotate("text", x = mean(ESN2IGs_P[is.finite(ESN2IGs_P)]), y = 6, label = "PEGs")+
  annotate("segment", x=mean(ESN2IGs_M[is.finite(ESN2IGs_M)]), xend=mean(ESN2IGs_M[is.finite(ESN2IGs_M)]), y=5, yend=0, color="tomato", size=2, arrow=arrow())+
  annotate("text", x = mean(ESN2IGs_M[is.finite(ESN2IGs_M)]), y = 6, label = "MEGs")

## plot
sub_g<-induced.subgraph(g,vids = c(imprinted_genes_in_network,ESN))
#
sub_g <- induced.subgraph(sub_g, names(membership(components(sub_g))[membership(components(sub_g))==1]))
plot_personal_graph<-function(graph_object){
  V(sub_g)$color[names(V(sub_g))%in%maternal_imprint]<-"tomato" # vertex color 
  V(sub_g)$color[names(V(sub_g))%in%paternal_imprint]<-"steelblue"
  V(sub_g)$color[names(V(sub_g))%in%ESN]<-"orange"
  
  V(sub_g)$size <- 2              # vertex size
  V(sub_g)$label.cex <- 0.8          # vertex label size
  V(sub_g)$label.color<-"black"
  E(sub_g)$color <- "gray"        # edge color 
  E(sub_g)$width=1
  V(sub_g)$label<-NA
  
  plot.igraph(sub_g,layout=layout.fruchterman.reingold,
              vertex.frame.color= "white")
  legend('topleft',
         legend=c("Maternal","Paternal","Essential"),
         pch=19, #shape
         box.lty=2, # 
         pt.cex= 3, #lines size 
         cex=1, #box size
         col=c("tomato","steelblue","orange"),
         y.intersp=1.5)
}
plot_personal_graph(sub_g)

## partners overlaped with essential genes
intersect_with_partners<-c()
for (r in seq(1:1000)){
  message(r)
  intersect_with_partners<-c(intersect_with_partners,length(intersect(ESN,sample(names(V(g)),vcount(IGPN_sub_imp))))/vcount(IGPN_sub_imp))
}
real<-length(intersect(ESN,names(V(IGPN_sub_imp))))/vcount(IGPN_sub_imp)

ggplot(data.frame(overlap=intersect_with_partners),aes(x=overlap))+
  geom_histogram(bins=100,fill="steelblue")+
  theme_bw(base_size = 20)+
  annotate("segment", x=real, xend=real, y=10, yend=0, color="red", size=2, arrow=arrow())+
  annotate("text", x = real, y = 12, label = "Observed")+ylab("Frequency")+xlab("Percentage of genes being essential genes")
  

#### overlap with no phenotype
intersect_with_partners<-c()
for (r in seq(1:1000)){
  intersect_with_partners<-c(intersect_with_partners,length(intersect(NO_PHE5NOTYPE,sample(names(V(g)),vcount(IGPN_sub_imp)))) )
}
real<-length(intersect(NO_PHE5NOTYPE,names(V(IGPN_sub_imp))))

ggplot(data.frame(overlap=intersect_with_partners),aes(x=overlap))+
  geom_histogram(bins=20,fill="steelblue")+
  theme_bw(base_size = 20)+
  annotate("segment", x=real, xend=real, y=25, yend=0, color="red", size=2, arrow=arrow())+
  annotate("text", x = real, y = 30, label = "Observed")+ylab("Count")+xlab("Overlaps with no-phenotype genes")

#####
### 5.1 essential degree comparsion (Figure S7)
degree_essential <- degree(g,v = ESN) ## long time
mean(degree_essential)
#degree_random <- degree(g,v = sample(names(V(g)),length(essential_genes),replace = FALSE)) # long time

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(ESN),replace = FALSE)
  degree_random<-degree(g,rand_genes)
  return(mean(as.numeric(degree_random)))
}

library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 60))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins=200)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean degree")+
  annotate("segment", x=78.8, xend=78.9, y=5, yend=0, color="red", size=2, arrow=arrow())+
  annotate("text", x = 78.8, y = 7, label = "Essential genes")

##### betweenness of essential to random
between_essential <- IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%ESN] ## long time
mean(between_essential)
#degree_random <- degree(g,v = sample(names(V(g)),length(essential_genes),replace = FALSE)) # long time

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(ESN),replace = FALSE)
  betweenness_random<-IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%rand_genes]
  return(mean(betweenness_random))
}

library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 60))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 200)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean betweenness")+
  annotate("segment", x=mean(between_essential), xend=mean(between_essential), y=5, yend=0, color="red", size=2, arrow=arrow())+
  annotate("text", x = mean(between_essential), y =6, label = "Essential genes")

### the percentage of partners being essential 
partners_non_essential <-




### 6.pathogen genes analysis
### pathogen genes analysis

##### 6.1 indivudule kegg pathway validation 
require(pathview)
library(org.At.tair.db)
xx<- as.list(org.At.tairPATH2TAIR)
focused_kegg_Pathway<-"04626"
#sub_g<-induced_subgraph(g,vids = intersect(xx[[focused_kegg_Pathway]],names(V(g))))
#ath.dat.kegg <- sim.mol.data(mol.type="gene",id.type="tair",species="ath",nmol=3000)
#pv.out <- pathview(gene.data = ath.dat.kegg, gene.idtype="tair",pathway.id = focused_kegg_Pathway, species = "ath", out.suffix = "ath.kegg",kegg.native = T, same.layer=T)
pathogen_degree<- degree(g,v = intersect(xx[[focused_kegg_Pathway]],names(V(g))))
pathogen_genes<-intersect(xx[[focused_kegg_Pathway]],names(V(g)))

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(pathogen_genes),replace = FALSE)
  degree_random<-degree(g,rand_genes)
  return(mean(as.numeric(degree_random)))
}

library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 7))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue")+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean degree")+
  annotate("segment", x=102, xend=102, y=25, yend=0, color="black", size=2, arrow=arrow())

#random_degree <- degree(g,v = sample(names(V(g)),length(names(pathogen_degree)),replace = FALSE))
#boxplot(list(pathogen_degree,random_degree),names=c("Pathogen genes","Random"))

##6.2betweenness
between_pathogen <- IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%pathogen_genes] ## long time
mean(between_pathogen)
#degree_random <- degree(g,v = sample(names(V(g)),length(essential_genes),replace = FALSE)) # long time

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(pathogen_genes),replace = FALSE)
  betweenness_random<-IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%rand_genes]
  return(mean(betweenness_random))
}

library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 7))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 50)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean betweenness")+
  annotate("segment", x=mean(between_pathogen), xend=mean(between_pathogen), y=25, yend=0, color="red", size=2, arrow=arrow())





###6.3 eigenvector
eigenvector <-eigen_centrality(g,weights = NA)
pathogen_eigenvector<-eigenvector$vector[names(pathogen_degree)]

fun <- function(x){
  require(igraph)
  return(mean(sample(size = 1000,x = eigenvector$vector)))
}

library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 7))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue")+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Eigenvector")+
  annotate("segment", x=102, xend=102, y=25, yend=0, color="black", size=2, arrow=arrow())

#random_eigenvector <-eigenvector$vector[names(betweenness_random)]
#boxplot(list(log(pathogen_eigenvector,2),log(random_eigenvector,2)),names=c("Pathogen genes","Random"))


#####7.15  imprinted genes and pathogen genes
#### imprinted genes and pathogen genes shortest path
PATHOGEN2IGs <-shortest.paths(graph = g,v = imprinted_genes_in_network,pathogen_genes,weights = NA)
PATHOGEN2IGs<-PATHOGEN2IGs[upper.tri(PATHOGEN2IGs)]

PATHOGEN2IGs_P <-shortest.paths(graph = g,v = paternal_imprint,pathogen_genes,weights = NA)
PATHOGEN2IGs_P <- PATHOGEN2IGs_P[upper.tri(PATHOGEN2IGs_P)]

PATHOGEN2IGs_M <-shortest.paths(graph = g,v = maternal_imprint,pathogen_genes,weights = NA)
PATHOGEN2IGs_M<-PATHOGEN2IGs_M[upper.tri(PATHOGEN2IGs_M)]

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(pathogen_genes),replace = FALSE)
  d_random<-shortest.paths(g,imprinted_genes_in_network, rand_genes,weights = NA)
  d_random<-d_random[upper.tri(d_random)]
  return(mean(d_random[is.finite(d_random)]))
}                              
library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 7))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 50)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean shortest path to pathogen related genes")+
  #xlim(3,4)+
  annotate("segment", x=mean(PATHOGEN2IGs[is.finite(PATHOGEN2IGs)]), xend=mean(PATHOGEN2IGs[is.finite(PATHOGEN2IGs)]), y=25, yend=0, color="black", size=1, arrow=arrow(length = unit(0.3, "cm")))+
  annotate("text", x = mean(PATHOGEN2IGs[is.finite(PATHOGEN2IGs)]), y = 30, label = "IGs")+
  annotate("segment", x=mean(PATHOGEN2IGs_P[is.finite(PATHOGEN2IGs_P)]), xend=mean(PATHOGEN2IGs_P[is.finite(PATHOGEN2IGs_P)]), y=25, yend=0, color="blue", size=1, arrow=arrow(length = unit(0.3, "cm")))+
  annotate("text", x = mean(PATHOGEN2IGs_P[is.finite(PATHOGEN2IGs_P)]), y = 30, label = "PEGs")+
  annotate("segment", x=mean(PATHOGEN2IGs_M[is.finite(PATHOGEN2IGs_M)]), xend=mean(PATHOGEN2IGs_M[is.finite(PATHOGEN2IGs_M)]), y=25, yend=0, color="tomato", size=1, arrow=arrow(length = unit(0.3, "cm")))+
  annotate("text", x = mean(PATHOGEN2IGs_M[is.finite(PATHOGEN2IGs_M)]), y = 30, label = "MEGs")

### direct interaction pathogen genes
imprint_neighbours<-unique(names(unlist(neighborhood(graph = g,order = 1,nodes = imprinted_genes_in_network))))
imprint_neighbours_pathogen<-intersect(imprint_neighbours,pathogen_genes)
sub_gg<-induced_subgraph(graph = g,vids = imprinted_genes,imprint_neighbours_pathogen)

plot_personal_graph<-function(graph_object){
  V(graph_object)$color[names(V(graph_object))%in%maternal_imprint]<-"tomato" # vertex color 
  V(graph_object)$color[names(V(graph_object))%in%paternal_imprint]<-"steelblue"
  V(graph_object)$color[!names(V(graph_object))%in%imprinted_genes_in_network]<-"black"
  V(graph_object)$size <- log(degree(graph_object),2)*2              # vertex size
  V(graph_object)$label.cex <- 0.5          # vertex label size
  V(graph_object)$label.color<-"black"
  E(graph_object)$color <- "gray"        # edge color 
  E(graph_object)$width=1
  #V(graph_object)$label=NA
  
  plot.igraph(graph_object,layout=layout_as_tree,vertex.frame.color= "white")
  legend('topleft',
         legend=c("Maternal","Paternal","Pathogen"),
         pch=19, #shape
         box.lty=2, # 
         pt.cex= 3, #lines size 
         cex=1, #box size
         col=c("tomato","steelblue","black"),
         y.intersp=1.5)
}
plot_personal_graph(graph_object = sub_gg)

## maximum component size
nodes_in_maximum<-names(components(sub_gg)$membership)[components(sub_gg)$membership==1]
sub_ggg<-induced_subgraph(graph = g,vids = nodes_in_maximum)

###7.155 immune related genes (new from network)
immune_related_genes <- read.table("myresearch/network/data/900immune_related.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)
immune_related_genes <-unique(c(immune_related_genes$ida,immune_related_genes$idb))
immune_related_genes_in_network<-intersect(immune_related_genes,names(V(g)))


### degree simulation

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(immune_related_genes_in_network),replace = FALSE)
  degree_random<-degree(g,rand_genes)
  return(mean(as.numeric(degree_random)))
}

library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 60))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 200)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean degree")+
  annotate("segment", x=62, xend=62, y=5, yend=0, color="black", size=2, arrow=arrow())+
  annotate("text", x = 62, y = 6, label = "Immune genes")
### betweenness 

between_pathogen <- IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%immune_related_genes_in_network] ## long time
mean(between_pathogen)
#degree_random <- degree(g,v = sample(names(V(g)),length(essential_genes),replace = FALSE)) # long time

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(immune_related_genes_in_network),replace = FALSE)
  betweenness_random<-IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%rand_genes]
  return(mean(betweenness_random))
}

library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 60))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 200)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean betweenness")+
  annotate("segment", x=mean(between_pathogen), xend=mean(between_pathogen), y=5, yend=0, color="red", size=2, arrow=arrow())+
  annotate("text", x = mean(between_pathogen), y = 6, label = "Immune genes")

###### imprinted genes and pathogen genes shortest path
PATHOGEN2IGs <-shortest.paths(graph = g,v = imprinted_genes_in_network,immune_related_genes_in_network,weights = NA)
PATHOGEN2IGs<-PATHOGEN2IGs[upper.tri(PATHOGEN2IGs)]

PATHOGEN2IGs_P <-shortest.paths(graph = g,v = paternal_imprint,immune_related_genes_in_network,weights = NA)
PATHOGEN2IGs_P <- PATHOGEN2IGs_P[upper.tri(PATHOGEN2IGs_P)]

PATHOGEN2IGs_M <-shortest.paths(graph = g,v = maternal_imprint,immune_related_genes_in_network,weights = NA)
PATHOGEN2IGs_M<-PATHOGEN2IGs_M[upper.tri(PATHOGEN2IGs_M)]

fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(immune_related_genes_in_network),replace = FALSE)
  d_random<-shortest.paths(g,imprinted_genes_in_network, rand_genes,weights = NA)
  d_random<-d_random[upper.tri(d_random)]
  return(mean(d_random[is.finite(d_random)]))
}                              
library(parallel)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 60))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 200)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean shortest path to pathogen related genes")+
  #xlim(3,4)+
  annotate("segment", x=mean(PATHOGEN2IGs[is.finite(PATHOGEN2IGs)]), xend=mean(PATHOGEN2IGs[is.finite(PATHOGEN2IGs)]), y=5, yend=0, color="black", size=1, arrow=arrow(length = unit(0.3, "cm")))+
  annotate("text", x = mean(PATHOGEN2IGs[is.finite(PATHOGEN2IGs)]), y = 6, label = "IGs")+
  annotate("segment", x=mean(PATHOGEN2IGs_P[is.finite(PATHOGEN2IGs_P)]), xend=mean(PATHOGEN2IGs_P[is.finite(PATHOGEN2IGs_P)]), y=5, yend=0, color="blue", size=1, arrow=arrow(length = unit(0.3, "cm")))+
  annotate("text", x = mean(PATHOGEN2IGs_P[is.finite(PATHOGEN2IGs_P)]), y = 6, label = "PEGs")+
  annotate("segment", x=mean(PATHOGEN2IGs_M[is.finite(PATHOGEN2IGs_M)]), xend=mean(PATHOGEN2IGs_M[is.finite(PATHOGEN2IGs_M)]), y=5, yend=0, color="tomato", size=1, arrow=arrow(length = unit(0.3, "cm")))+
  annotate("text", x = mean(PATHOGEN2IGs_M[is.finite(PATHOGEN2IGs_M)]), y = 6, label = "MEGs")

### visulazation
### direct interaction pathogen genes
imprint_neighbours<-unique(names(unlist(neighborhood(graph = g,order = 1,nodes = imprinted_genes_in_network))))
imprint_neighbours_pathogen<-intersect(imprint_neighbours,immune_related_genes_in_network)
sub_gg<-induced_subgraph(graph = g,vids = c(imprinted_genes_in_network,imprint_neighbours_pathogen))

plot_personal_graph<-function(graph_object){
  V(graph_object)$color[names(V(graph_object))%in%maternal_imprint]<-"tomato" # vertex color 
  V(graph_object)$color[names(V(graph_object))%in%paternal_imprint]<-"steelblue"
  V(graph_object)$color[!names(V(graph_object))%in%imprinted_genes_in_network]<-"lightgreen"
  V(graph_object)$size <- log(degree(graph_object),2)*2              # vertex size
  V(graph_object)$label.cex <- 0.5          # vertex label size
  V(graph_object)$label.color<-"black"
  E(graph_object)$color <- "gray"        # edge color 
  E(graph_object)$width=1
  #V(graph_object)$label=NA
  
  plot.igraph(graph_object,layout=layout_as_tree,vertex.frame.color= "white")
  legend('topleft',
         legend=c("Maternal","Paternal","Pathogen"),
         pch=19, #shape
         box.lty=2, # 
         pt.cex= 3, #lines size 
         cex=1, #box size
         col=c("tomato","steelblue","lightgreen"),
         y.intersp=1.5)
}
nodes_in_maximum<-names(components(sub_gg)$membership)[components(sub_gg)$membership==1]
sub_ggg<-induced_subgraph(graph = g,vids = nodes_in_maximum)
#plot_personal_graph(graph_object = sub_ggg)

## 7.16: cluster_one parse with kegg
#7.16: cluster_one parse with kegg
require(org.At.tair.db)
xx<- as.list(org.At.tairPATH2TAIR)

dat <- readLines("4network_analysis_result/2cluster_one/cluster_one_out.txt")
names(dat)<-seq(1,length(dat))
dat <- strsplit(dat, "\t")
module_genes<-unique(unlist(dat))

for (name in names(xx)){ # for loop kegg
  message(c(name,",",length(xx[[name]]),",",length(intersect(xx[[name]],module_genes))))
}

for (name in names(dat)){
  if (length(intersect(imprinted_genes_in_network,dat[[name]]))>0){
    message(name,",",length(intersect(imprinted_genes_in_network,dat[[name]])))
  }
}

## focues group
module_genes_focus<-dat$`908`
for (name in names(xx)){ # for loop kegg
  if (length(intersect(module_genes_focus,xx[[name]]))>0){
    message(c(name,",",length(intersect(module_genes_focus,xx[[name]]))))
  }
}


## 7.17: methylation related genes
#7.16: methylation related genes
methy_genes<-read.csv("other_resources/methylation_related_genes.csv",header = FALSE,stringsAsFactors = FALSE)
IGs2Methy<-shortest.paths(graph = g,v = methy_genes$V2,to = imprinted_genes_in_network,weights = NA)
IGs2Methy<-IGs2Methy[is.finite(IGs2Methy)]

PEGs2Methy<-shortest.paths(graph = g,v = methy_genes$V2,to = paternal_imprint,weights = NA)
PEGs2Methy<-PEGs2Methy[is.finite(PEGs2Methy)]

MEGs2Methy<-shortest.paths(graph = g,v = methy_genes$V2,to = maternal_imprint,weights = NA)
MEGs2Methy<-MEGs2Methy[is.finite(MEGs2Methy)]


fun <- function(x){
  require(igraph)
  rand_genes<-sample(names(V(g)),length(methy_genes),replace = FALSE)
  d_random<-shortest.paths(g,imprinted_genes_in_network, rand_genes,weights = NA)
  return(mean(as.numeric(d_random)[is.finite(as.numeric(d_random))]))
}                              
library(parallel)
library(ggplot2)
system.time(simulations<-mclapply(1:1000,fun,mc.cores = 7))
simulations<-unlist(simulations)
ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 50)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean shortest path")+
  annotate("segment", x=mean(IGs2Methy), xend=mean(IGs2Methy), y=25, yend=0, color="black", size=2, arrow=arrow())+
  annotate("text", x = mean(IGs2Methy), y = 30, label = "IGs")+
  annotate("segment", x=mean(PEGs2Methy[is.finite(PEGs2Methy)]), xend=mean(PEGs2Methy[is.finite(PEGs2Methy)]), y=25, yend=0, color="blue", size=2, arrow=arrow())+
  annotate("text", x = mean(PEGs2Methy[is.finite(PEGs2Methy)]), y = 30, label = "PEGs")+
  annotate("segment", x=mean(MEGs2Methy[is.finite(MEGs2Methy)]), xend=mean(MEGs2Methy[is.finite(MEGs2Methy)]), y=25, yend=0, color="tomato", size=2, arrow=arrow())+
  annotate("text", x = mean(MEGs2Methy[is.finite(MEGs2Methy)]), y = 30, label = "MEGs")

ggplot(as.data.frame(simulations),aes(x=simulations))+geom_histogram(fill="steelblue",bins = 50)+
  theme_bw()+
  theme(text = element_text(size=20))+
  ylab("Frequency")+
  xlab("Mean shortest path")+
  annotate("segment", x=median(IGs2Methy), xend=median(IGs2Methy), y=25, yend=0, color="black", size=2, arrow=arrow())+
  annotate("text", x = median(IGs2Methy), y = 30, label = "IGs")+
  annotate("segment", x=median(PEGs2Methy[is.finite(PEGs2Methy)]), xend=median(PEGs2Methy[is.finite(PEGs2Methy)]), y=25, yend=0, color="blue", size=2, arrow=arrow())+
  annotate("text", x = median(PEGs2Methy[is.finite(PEGs2Methy)]), y = 30, label = "PEGs")+
  annotate("segment", x=median(MEGs2Methy[is.finite(MEGs2Methy)]), xend=median(MEGs2Methy[is.finite(MEGs2Methy)]), y=25, yend=0, color="tomato", size=2, arrow=arrow())+
  annotate("text", x = median(MEGs2Methy[is.finite(MEGs2Methy)]), y = 30, label = "MEGs")

## 7.18 paralogs genes
##paralogs genes
paralogs_genes <-read.table("/home/wuzefeng/MyResearch/networks/2network_prediction/other_resources/paralogs/90IGs_paralogs.txt",header = TRUE,stringsAsFactors = FALSE,sep = "\t")
paralogs_genes<-unique(paralogs_genes$paralogue_genes)[unique(paralogs_genes$paralogue_genes) %in% names(V(g))][!unique(paralogs_genes$paralogue_genes)[unique(paralogs_genes$paralogue_genes) %in% names(V(g))]%in% imprinted_genes_in_network]
## degree comparison
require(ggsignif)
df<-data.frame(degree=c(degree(g,imprinted_genes_in_network),degree(g,paralogs_genes)),
               class=c(rep("IGs",length(imprinted_genes_in_network)),rep("IGPs",length(paralogs_genes))),stringsAsFactors = FALSE)
df$class<-factor(df$class,levels = c("IGs","IGPs"))
p_degree<-ggplot(df,aes(y=degree,x=class,fill=class))+geom_boxplot()+theme_bw(base_size = 20)+theme(legend.position = "none")+geom_signif(comparisons = list(c("IGs","IGPs")),test.args = "greater")+xlab("")+ylab("Degree")+scale_fill_manual(values = c("#E69F00", "#56B4E9"))

#### betweenness 
df<-data.frame(betweenness=c(imprinted_betweenness,IGFN_betweenness$betweenness[IGFN_betweenness$genes%in%paralogs_genes]),
               class=c(rep("IGs",length(imprinted_genes_in_network)),rep("IGPs",length(paralogs_genes))),stringsAsFactors = FALSE)
df$class<-factor(df$class,levels = c("IGs","IGPs"))
p_betweenness<-ggplot(df,aes(y=betweenness,x=class,fill=class))+geom_boxplot()+theme_bw(base_size = 20)+geom_signif(comparisons = list(c("IGs","IGPs")),test.args = "greater")+xlab("")+ylab("Betweenness")+theme(legend.position = "none")+scale_fill_manual(values = c("#E69F00", "#56B4E9"))

### closeness
df<-data.frame(closeness=c(closeness(g,imprinted_genes_in_network),closeness(g,paralogs_genes)),
               class=c(rep("IGs",length(imprinted_genes_in_network)),rep("IGPs",length(paralogs_genes))),stringsAsFactors = FALSE)
df$class<-factor(df$class,levels = c("IGs","IGPs"))
p_closeness<-ggplot(df,aes(y=closeness,x=class,fill=class))+geom_boxplot()+theme_bw(base_size = 20)+geom_signif(comparisons = list(c("IGs","IGPs")),test.args = "greater")+xlab("")+ylab("Closeness")+theme(legend.position = "none")+scale_fill_manual(values = c("#E69F00", "#56B4E9"))

library(gridExtra)
grid.arrange(p_degree, p_betweenness, p_closeness,nrow=1, ncol=3)


### 8.1 famod plot high confidence
edge_list.df<-as.data.frame(as_edgelist(IGPN),stringsAsFactors = FALSE)
colnames(edge_list.df)<-c("X1","X2")
gene_list<-unique(c(edge_list.df$X1,edge_list.df$X2)) #2746
gene_list_numbers <-seq(0,length(gene_list)-1)
gene_name2number<-data.frame(gene_list,gene_list_numbers,stringsAsFactors = FALSE)
edge_list.df$g1<-gene_name2number$gene_list_numbers[match(edge_list.df$X1,gene_name2number$gene_list)]
edge_list.df$g2<-gene_name2number$gene_list_numbers[match(edge_list.df$X2,gene_name2number$gene_list)]
imprinted_genes$color<-ifelse(imprinted_genes$V2=="m",1,2)
edge_list.df$c1<- imprinted_genes$color[match(edge_list.df$X1,imprinted_genes$V1)]
edge_list.df$c2<- imprinted_genes$color[match(edge_list.df$X2,imprinted_genes$V1)]
edge_list.df$c1[is.na(edge_list.df$c1)]=0
edge_list.df$c2[is.na(edge_list.df$c2)]=0
fanmond_input<-edge_list.df[,c(3,4,5,6)]
write.table(fanmond_input,"4network_analysis_result/2low_confidence_imp_result/network_with_low_imp.txt",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
## plot

m<-read.table("/home/wuzefeng/MyResearch/networks/2network_prediction/4network_analysis_result/2high_confidence_imp_result/imp_fanmod/motif.txt",sep="\t",header = TRUE)
m$num<-seq(1,nrow(m))
m$nz<-m$Z.Score/sqrt(sum(m$Z.Score^2))
ggplot(m,aes(x=num,y=nz,group=1))+ylim(-0.8,0.8)+
  geom_vline(xintercept = 0)+
  theme_minimal(base_size = 20)+
  scale_x_discrete(limits=seq(1,28),labels=seq(1,28))+xlab("")+
  ylab("Normalized Z-score")+
  geom_rect(aes(xmin=m$num-0.5,
                xmax=m$num+0.5,
                ymin=-Inf,
                ymax=Inf),
                fill = rep(c("gray70","#F8F8FF"),14))+
  geom_line(color="steelblue")+
  geom_point(size=6, shape=20,color="steelblue")+
  geom_hline(yintercept = 0)

### enrich motif instance extract
### motif instance extract
system("grep 011111112 fanmod_input.txt.csv.dump | awk 'BEGIN{FS=","}{print $2,$3,$4}' > 2motif18/motif18.txt")
focused_motif<- read.table("../2network_prediction/4network_analysis_result/2high_confidence_imp_result/imp_fanmod/2motif18/motif18.txt")
focused_motif$V4<-gene_name2number$gene_list[match(focused_motif$V1,gene_name2number$gene_list_numbers)]
focused_motif$V5<-gene_name2number$gene_list[match(focused_motif$V2,gene_name2number$gene_list_numbers)]
focused_motif$V6<-gene_name2number$gene_list[match(focused_motif$V3,gene_name2number$gene_list_numbers)]
class_genes<-function(x){
  if (x %in% paternal_imprint){
    return("P")
  }
  if (x %in% maternal_imprint){
    return("M")
  }
  if (! x %in% paternal_imprint & ! x %in% paternal_imprint){
    return("N")
  }
}
focused_motif$V7<-unlist(lapply(focused_motif$V4,class_genes))
focused_motif$V8<-unlist(lapply(focused_motif$V5,class_genes))
focused_motif$V9<-unlist(lapply(focused_motif$V6,class_genes))

focused_motif_sort<-as.data.frame(t(apply(focused_motif,1,function(x)x[c(4,5,6)][order(x[c(7,8,9)])])),stringsAsFactors = FALSE)
colnames(focused_motif_sort)<-sort(c("M","N","P"))
write.table(focused_motif_sort,file="4network_analysis_result/2high_confidence_imp_result/imp_fanmod/2motif18/2motif_18.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)


### low_confidence imp export famond 
edge_list.df<-as.data.frame(as_edgelist(IGPN),stringsAsFactors = FALSE)
colnames(edge_list.df)<-c("X1","X2")
gene_list<-unique(c(edge_list.df$X1,edge_list.df$X2)) #2746
gene_list_numbers <-seq(0,length(gene_list)-1)
gene_name2number<-data.frame(gene_list,gene_list_numbers,stringsAsFactors = FALSE)
edge_list.df$g1<-gene_name2number$gene_list_numbers[match(edge_list.df$X1,gene_name2number$gene_list)]
edge_list.df$g2<-gene_name2number$gene_list_numbers[match(edge_list.df$X2,gene_name2number$gene_list)]
imprinted_genes<-read.table("~/MyResearch/Imprinting_prediction/imprint_gene_list/6.1imprinted.list_ara",stringsAsFactors = FALSE)
imprinted_genes$color<-ifelse(imprinted_genes$V3=="m",1,2)
edge_list.df$c1<- imprinted_genes$color[match(edge_list.df$X1,imprinted_genes$V1)]
edge_list.df$c2<- imprinted_genes$color[match(edge_list.df$X2,imprinted_genes$V1)]
edge_list.df$c1[is.na(edge_list.df$c1)]=0
edge_list.df$c2[is.na(edge_list.df$c2)]=0
fanmond_input<-edge_list.df[,c(3,4,5,6)]
write.table(fanmond_input,"4network_analysis_result/2low_confidence_imp_result/network_with_low_imp.txt",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
## plot
m1<-read.table("/home/wuzefeng/MyResearch/networks/2network_prediction/4network_analysis_result/2low_confidence_imp_result/famod/2motfi_enrich.txt",sep="\t",header = TRUE)
m1$nz<-m1$Z.Score/sqrt(sum(m1$Z.Score^2))
m1<-m1[match(m$Adj,m1$Adj),]
m1$num<-seq(1,nrow(m1))

ggplot(m1,aes(x=num,y=nz,group=1))+ylim(-0.8,0.8)+
  geom_vline(xintercept = 0)+
  theme_minimal(base_size = 20)+
  scale_x_discrete(limits=seq(1,28),labels=seq(1,28))+xlab("")+
  ylab("Normalized Z-score")+
  geom_rect(aes(xmin=m1$num-0.5,
                xmax=m1$num+0.5,
                ymin=-Inf,
                ymax=Inf),
            fill = rep(c("gray70","#F8F8FF"),14))+
  geom_line(color="steelblue")+
  geom_point(size=6, shape=20,color="steelblue")+
  geom_hline(yintercept = 0)


### motif analysis
# count interecting peg and meg for each partner 
motif_peg <-c()
motif_meg <-c()


peg_threshold = 1
meg_threshold = 1

partners_peg_meg<-list()

for (p in partners){
  neib<-names(neighbors(graph = g,v = p))
  neib_peg <-intersect(paternal_imprint, neib)
  neib_meg <-intersect(maternal_imprint, neib)
  message("Neighbors number is: ",length(neib))
  message("-----PEGs Neighbors number is: ",length(neib_peg))
  message("-----MEGs Neighbors number is: ",length(neib_meg))      
  
  motif_peg<-c(motif_peg,length(neib_peg))
  motif_meg <-c(motif_meg, length(neib_meg))
  
  if (length(neib_meg)>=1 & length(neib_peg)>=1){
    edge_combins <- expand.grid(neib_meg,neib_peg)
    colnames(edge_combins)<-c("MEG","PEG")
    edge_combins$Partner<-p
    partners_peg_meg[[p]]<-edge_combins
  }
}
partners_peg_meg = do.call(rbind, partners_peg_meg)
partners_peg_meg$link <-apply(partners_peg_meg,1,function(x) ifelse(are.connected(g,x[1],x[2]),1,0)) # test link or not link between peg and meg

selected_partners<-partners_peg_meg%>%group_by(Partner,add = TRUE)%>%summarise(total = n(), type = n_distinct(link))%>%filter(type>1)  # keep both peg-meg linked and peg-med unlinked two types
selected_partners_peg_meg<-subset(partners_peg_meg,partners_peg_meg$Partner%in%selected_partners$Partner)

### import dn/ds value

#sta<-selected_partners_peg_meg%>%group_by(Partner,link)%>%summarise(n=n()) # numnber of linked and non-linked PEG-MEG for each partner
#plot(sta$n[sta$link==0],sta$n[sta$link!=0])

dnds<-read.table("myresearch/network/data/TAIR2lyrata.ds.mart_export.txt",sep="\t",header = TRUE,stringsAsFactors = FALSE)
dnds<-dnds[complete.cases(dnds),]
dnds<-unique(dnds)
dnds<-dnds %>% group_by(Gene.stable.ID) %>% slice(which.min(dS.with.Arabidopsis.lyrata)) # select minumum ds for same genes 
dnds$dnds<-dnds$dN.with.Arabidopsis.lyrata/dnds$dS.with.Arabidopsis.lyrata
rownames(dnds)<-dnds$Gene.stable.ID
selected_partners_peg_meg$delta_dnds<-apply(selected_partners_peg_meg,1,function(x) abs(dnds[x[1],]$dnds-dnds[x[2],]$dnds))
selected_partners_peg_meg<-selected_partners_peg_meg[complete.cases(selected_partners_peg_meg),]

require(ggsignif)
ggplot(selected_partners_peg_meg,aes(x=as.factor(link),y=delta_dnds))+
  geom_boxplot()+
  theme_bw()+
  theme(text = element_text(size = 20))+
  xlab("class")+
  scale_x_discrete(labels=c("Non-linked PEG-MEG", "Linked PEG-MEG"))+
  geom_signif(comparisons = list(c("0","1")),test = "wilcox.test",map_signif_level = "FALSE")






### HIC data analysis
### hic data analysis
# 1.1 find interacting impritned gene paris by katyoplote
library(GenomicInteractions)
require(GenomicFeatures)
require(karyoploteR)
tr<-makeTxDbFromGFF("~/MyResearch/genome.db/TAIR/gtf/Arabidopsis_thaliana.TAIR10.31.gtf",format = "gtf")
genes<-genes(tr)

GenomeInteract<-function(interaction_file="other_resources/HIC/sigInteractions_10kb.csv"){
  ara_intra_hic<-read.csv(interaction_file)
  region1<-GRanges(seqnames = ara_intra_hic$chrom_left,ranges = IRanges(start = ara_intra_hic$start_left,end = ara_intra_hic$end_left))
  region2<-GRanges(seqnames = ara_intra_hic$chrom_right,ranges = IRanges(start = ara_intra_hic$start_right,end = ara_intra_hic$end_right))
  gi <- GInteractions(region1, region2)
  gi$reads<-ara_intra_hic$Reads.counts
  gi$p_value<-ara_intra_hic$p.value
  return(gi)
}
GI_intra<-GenomeInteract(interaction_file="other_resources/HIC/sigInteractions_10kb.csv") 
GI_inter <-GenomeInteract(interaction_file = "other_resources/HIC/sigInterChrInteractions_10kb.csv")
GI_all<- c(GI_intra,GI_inter) # all hic interaction 

Interaction_analysis<-function(gene_list=imprinted_genes,Hic_interacts=GI_all){
  IGs_ranges<-genes(tr)[gene_list]
  
  ##### genes overlaps with interacting paris
  overlaps1<-findOverlaps(query = IGs_ranges,subject = Hic_interacts,use.region="first") # left overlap
  message(c("Left overlapped gene number:",length(unique(queryHits(overlaps1)))))
  
  overlaps2<-findOverlaps(query = IGs_ranges,subject = Hic_interacts,use.region="second") # right overlaps
  message(c("Right overlapped gene number:",length(unique(queryHits(overlaps2)))))
  ## find target interaction pairs
  interac_freq<-sort(table(c(unique(subjectHits(overlaps1)),unique(subjectHits(overlaps2)))),decreasing = TRUE) 
  interac_freq<-interac_freq[interac_freq>1]
  message(c("Interaction gene pairs:",length(interac_freq)))
  
  ### visuliaztion (mapping interacting genes into interaction pairs)
  ara_pos<-read.table("/home/wuzefeng/MyResearch/genome.db/TAIR/dna/chrom_size")
  ara_pos$start=1
  colnames(ara_pos)<-c("chr","end","start")
  ara_pos<-ara_pos[,c(1,3,2)]
  ara_pos<-ara_pos[order(ara_pos$chr),]
  kp <- plotKaryotype(genome = toGRanges(ara_pos),plot.type = 1)
  
  ## get interaction genes
  for (name in names(interac_freq)){
    target_IG_left <- IGs_ranges[queryHits(overlaps1)[subjectHits(overlaps1)==name]]
    target_IG_right <- IGs_ranges[queryHits(overlaps2)[subjectHits(overlaps2)==name]]
    
    genes_interacts<-data.frame(left = target_IG_left,right = target_IG_right)
    colnames(genes_interacts)[7:9]<-c("link.chr","link.start","link.end")
    #print(genes_interacts)
    genes_interacts<-toGRanges(genes_interacts[,c(1,2,3,7,8,9,11,6,12)])
    strand(genes_interacts)<-strand(target_IG_left)
    
    kpPlotLinks(kp,genes_interacts,arch.height = 0.8,col="blue")
    kpText(kp,data = target_IG_left,labels = target_IG_left$gene_id,y = seq(target_IG_left$gene_id)*0.1,col="tomato",cex=0.8)
    kpText(kp,data = target_IG_right,labels = target_IG_right$gene_id,y = -(seq(target_IG_right$gene_id))*0.5,col="tomato",cex=0.8)
  }
  return(length(interac_freq))
}
Interaction_analysis(gene_list = imprinted_genes,Hic_interacts=GI_all)

Interaction_analysis_extend_promoter<-function(gene_list=imprinted_genes,Hic_interacts=GI_all){
  IGs_GR<-genes(tr)[gene_list]
  IGs_ranges<-punion(promoters(IGs_GR,upstream = 1000),IGs_GR)
  IGs_ranges$gene_id<-names(IGs_ranges)
  
  ##### genes overlaps with interacting paris
  overlaps1<-findOverlaps(query = IGs_ranges,subject = Hic_interacts,use.region="first")
  message(c("Left overlapped gene number:",length(unique(queryHits(overlaps1)))))
  
  overlaps2<-findOverlaps(query = IGs_ranges,subject = Hic_interacts,use.region="second") 
  message(c("Right overlapped gene number:",length(unique(queryHits(overlaps2)))))
  ## find target interaction pairs
  interac_freq<-sort(table(c(unique(subjectHits(overlaps1)),unique(subjectHits(overlaps2)))),decreasing = TRUE) 
  interac_freq<-interac_freq[interac_freq>1]
  message(c("Interaction gene pairs:",length(interac_freq)))
  
  ### visuliaztion (mapping interacting genes into interaction pairs)
  ara_pos<-read.table("/home/wuzefeng/MyResearch/genome.db/TAIR/dna/chrom_size")
  ara_pos$start=1
  colnames(ara_pos)<-c("chr","end","start")
  ara_pos<-ara_pos[,c(1,3,2)]
  ara_pos<-ara_pos[order(ara_pos$chr),]
  kp <- plotKaryotype(genome = toGRanges(ara_pos),plot.type = 1)
  
  interac_pairs<-list()
  ## get interaction genes
  for (name in names(interac_freq)){
    message(name)
    target_IG_left <- IGs_ranges[queryHits(overlaps1)[subjectHits(overlaps1)==name]]
    target_IG_right <- IGs_ranges[queryHits(overlaps2)[subjectHits(overlaps2)==name]]
    
    genes_interacts<-data.frame(left = target_IG_left,right = target_IG_right)
    interac_pairs[[name]]<-genes_interacts[,c(6,12)]
    colnames(genes_interacts)[7:9]<-c("link.chr","link.start","link.end")
    genes_interacts<-toGRanges(genes_interacts[,c(1,2,3,7,8,9,11,6,12)])
    strand(genes_interacts)<-strand(target_IG_left)
    
    kpPlotLinks(kp,genes_interacts,arch.height = 0.8,col="blue")
    kpText(kp,data = target_IG_left,labels = target_IG_left$gene_id,y = seq(target_IG_left$gene_id)*0.1,col="tomato",cex=0.8)
    kpText(kp,data = target_IG_right,labels = target_IG_right$gene_id,y = -(seq(target_IG_right$gene_id))*0.5,col="tomato",cex=0.8)
  }
  return(interac_pairs)
}
Interaction_analysis_extend_promoter(gene_list = imprinted_genes,Hic_interacts=GI_all)


#simulation to interacting pair number
Interaction_simulation_analysis<-function(gene_list=imprinted_genes,Hic_interacts=GI_all){
  IGs_ranges<-genes(tr)[gene_list]
  
  ##### genes overlaps with interacting paris
  overlaps1<-findOverlaps(query = IGs_ranges,subject = Hic_interacts,use.region="first")
  message(c("Left overlapped gene number:",length(unique(queryHits(overlaps1)))))
  
  overlaps2<-findOverlaps(query = IGs_ranges,subject = Hic_interacts,use.region="second") 
  message(c("Right overlapped gene number:",length(unique(queryHits(overlaps2)))))
  ## find target interaction pairs
  interac_freq<-sort(table(c(unique(subjectHits(overlaps1)),unique(subjectHits(overlaps2)))),decreasing = TRUE) 
  interac_freq<-interac_freq[interac_freq>1]
  message(c("Interaction gene pairs:",length(interac_freq)))
  return(length(interac_freq))
}
simulations_genes_pairs<-c()
for(m in 1:1000){
  message(m)
  simulations_genes_pairs<-c(simulations_genes_pairs,Interaction_simulation_analysis(gene_list = sample(names(genes(tr)),length(imprinted_genes)),Hic_interacts=GI_all))
}

ggplot(data.frame(genes_pairs=simulations_genes_pairs,stringsAsFactors = FALSE),aes(x=genes_pairs))+
  geom_histogram(bins = 100)+
  theme_bw()+
  theme(text = element_text(size = 20))+
  annotate(geom = "segment",x = 8,xend = 8,y = 25,yend = 0,arrow=arrow(),color="red",size=2)+
  annotate(geom="text",label="Mean", x = 8, y = 30,size=5)+
  ylab("Frequency")+xlab("Interacting gene pairs")



## 1.2 (option) find interacting imprinted gene paris by circulize (circos)
ara_pos<-read.table("/home/wuzefeng/MyResearch/genome.db/TAIR/dna/chrom_size")
ara_pos$start=1
colnames(ara_pos)<-c("chr","end","start")
ara_pos<-ara_pos[,c(1,3,2)]
ara_pos<-ara_pos[order(ara_pos$chr),]

#####(i) initiate circos (order aware !!! (0:inner--->1:outer))
library(circlize)
circos.initializeWithIdeogram(ara_pos,plotType = NULL)

##### epxression 
## expression 
library(stringr)
IG_bed = genes[gene_label_list]
IG_bed<-as.data.frame(IG_bed)[,c(1,2,3)]

### (i) expression heatmap
tissues<-unique(str_sub(colnames(rna_expression),start = 1,end = -2))
tissue_gene_expression<-list()
for (tissue in tissues){
  tissue_gene_expression[[tissue]]<-rowSums(rna_expression[,grep(tissue,colnames(rna_expression))])/3
}
gene_expression_bed<-do.call(cbind.data.frame, tissue_gene_expression)
gene_expression_bed_IGs<-gene_expression_bed[rownames(IG_bed),]

library(BBmisc)
gene_expression_bed_IGs<-t(apply(gene_expression_bed_IGs,1,function(x) normalize(x,range=c(-1,1),method="range")))
expression_IG_bed <- cbind(IG_bed,gene_expression_bed_IGs)


col_fun = colorRamp2(c(-1, 0, 1), c("RED", "black", "green"))
circos.genomicHeatmap(expression_IG_bed, col = col_fun, side = "outside", border = "white")

#######(ii) get interaction paris (bed file and corss links)
Interaction_analysis<-function(gene_list=imprinted_genes,Hic_interacts=GI_all){
  IGs_ranges<-genes(tr)[gene_list]
  
  ##### genes overlaps with interacting paris
  overlaps1<-findOverlaps(query = IGs_ranges,subject = Hic_interacts,use.region="first") # left overlap
  message(c("Left overlapped gene number:",length(unique(queryHits(overlaps1)))))
  
  overlaps2<-findOverlaps(query = IGs_ranges,subject = Hic_interacts,use.region="second") # right overlaps
  message(c("Right overlapped gene number:",length(unique(queryHits(overlaps2)))))
  ## find target interaction pairs
  interac_freq<-sort(table(c(unique(subjectHits(overlaps1)),unique(subjectHits(overlaps2)))),decreasing = TRUE) 
  interac_freq<-interac_freq[interac_freq>1]
  message(c("Interaction gene pairs:",length(interac_freq)))
  
  
  label_bed_datalist = list()
  
  ## get each interaction genes paris or links
  for (name in names(interac_freq)){
    target_IG_left <- IGs_ranges[queryHits(overlaps1)[subjectHits(overlaps1)==name]]
    target_IG_right <- IGs_ranges[queryHits(overlaps2)[subjectHits(overlaps2)==name]]
    
    ## interactinf imprinted genes paris to two bed files
    genes_interacts<-data.frame(left = target_IG_left,right = target_IG_right) # data frame is usful when one to multiple
    bed_left <- genes_interacts[,1:6]
    bed_right <- genes_interacts[,7:12]
    
    
    ## add links to circos
    colnames(bed_left)<-colnames(bed_right)<-c("chr","start","end","value","strand","gene_ID")
    circos.genomicLink(bed_left, bed_right, col = rand_color(nrow(bed_left), transparency = 0), border = "black",rou = 0.3)
    
    ## add label to bed file. (not loop to add label,because for each add, it will add a track)
    label_bed_left<-unique(bed_left)
    label_bed_right<-unique(bed_right)
    label_bed_datalist[[name]]<-rbind(label_bed_left,label_bed_right)
  
  }
  label_bed = do.call(rbind, label_bed_datalist)
  #print(head(label_bed))
  circos.genomicLabels(bed = unique(label_bed),labels.column = 6,side = "outside",col=ifelse(unique(label_bed)$gene_ID%in%maternal_imprint,"tomato","blue"),labels_height = 0.19)
  return(unique(label_bed$gene_ID))
}
gene_label_list<-Interaction_analysis()

##(iii) chromosomes
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  #circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(7))
  circos.text(mean(xlim), mean(ylim), chr, cex = 1, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA, bg.col=c("#8B6508", "#DEB887", "#458B74", "#483D8B", "#104E8B", "#CD9B1D", "#6495ED")
)
### gene density
circos.genomicDensity(as.data.frame(genes), window.size = 1e6,track.height = 0.15,bg.border="black")




#1.2
### hic tf-->impritned genes
tf_hic<-read.table("~/Desktop/imprinted_TF.txt",sep=",",stringsAsFactors = FALSE,header = TRUE)

library(gplots)
library(VennDiagram)
venn.diagram(
  x = list(
    IFGN = tf_hic$neigs,
    HiC = rownames(mmm)
  ),
  imagetype = "tiff",
  filename = "./3D-quadruple_Venn.tiff",#-outfile-name--#
  col = "black",                      #----border_colour#
  lty='dotted',
  lwd = 4,
  fill = c("cornflowerblue", "green"),#filling colour,length=length(x)#
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4"),#digit colour#
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen"),
  cat.cex = 1.5,
  #cat.dist = c(0.25,0.25),
  cat.fontfamily = "serif",
  cat.dist = c(0.051,0.051),
  main="TFs associated with imprinted genes",
  main.pos= c(0.5,0.53),
  margin=4)



##1.3  DNase-seq footprint and fimo TF-target gene network and hic enhancer 
#### hiv plot for enhancer
promoter_tfs_of_dnase<-read.csv("other_resources/LZH_dnase/impGenePromoterTF.csv")
promoter_tfs_of_dnase$type<-"promoter"
colnames(promoter_tfs_of_dnase)<-c("a","b","type")

enhancer_tfs_of_dnase<-read.csv("other_resources/LZH_dnase/impGeneEnhancerTF.csv")
enhancer_tfs_of_dnase<-enhancer_tfs_of_dnase[,c(2,1)]
enhancer_tfs_of_dnase$type<-"enhancer"
colnames(enhancer_tfs_of_dnase)<-c("a","b","type")
df<-rbind(promoter_tfs_of_dnase,enhancer_tfs_of_dnase)

library(ggraph)
pe_tf <- graph_from_data_frame(df)
#ggraph(pe_tf, layout = 'kk',circular=TRUE) + geom_edge_link(aes(colour = type))

###plot by igraph
par(bg="white",mar=c(2,2,2,2))

V(pe_tf)$color[names(V(pe_tf))%in%promoter_tfs_of_dnase$a]<-"tomato" # vertex color 
V(pe_tf)$color[names(V(pe_tf))%in%promoter_tfs_of_dnase$b]<-"steelblue"

V(pe_tf)$size <- 8              # vertex size
V(pe_tf)$label.cex <- 0.8         # vertex label size
V(pe_tf)$label.color<-"black"
E(pe_tf)$color<-ifelse(E(pe_tf)$type=="promoter","orange","gray")        # edge color 

E(pe_tf)$width=3
plot(pe_tf,layout=layout.fruchterman.reingold,
     vertex.frame.color= NA,
     vertex.color=V(pe_tf)$color,
     vetex.shape=V(pe_tf)$shape)

legend('topleft',
       #bg="white",
       text.col="black",
       legend=c("Paternal","Maternal","Correlated","Anti-correlated"),
       pch=c(19,19,NA,NA), #shape
       lty = c(0.5,0.5,1,1),
       box.lty = 2, # 
       pt.cex = 3, #lines size 
       cex=1, #box size
       col=c("tomato","steelblue","orange","gray"),
       y.intersp=1.5)


### hive plot
## hive plot for TFs in enhancer, promoter and imprinted genes
library(HiveR)
promoter_tfs_of_dnase<-read.csv("other_resources/LZH_dnase/impGenePromoterTF.csv")
promoter_tfs_of_dnase$promoterTFID<-paste("p",promoter_tfs_of_dnase$promoterTFID,sep = "")
promoter_tfs_of_dnase$weight<-1
colnames(promoter_tfs_of_dnase)<-c("tf","target","weight")

enhancer_tfs_of_dnase<-read.csv("other_resources/LZH_dnase/impGeneEnhancerTF.csv")
enhancer_tfs_of_dnase<-enhancer_tfs_of_dnase[,c(2,1)]
enhancer_tfs_of_dnase$EnhancerTFID<-paste("e",enhancer_tfs_of_dnase$EnhancerTFID,sep="")
enhancer_tfs_of_dnase$weight<-1
colnames(enhancer_tfs_of_dnase)<-c("tf","target","weight")

df<-rbind(promoter_tfs_of_dnase,enhancer_tfs_of_dnase)
hiv<-edge2HPD(df)
hiv$nodes$axis[startsWith(hiv$nodes$lab,"p")]<-2
hiv$nodes$axis[startsWith(hiv$nodes$lab,"e")]<-3
hiv$nodes$axis<-as.integer(hiv$nodes$axis)

## node position or radium
#hiv$nodes$radius<-as.numeric(sample(seq(1:1000),length(hiv$nodes$radius),replace = TRUE))
hiv$nodes$radius<-as.numeric(c(seq(1:length(hiv$nodes$axis[startsWith(hiv$nodes$lab,"p")])),
                            seq(1:length(hiv$nodes$axis[startsWith(hiv$nodes$lab,"e")])),
                            seq(1:length(hiv$nodes$axis[startsWith(hiv$nodes$lab,"AT")]))))

# node color
#hiv$nodes$color[startsWith(hiv$nodes$lab,"p")]<-"red"
#hiv$nodes$color[startsWith(hiv$nodes$lab,"e")]<-"blue"
#hiv$nodes$color[!startsWith(hiv$nodes$lab,"p")&!startsWith(hiv$nodes$lab,"e")]<-"yellow"

## axis color
hiv$axis.cols<-"black" #c("#E41A1C", "#377EB8", "#4DAF4A")
## plot
plotHive(hiv,axLabs = c("Imprinted genes", "TFs in promoter","TFs in enhancer"),
         ch = 5,
         dr.nodes = FALSE, # whether show the node
         bkgnd = "white",
         axLab.gpar=gpar(col="black"),
         axLab.pos = c(20, 80, 80),
         rot = c(0, 25, -25))





set.seed(1000)
options(digits = 5)
##import data
data1 <-read.table("3prediction_result/tissue_naive_prediction.txt",header = TRUE,stringsAsFactors = FALSE)
data2 <-read.table("3prediction_result/tissue_naive_prediction2.txt",header = FALSE,stringsAsFactors = FALSE)
colnames(data2)<-colnames(data1)
data<-rbind(data1,data2)
# data$pos<-log(data$pos/data$neg,2)
data<-data[,c(1,2,4)]
colnames(data)<-c("X1","X2","weight")
data<-data[!duplicated(data[c(1,2)]),] ## exsit duplciations and remove


### simplify network
library(igraph)
g<-graph.data.frame(d = data, directed = FALSE)
g<-simplify(g,remove.multiple = TRUE,remove.loops = TRUE) #22248 nodes * 578573 edges


### network validation using tair interaction data
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
  write.table(file="out",data.frame(gene=x,experiment=dim(focus_links)[1],predicted=length(intersect(names(V(sub_g)),unique(as.vector(t(as.matrix(focus_links))))))-1),sep="\t",append=TRUE,col.names = !file.exists("out"),row.names = FALSE,quote = FALSE)
  focus_links<-subset(focus_links,focus_links$V2%in%names(V(sub_g))&focus_links$V1%in%names(V(sub_g))) 
  }
} 


#### individule target gene validation and showing  and plot
tair_inter<-read.table("/home/wuzefeng/MyResearch/networks/PPI/Arabidopsis/TairProteinInteraction.20090527.txt",header = TRUE,stringsAsFactors = FALSE,sep="\t")
focused_genes<-"AT1G02340" # AT1G02340 #AT2G18790
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
  E(sub_g)$width=1+E(sub_g)$weight
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



## kegg validation
### network validation using kegg pathway
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


#### reanalysis kegg data with permutations network
permute_g<-g
V(permute_g)$name<-sample(names(V(g)))

for (col in 1:ncol(combins)){
  gene_list1<-xx[[combins[,col][1]]]
  gene_list2<-xx[[combins[,col][2]]]
  
  group_in_possible<-choose(length(gene_list1),2)+choose(length(gene_list2),2)-choose(length(intersect(gene_list1,gene_list2)),2)
  group_in_prediction <- ecount(induced_subgraph(permute_g,intersect(gene_list1,names(V(permute_g)))) %u% induced_subgraph(permute_g,intersect(gene_list2,names(V(permute_g)))))
  
  between_group_possible <-length(gene_list1)*length(gene_list2)-length(intersect(gene_list1,gene_list2))^2+choose(length(intersect(gene_list1,gene_list2)),2)
  between_group_prediction <-ecount(induced_subgraph(permute_g,intersect(unique(c(gene_list1,gene_list2)),names(V(permute_g)))))-group_in_prediction
  
  print(c(col, group_in_possible,group_in_prediction,between_group_possible,between_group_prediction))
  message(fisher.test(matrix(c(group_in_prediction,group_in_possible,between_group_prediction,between_group_possible),nrow=2),alternative="greater")$p.value)
  write.table(file="KEGG_validation_random_network",data.frame(pathway=paste(names(xx)[combins[,col][1]],names(xx)[combins[,col][2]],sep = "_"),
                                                predicted_in_pathway = group_in_prediction,
                                                in_pathway = group_in_possible,
                                                prediction_between_pathway = between_group_prediction,
                                                between_pathway = between_group_possible,
                                                p_value = fisher.test(matrix(c(group_in_prediction,group_in_possible,between_group_prediction,between_group_possible),nrow=2),alternative="greater")$p.value),
                                                sep="\t",append=TRUE,col.names = !file.exists("KEGG_validation_random_network"),
                                                row.names = FALSE,quote = FALSE)
}


##### indivudule kegg pathway validation 
require(pathview)
library(org.At.tair.db)
xx<- as.list(org.At.tairPATH2TAIR)
focused_kegg_Pathway<-"04626"
sub_g<-induced_subgraph(g,vids = intersect(xx[[focused_kegg_Pathway]],names(V(g))))
ath.dat.kegg <- sim.mol.data(mol.type="gene",id.type="tair",species="ath",nmol=3000)
pv.out <- pathview(gene.data = ath.dat.kegg, gene.idtype="tair",pathway.id = focused_kegg_Pathway, species = "ath", out.suffix = "ath.kegg",kegg.native = T, same.layer=T)


#### compare wit aranet2
ara_net2 <-read.table("4network_analysis/AraNet.txt",header = FALSE,stringsAsFactors = FALSE)
colnames(ara_net2)<-c("X1","X2","weight")
ara_net2<-ara_net2[!duplicated(ara_net2[c(1,2)]),]
library(igraph)
ara_net_g<-graph.data.frame(d = ara_net2, directed = FALSE)
ara_net_g<-simplify(ara_net_g,remove.multiple = TRUE,remove.loops = TRUE) #22894 *895000 

intersect_g<-ara_net_g%s%g #91987
intersect_permute_g<-permute_g%s%ara_net_g #1957


### analysis flowering geneset
##nalysis flowering geneset

flowering_genes<-read.csv("~/Desktop/Flowering Interactive Database FLOR-ID - Flowering time.csv",header = TRUE,stringsAsFactors = FALSE)
commen_flower_genes<- intersect(flowering_genes$Gene.details,names(V(g)))

#short distence among these flowering genes
d_flower<-shortest.paths(g,commen_flower_genes, commen_flower_genes,weights = NA)
hist(d)
#random smapeld genes
rand_genes<-sample(names(V(g)),286,replace = FALSE)
d_random<-shortest.paths(g,rand_genes, rand_genes,weights = NA)

## distence among these genes
d_flower<-distances(g,commen_flower_genes, commen_flower_genes,weights = NA)
d_random<-distances(g,rand_genes, rand_genes,weights = NA)


### essential genes analysis
### essential genes analysis

essential_genes<-read.table("other_sourse/essential_genes.csv",stringsAsFactors = FALSE)$V1
essential_genes<-intersect(essential_genes,names(V(g)))

betweenness_essential <- betweenness(g,v = essential_genes,directed = FALSE,normalized = TRUE) ## long time
betweenness_random <- betweenness(g,v = sample(names(V(g)),length(essential_genes),replace = FALSE),directed = FALSE,normalized = TRUE) # long time


## sub graph for imprinted genes
### import imprinted genes
imprinted_genes <-read.table("/home/wuzefeng/MyResearch/Imprinting_prediction/imprint_gene_list/TAIR10/ara_imp_by_paper_num/imp2+.list",stringsAsFactors = FALSE)
imprinted_genes<-unique(imprinted_genes$V1)
all_neibours <-unique(unlist(lapply(imprinted_genes,function(x) if (x %in% names(V(g))){names(neighbors(g,x))})))
sub_g<-induced_subgraph(g,all_neibours)




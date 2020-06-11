#!/bin/R
trainning_set_unsort <-read.table("/home/wuzefeng/MyResearch/networks/GO_annotation/arabidopsis/1GO_sim_18.3.15_out/4gold_standard.gene_pairs.txt")  ### training gene pair contain positive and negitive pairs 

trainning_set<-data.frame(t(apply(trainning_set_unsort[c(-3,-4)],1,sort)))
trainning_set$weight<-trainning_set_unsort$V3# sort two columns   # X1,X2,X3
trainning_set$class<-trainning_set_unsort$V4   # 2196935 instances
ind <- duplicated(trainning_set[,1:2])
trainning_set<-trainning_set[!ind,]  # 1796174 neg +196491 pos (very important) # dim(1992655,4)

#trainning_set <-trainning_set[sample(seq(1:nrow(trainning_set)),100),] # sampling
rm(trainning_set_unsort)
gc()

## 
library(ggplot2)
dd<-rbind(trainning_set[1:100,], tail(trainning_set,100))
p<-ggplot(trainning_set,aes(x=weight,fill=class))+geom_density(alpha=0.5)+theme(text=element_text(size=20))



fisher.z<-function(m){z <- 0.5*(log(1+m+0.0001)-log(1-m+0.0001));return(z)}  # fisher.z transform 

######################### function to prepare data from different vidence sources

Data_Preprocess<-function(file_name){    # rows are genes and group information  are columns elements (quick)
  data<-read.table(file_name,header = TRUE)
  data<-t(data)
  
  uniquelength <- apply(data,2,function(x) length(unique(x))) # drop genes columns with same values (can not get the corrrleation)
  data <- subset(data, select=uniquelength>1)
  
  gene_list<-vector(mode="list",length=ncol(data))  # convert gene names in matrix to numbers 
  names(gene_list)<-colnames(data)
  colnames(data)<-seq(1:ncol(data))
  
  genelist_data<-list(gene_list,data)  ########### return matrix and gene list
  return(genelist_data)
} 

Matrix_Transform<-function(data){  # colnames are gene names    ? problems
  cor_matrix<-cor(data)
  #cor_matrix<- fisher.z(cor_matrix)
  #mean<-mean(cor_matrix[upper.tri(cor_matrix,diag = FALSE)]) #mean vlaue
  #cor_matrix <- (cor_matrix-mean)/sd(cor_matrix[upper.tri(cor_matrix,diag = FALSE)]) #mean substrated and divided by sd
  #cor_matrix <- abs(cor_matrix)
  return(cor_matrix)
}

################################## gene distance function 
library(GenomicFeatures)
tx<-makeTxDbFromGFF("/home/wuzefeng/MyResearch/genome.db/TAIR/gtf/Arabidopsis_thaliana.TAIR10.31.gtf")
genes<-genes(tx)
genes$number<-1:length(genes)

Gene_Distance_By_number <- function(gene1,gene2){
  if (!gene1%in% names(genes)|!gene2%in% names(genes))
    {return (NA)}
  chro_info1<-substr(gene1,1,4)
  chro_info2<-substr(gene2,1,4)
  if (chro_info2==chro_info1){
    gene_dis = 1/abs(genes[gene1]$number-genes[gene2]$number)
    return(gene_dis)
  }
  else return(0)
}

########## gene distance
trainning_set$distance<-apply(trainning_set,1,function(x) {Gene_Distance_By_number(x[1],x[2])}) # not significant (too slow)
p<-ggplot(trainning_set,aes(y=distance_log,x=class))+geom_boxplot()+xlab("Gene class")+ylab("Distance of gene pair(-log2)")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
###########  add weights (from gene expression data correation,TF corr, etal ) for each pair of genes in training set(f1)
express_data<-Data_Preprocess("/home/wuzefeng/MyResearch/RNA_Seq/Plants/Arabidopsis/6cluster_out/cluster_82.exp") #expression data
express_gene_list <-express_data[[1]]
express_cor<-Matrix_Transform(express_data[[2]])  #slow (5min)
trainning_set$pcc<-apply(trainning_set,1,function(x) express_cor[match(x[1],names(express_gene_list)),match(x[2],names(express_gene_list))]) 
#p<-ggplot(na.omit(trainning_set),aes(y=pcc,x=class))+geom_boxplot()+xlab("")+theme(text=element_text(size=20))+ylab("z-score")
#library(ggsignif)
#p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE,step_increase = 0.1)

trainning_set$pcc<-cut(trainning_set$pcc,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])  # discreation data
rm(express_cor)
gc()

tf_bind_data<-Data_Preprocess("/home/wuzefeng/MyResearch/motif_dbs/5jaspar2018/1fimo_out/Tair2_p_0.01/gene2motif_0_1.txt")  # tf binding # updataed to p-0.001, and jaspar2018
tf_gene_list<-tf_bind_data[[1]]
tf_cor<-Matrix_Transform(tf_bind_data[[2]]) #slow
trainning_set$tf<-apply(trainning_set,1,function(x) tf_cor[match(x[1],names(tf_gene_list)),match(x[2],names(tf_gene_list))]) #slow
p<-ggplot(na.omit(trainning_set),aes(y=tf,x=class))+geom_boxplot()+xlab("")+theme(text=element_text(size=20))+ylab("z-score")
#library(ggsignif)
#p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE,step_increase = 0.1)
trainning_set$tf<-cut(trainning_set$tf,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
rm(tf_cor)
gc()

###arrary data(optional)
arrary_data<-Data_Preprocess("/home/wuzefeng/MyResearch/microarray/AT40_RMA_MAS/5ara_express") ###array data
arrary_gene_list <-arrary_data[[1]]
arrary_cor<-Matrix_Transform(arrary_data[[2]])
trainning_set$microarrary<-apply(trainning_set,1,function(x) arrary_cor[match(x[1],names(arrary_gene_list)),match(x[2],names(arrary_gene_list))]) 
rm(arrary_cor)
gc()
#
source("/home/wuzefeng/R/ggplot_publication_theme.R")
#p<-ggplot(trainning_set,aes(x=V6,fill=V3))+geom_density(alpha=0.1)
trainning_set$V6<-cut(trainning_set$V6,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])

############## import PPI-data

data_biogrid<-read.table("/home/wuzefeng/MyResearch/networks/PPI/Arabidopsis/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-3.4.134.tab2.txt",header = TRUE,sep = "\t", check.names=FALSE)
data_biogrid<-data.frame(data_biogrid$`Systematic Name Interactor A`,data_biogrid$`Systematic Name Interactor B`)
colnames(data_biogrid)<-c("X1","X2")
indx <- !duplicated(t(apply(data_biogrid, 1, sort))) # finds non - duplicates in sorted rows
data_biogrid<-data_biogrid[indx, ]
data_biogrid$X1<-as.character(data_biogrid$X1)
data_biogrid$X2<-as.character(data_biogrid$X2)
data_biogrid<-subset(data_biogrid,data_biogrid$X1!=data_biogrid$X2)
data_biogrid<-as.data.frame(t(apply(data_biogrid,1,sort)))
colnames(data_biogrid)<-c("X1","X2")
data_biogrid$biogrid<-1
trainning_set<-merge(trainning_set,data_biogrid,by = c("X1","X2"),all.x = TRUE)
p<-ggplot(na.omit(trainning_set),aes(x=as.factor(biogrid),fill=class))+geom_bar(position = "dodge")+xlab("")+theme(text=element_text(size = 20))
######################### import phylogentics files

phy1<-Data_Preprocess("/home/wuzefeng/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair.ensembl_bacterial.gene2species")
phy2<-Data_Preprocess("/home/wuzefeng/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair_ensembl.fungi.gene2species")
phy3<-Data_Preprocess("/home/wuzefeng/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair.ensembl.gene2species")
phy4<-Data_Preprocess("/home/wuzefeng/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair_metazoa.gene2species")
phy5<-Data_Preprocess("/home/wuzefeng/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair_phytozome.gene2species")

profile_gene_list <-phy1[[1]]
profile_cor<-Matrix_Transform(phy1[[2]])
trainning_set$profile_phy1<-apply(trainning_set,1,function(x) profile_cor[match(x[1],names(profile_gene_list)),match(x[2],names(profile_gene_list))]) # exist na
trainning_set$profile_phy1<-cut(trainning_set$profile_phy1,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
rm(profile_cor)
gc()

p<-ggplot(na.omit(trainning_set),aes(y=profile_phy1,x=class))+geom_boxplot()+xlab("")+theme(text=element_text(size=20))+ylab("Phylogenetic profile similarity (bacterial)")
#library(ggsignif)
#p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE,step_increase = 0.1)



profile_gene_list <-phy2[[1]]
profile_cor<-Matrix_Transform(phy2[[2]])
trainning_set$profile_phy2<-apply(trainning_set,1,function(x) profile_cor[match(x[1],names(profile_gene_list)),match(x[2],names(profile_gene_list))]) 
#trainning_set$profile_phy2<-cut(trainning_set$profile_phy2,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
p<-ggplot(na.omit(trainning_set),aes(y=profile_phy2,x=class))+geom_boxplot()+xlab("")+theme(text=element_text(size=20))+ylab("Phylogenetic profile similarity (fungi)")
rm(profile_cor)
gc()

profile_gene_list <-phy3[[1]] 
profile_cor<-phy3[[2]]
profile_cor<-Matrix_Transform(profile_cor)
trainning_set$profile_phy3<-apply(trainning_set,1,function(x) profile_cor[match(x[1],names(profile_gene_list)),match(x[2],names(profile_gene_list))]) 
#trainning_set$profile_phy3<-cut(trainning_set$profile_phy3,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
#p<-ggplot(na.omit(trainning_set),aes(y=profile_phy3,x=class))+geom_boxplot()+xlab("")+theme(text=element_text(size=20))+ylab("Phylogenetic profile similarity (ensembl)")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE,step_increase = 0.1)
rm(profile_cor)
gc()


profile_gene_list <-phy4[[1]]
profile_cor<-phy4[[2]]
profile_cor<-Matrix_Transform(profile_cor)
trainning_set$profile_phy4<-apply(trainning_set,1,function(x) profile_cor[match(x[1],names(profile_gene_list)),match(x[2],names(profile_gene_list))]) 
#trainning_set$profile_phy4<-cut(trainning_set$profile_phy4,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
#p<-ggplot(na.omit(trainning_set),aes(y=profile_phy4,x=class))+geom_boxplot()+xlab("")+theme(text=element_text(size=20))+ylab("Phylogenetic profile similarity(metazoa)")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE,step_increase = 0.1)
rm(profile_cor)
gc()

profile_gene_list <-phy5[[1]]
profile_cor <- phy5[[2]]
profile_cor<-Matrix_Transform(profile_cor)
trainning_set$profile_phy5<-apply(trainning_set,1,function(x) profile_cor[match(x[1],names(profile_gene_list)),match(x[2],names(profile_gene_list))]) 
#trainning_set$profile_phy5<-cut(trainning_set$profile_phy5,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
p<-ggplot(na.omit(trainning_set),aes(y=profile_phy5,x=class))+geom_boxplot()+xlab("")+theme(text=element_text(size=20))+ylab("Phylogenetic profile similarity (phytozome)")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE,step_increase = 0.1)
rm(profile_cor)
gc()



################## reactome (database)interaction
rec<-read.table("reactome/3pair.int.txt")
rec<-as.data.frame(t(apply(rec,1,sort)))
colnames(rec)<-c("X1","X2")
rec$reactome<-1
rec<-unique(rec)
trainning_set<-merge(trainning_set,rec,by = c("X1","X2"),all.x = TRUE)
p<-ggplot(subset(trainning_set,trainning_set$reactome==1),aes(x=as.factor(reactome),fill=class))+geom_bar(position = "dodge")+theme_Publication()+xlab("")

################# domain interaction
gene_domain_intact<-read.table("/home/wuzefeng/MyResearch/networks/UniDomInt/ara/4ara_int.txt",header = FALSE)
score<-gene_domain_intact$V3
#gene_domain_intact<-data.frame(t(apply(gene_domain_intact[-3],1,sort))) # to further match 
gene_domain_intact$domain<-score
gene_domain_intact<-gene_domain_intact[-4]
colnames(gene_domain_intact)<-c("X1","X2","domain")
gene_domain_intact<-as.data.frame(gene_domain_intact %>% group_by(X1, X2) %>% summarise_all(funs(mean)))
#write.table("UniDomInt/ara/3ara_int.order.txt",col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
trainning_set<-merge(trainning_set,gene_domain_intact,by = c("X1","X2"),all.x = TRUE)
#p<-ggplot(na.omit(trainning_set),aes(y=domain,x=class))+geom_boxplot()+xlab("")+theme(text=element_text(size=20))+ylab("Domain similarities")
#p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE,step_increase = 0.1)+ylab("Domain interaction score")
trainning_set$domain<-cut(trainning_set$domain,breaks = c(0,0.2,0.4,0.6,0.8,1),labels = letters[1:5])

trainning_set_orange<-data.frame(item=paste(trainning_set$X1,trainning_set$X2,sep = "_"),trainning_set$V4,trainning_set$V5,trainning_set$V6,trainning_set$biogrid,trainning_set$profile_phy1,trainning_set$profile_phy2,trainning_set$profile_phy3,trainning_set$profile_phy4,trainning_set$profile_phy5,trainning_set$cor,trainning_set$score,trainning_set$class)
write.table(trainning_set_orange,"/home/wuzefeng/MyResearch/networks/1feature_plot/2018.4.16/train_data2.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)


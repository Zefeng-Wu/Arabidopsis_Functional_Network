#!/bin/R

######################### function to prepare data from different vidence sources

Data_Preprocess<-function(file_name){    # rows are genes and group information are columns.
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
  return(cor_matrix)
}

######################### function to prepare data from different sources

Gene_RNA_seq_pairs<-function(data){  # use rna-seq pairs as "training set"
  cor_matrix<-cor(t(data))
  return(cor_matrix)
}


### other source of data

## feature1. biogrid data 

data_biogrid<-read.table("~/MyResearch/networks/PPI/Arabidopsis/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-3.4.134.tab2.txt",header = TRUE,sep = "\t", check.names=FALSE,stringsAsFactors = FALSE)
data_biogrid<-data.frame(data_biogrid$`Systematic Name Interactor A`,data_biogrid$`Systematic Name Interactor B`)
colnames(data_biogrid)<-c("X1","X2")
indx <- !duplicated(t(apply(data_biogrid, 1, sort))) # finds non - duplicates in sorted rows
data_biogrid<-data_biogrid[indx, ]
data_biogrid$X1<-as.character(data_biogrid$X1)
data_biogrid$X2<-as.character(data_biogrid$X2)
data_biogrid<-subset(data_biogrid,data_biogrid$X1!=data_biogrid$X2) # drop loop
data_biogrid<-as.data.frame(t(apply(data_biogrid,1,sort)))
colnames(data_biogrid)<-c("X1","X2")
data_biogrid$biogrid<-1

## feature2 (domain-domain interaction)
domin_int<-read.table("~/MyResearch/networks/UniDomInt/ara/4ara_int.txt",stringsAsFactors = TRUE) # have sorted

## feature3 (ractome interaction)
rec<-read.table("~/MyResearch/networks/reactome/3pair.int.txt")
rec<-as.data.frame(t(apply(rec,1,sort)))
indx <- !duplicated(rec) # finds non - duplicates in sorted rows
rec<-rec[indx, ]
colnames(rec)<-c("X1","X2")
rec$reactome<-1

## feature 4 (phylogenetic profile)

phy1<-Data_Preprocess("~/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair.ensembl_bacterial.gene2species")
phy2<-Data_Preprocess("~/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair_ensembl.fungi.gene2species")
phy3<-Data_Preprocess("~/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair.ensembl.gene2species")
phy4<-Data_Preprocess("~/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair_metazoa.gene2species")
phy5<-Data_Preprocess("~/MyResearch/networks/phylogenomes/blast_out/phylogenetic_profile_out/gene2species/tair_phytozome.gene2species")

profile_gene_list1 <-phy1[[1]]
profile_cor1<-Matrix_Transform(phy1[[2]])
gc()
profile_gene_list2 <-phy2[[1]]
profile_cor2<-Matrix_Transform(phy2[[2]])
gc()
profile_gene_list3 <-phy3[[1]]
profile_cor3<-Matrix_Transform(phy3[[2]])
gc()
profile_gene_list4 <-phy4[[1]]
profile_cor4<-Matrix_Transform(phy4[[2]])
gc()
profile_gene_list5<-phy5[[1]]
profile_cor5<-Matrix_Transform(phy5[[2]])
gc()

##### feature5 (TFBS-data)
tf_bind_data<-Data_Preprocess("/home/wuzefeng/MyResearch/motif_dbs/5jaspar2018/1fimo_out/Tair2_p_0.01/gene2motif_0_1.txt")  # tf binding # updataed to p-0.001, and jaspar2018
tf_gene_list<-tf_bind_data[[1]]
tf_cor<-Matrix_Transform(tf_bind_data[[2]]) #slow


#### feature 6 (String database)
String_data<-read.table("/home/wuzefeng/MyResearch/networks/String_10.5/1sort_average_uniq.txt",stringsAsFactors = FALSE,header = TRUE)

##########  expression_data as first input to initation of prediction lines

express_data<-read.table("/home/wuzefeng/MyResearch/RNA_Seq/Plants/Arabidopsis/6cluster_out/cluster_82.exp",stringsAsFactors = FALSE) #expression data
express_cor<-Gene_RNA_seq_pairs(express_data)  #slow (5min)
express_cor.df<-data.frame(row=rownames(express_cor)[row(express_cor)[upper.tri(express_cor)]],col=colnames(express_cor)[col(express_cor)[upper.tri(express_cor)]], pcc=express_cor[upper.tri(express_cor)]) #(slow~20min)

for (m in seq(1,nrow(express_cor.df),1000000)){
  q<-m+999999
  if (q>nrow(express_cor.df)){q<-nrow(express_cor.df)}
  message(m,":",q)
  gc()
  sub_express_unsort <-express_cor.df[m:q,]
  sub_express<-data.frame(t(apply(sub_express_unsort[-3],1,sort)))
  sub_express$pcc<-sub_express_unsort$pcc
  colnames(sub_express)<-c("X1","X2","pcc")
  message("Sub expression dataframe done!")
  gc()
  ## reactome
  sub_express<-merge(sub_express,rec,by = c("X1","X2"),all.x = TRUE)
  message("reactome integration done!")
  ### domin interaction data
  colnames(domin_int)<-c("X1","X2","domain")
  sub_domin_int<-subset(domin_int,domin_int$X1%in%sub_express$X1&domin_int$X2%in%sub_express$X2)
  sub_domin_int$domain<-cut(sub_domin_int$domain,breaks = c(0,0.2,0.4,0.6,0.8,1),labels = letters[1:5])
  sub_express<-merge(sub_express,sub_domin_int,by = c("X1","X2"),all.x = TRUE)
  message("domin integration done!")
  ### phylogentic profile
  sub_express$profile_phy1<-apply(sub_express,1,function(x) profile_cor1[match(x[1],names(profile_gene_list1)),match(x[2],names(profile_gene_list1))]) # esist na
  sub_express$profile_phy1<-cut(sub_express$profile_phy1,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
  message("phy1 integration done!")
  sub_express$profile_phy2<-apply(sub_express,1,function(x) profile_cor2[match(x[1],names(profile_gene_list2)),match(x[2],names(profile_gene_list2))]) # esist na
  sub_express$profile_phy2<-cut(sub_express$profile_phy2,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
  message("phy2 integration done!")
  sub_express$profile_phy3<-apply(sub_express,1,function(x) profile_cor3[match(x[1],names(profile_gene_list3)),match(x[2],names(profile_gene_list3))]) # esist na
  sub_express$profile_phy3<-cut(sub_express$profile_phy3,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
  message("phy3 integration done!")
  sub_express$profile_phy4<-apply(sub_express,1,function(x) profile_cor4[match(x[1],names(profile_gene_list4)),match(x[2],names(profile_gene_list4))]) # esist na
  sub_express$profile_phy4<-cut(sub_express$profile_phy4,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
  message("phy4 integration done!")
  sub_express$profile_phy5<-apply(sub_express,1,function(x) profile_cor5[match(x[1],names(profile_gene_list5)),match(x[2],names(profile_gene_list5))]) # esist na
  sub_express$profile_phy5<-cut(sub_express$profile_phy5,breaks = c(0,1,2,3,4,5,6,Inf),labels = letters[1:7])
  message("phy5 integration done!")
  gc()
  ## ppi
  sub_express<-merge(sub_express,data_biogrid,by = c("X1","X2"),all.x = TRUE)
  message("ppi integration done!")
  prediction_input_file<-data.frame(item=paste(sub_express$X1,sub_express$X2,sep="_"), 
                                    sub_express$pcc,
                                    sub_express$biogrid,
                                    sub_express$profile_phy1, 
                                    sub_express$profile_phy2,
                                    sub_express$profile_phy3,
                                    sub_express$profile_phy4,
                                    sub_express$profile_phy5,
                                    sub_express$reactome,
                                    sub_express$domain)
  
  write.table(prediction_input_file,"prediction_set_continous.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE,append = TRUE)
  gc()
}

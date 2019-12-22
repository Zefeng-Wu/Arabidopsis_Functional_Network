data<-read.table("train_data2.txt",header = TRUE,stringsAsFactors = TRUE)
library(ggplot2)
library(ggsignif)
plot<-list()
for (m in colnames(data)[5:11]){
  plot[[m]]<-ggplot(data,aes_string(x="class",y=m))+geom_boxplot(outlier.color = "NA")+ylim(0,3)+geom_signif(comparisons = list(c("pos","neg")))#,y_position = 2.8
}
library(gridExtra)
do.call(grid.arrange,plot)

#####
p<-ggplot(data,aes(y=profile_phy1,x=class))+geom_boxplot(outlier.color = "NA")+xlab("")+theme(text=element_text(size=18))+ylab("Phylogenetic profile similarity (bacterial)")#+scale_fill_manual(values=c("#E69F00", "#56B4E9"))library(ggsignif)
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE)+scale_x_discrete(labels=c(pos = "GSP",  neg = "GSN"))


p<-ggplot(data,aes(y=profile_phy2,x=class))+geom_boxplot(outlier.color = "NA")+xlab("")+theme(text=element_text(size=18))+ylab("Phylogenetic profile similarity (fungi)")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE)+scale_x_discrete(labels=c(pos = "GSP",  neg = "GSN"))


p<-ggplot(data,aes(y=profile_phy3,x=class))+geom_boxplot(outlier.color = "NA")+xlab("")+theme(text=element_text(size=16))+ylab("Phylogenetic profile similarity (vertebrates)")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE)+scale_x_discrete(labels=c(pos = "GSP",  neg = "GSN"))

p<-ggplot(data,aes(y=profile_phy4,x=class))+geom_boxplot(outlier.color = "NA")+xlab("")+theme(text=element_text(size=18))+ylab("Phylogenetic profile similarity(metazoa)")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE)+scale_x_discrete(labels=c(pos = "GSP",  neg = "GSN"))


p<-ggplot(data,aes(y=profile_phy5,x=class))+geom_boxplot(outlier.color = "NA")+xlab("")+theme(text=element_text(size=18))+ylab("Phylogenetic profile similarity (plants)")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE)+scale_x_discrete(labels=c(pos = "GSP",  neg = "GSN"))


p<-ggplot(data,aes(y=tf,x=class))+geom_boxplot(outlier.color = "NA")+xlab("")+theme(text=element_text(size=18))+ylab("TF binding similarity (PCC)")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE)+scale_x_discrete(labels=c(pos = "GSP",  neg = "GSN"))


#expresison
p<-ggplot(data,aes(y=pcc,x=class))+geom_boxplot(outlier.color = "NA")+xlab("")+theme(text=element_text(size=18))+ylab("Gene expression similarity (PCC)")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE)+scale_x_discrete(labels=c(pos = "GSP",  neg = "GSN"))

##  string interaction score
p<-ggplot(data,aes(y=score,x=class))+geom_boxplot(outlier.color = "NA")+xlab("")+theme(text=element_text(size=18))+ylab("Interaction score")
p+geom_signif(comparisons = list(c("pos", "neg")), map_signif_level=TRUE)+scale_x_discrete(labels=c(pos = "GSP",  neg = "GSN"))


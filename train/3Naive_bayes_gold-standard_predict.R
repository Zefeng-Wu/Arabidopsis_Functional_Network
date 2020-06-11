#### method3: binninig (best)
####  quantile cut into bin to prediction # before or after imputation (discreation and indispensiable is preferred before naive bayes) (best)
train_data <- read.table("/home/wuzefeng/MyResearch/networks/1feature_plot/2018.4.16/train_data2.txt",sep="\t",stringsAsFactors = FALSE,header  = TRUE)
train_data <- train_data[,4:ncol(train_data)]
train_data$class<-as.factor(train_data$class)
require(OneR)
options(digits = 4)

### record bin standard by quantile
pcc_quantiles <- as.numeric(quantile(train_data$pcc,c(0.2,0.4,0.6,0.8,1),na.rm = TRUE))
tf_quantiles <- as.numeric(quantile(train_data$tf,c(0.2,0.4,0.6,0.8,1),na.rm = TRUE))
phy1_quantiles <-as.numeric(quantile(train_data$profile_phy1,c(0.2,0.4,0.6,0.8,1),na.rm = TRUE))
phy2_quantiles <-as.numeric(quantile(train_data$profile_phy2,c(0.2,0.4,0.6,0.8,1),na.rm = TRUE))
phy3_quantiles <-as.numeric(quantile(train_data$profile_phy3,c(0.2,0.4,0.6,0.8,1),na.rm = TRUE))
phy4_quantiles <-as.numeric(quantile(train_data$profile_phy4,c(0.2,0.4,0.6,0.8,1),na.rm = TRUE))
phy5_quantiles <-as.numeric(quantile(train_data$profile_phy5,c(0.2,0.4,0.6,0.8,1),na.rm = TRUE))
domain_quantiles<-as.numeric(quantile(train_data$domain,c(0.2,0.4,0.6,0.8,1),na.rm = TRUE))
string_quantiles <- as.numeric(quantile(train_data$score,c(0.2,0.4,0.6,0.8,1),na.rm = TRUE))


## binning
train_data$pcc<-bin(train_data$pcc,nbins = 4,labels = c(1,2,3,4),method = "content",na.omit = FALSE) # binning
train_data$tf<-bin(train_data$tf,nbins = 4,labels = c(1,2,3,4),method = "content",na.omit = FALSE)
train_data$profile_phy1<-bin(train_data$profile_phy1,nbins = 4,labels = c(1,2,3,4),method = "content",na.omit = FALSE)
train_data$profile_phy2<-bin(train_data$profile_phy2,nbins = 4,labels = c(1,2,3,4),method = "content",na.omit = FALSE)
train_data$profile_phy3<-bin(train_data$profile_phy3,nbins = 4,labels = c(1,2,3,4),method = "content",na.omit = FALSE)
train_data$profile_phy4<-bin(train_data$profile_phy4,nbins = 4,labels = c(1,2,3,4),method = "content",na.omit = FALSE)
train_data$profile_phy5<-bin(train_data$profile_phy5,nbins = 4,labels = c(1,2,3,4),method = "content",na.omit = FALSE)
train_data$domain<-bin(train_data$domain,nbins = 4,labels = c(1,2,3,4),method = "content",na.omit = FALSE)
train_data$score<-bin(train_data$score,nbins = 4,labels = c(1,2,3,4),method = "content",na.omit = FALSE)



###### binging function by quantile 
binnings <-function(df){
  ## pcc qunatiles
  df$pcc[df$pcc<=pcc_quantiles[2]]<-1
  df$pcc[df$pcc>pcc_quantiles[2]&df$pcc<=pcc_quantiles[3]]<-2
  df$pcc[df$pcc>pcc_quantiles[3]&df$pcc<=pcc_quantiles[4]]<-3
  df$pcc[df$pcc>pcc_quantiles[4]&df$pcc<1]<-4
  
  ### tf quantiles
  df$tf[df$tf<=tf_quantiles[2]]<-1
  df$tf[df$tf>tf_quantiles[2]&df$tf<=tf_quantiles[3]]<-2
  df$tf[df$tf>tf_quantiles[3]&df$tf<=tf_quantiles[4]]<-3
  df$tf[df$tf>tf_quantiles[4]&df$tf<1]<-4
  
  ### phy1 quantiles
  df$profile_phy1[df$profile_phy1<=phy1_quantiles[2]]<-1
  df$profile_phy1[df$profile_phy1>phy1_quantiles[2]&df$profile_phy1<=phy1_quantiles[3]]<-2
  df$profile_phy1[df$profile_phy1>phy1_quantiles[3]&df$profile_phy1<=phy1_quantiles[4]]<-3
  df$profile_phy1[df$profile_phy1>phy1_quantiles[4]&df$profile_phy1<1]<-4
  
  ##phy2 quantiles
  df$profile_phy2[df$profile_phy2<=phy2_quantiles[2]]<-1
  df$profile_phy2[df$profile_phy2>phy2_quantiles[2]&df$profile_phy2<=phy2_quantiles[3]]<-2
  df$profile_phy2[df$profile_phy2>phy2_quantiles[3]&df$profile_phy2<=phy2_quantiles[4]]<-3
  df$profile_phy2[df$profile_phy2>phy2_quantiles[4]&df$profile_phy2<1]<-4
  
  ##phy3 quantiles
  df$profile_phy3[df$profile_phy3<=phy3_quantiles[2]]<-1
  df$profile_phy3[df$profile_phy3>phy3_quantiles[2]&df$profile_phy3<=phy3_quantiles[3]]<-2
  df$profile_phy3[df$profile_phy3>phy3_quantiles[3]&df$profile_phy3<=phy3_quantiles[4]]<-3
  df$profile_phy3[df$profile_phy3>phy3_quantiles[4]&df$profile_phy3<1]<-4
  ##phy4 quantiles
  df$profile_phy4[df$profile_phy4<=phy4_quantiles[2]]<-1
  df$profile_phy4[df$profile_phy4>phy4_quantiles[2]&df$profile_phy4<=phy4_quantiles[3]]<-2
  df$profile_phy4[df$profile_phy4>phy4_quantiles[3]&df$profile_phy4<=phy4_quantiles[4]]<-3
  df$profile_phy4[df$profile_phy4>phy4_quantiles[4]&df$profile_phy4<1]<-4
  ##phy5 quantiles
  df$profile_phy5[df$profile_phy5<=phy5_quantiles[2]]<-1
  df$profile_phy5[df$profile_phy5>phy5_quantiles[2]&df$profile_phy5<=phy5_quantiles[3]]<-2
  df$profile_phy5[df$profile_phy5>phy5_quantiles[3]&df$profile_phy5<=phy5_quantiles[4]]<-3
  df$profile_phy5[df$profile_phy5>phy5_quantiles[4]&df$profile_phy5<1]<-4
  
  ##domain quantiles
  df$domain[df$domain<=domain_quantiles[2]]<-1
  df$domain[df$domain>domain_quantiles[2]&df$domain<=domain_quantiles[3]]<-2
  df$domain[df$domain>domain_quantiles[3]&df$domain<=domain_quantiles[4]]<-3
  df$domain[df$domain>domain_quantiles[4]&df$domain<1]<-4
  
  ##String quantiles
  df$score[df$score<=string_quantiles[2]]<-1
  df$score[df$score>string_quantiles[2]&df$score<=string_quantiles[3]]<-2
  df$score[df$score>string_quantiles[3]&df$score<=string_quantiles[4]]<-3
  df$score[df$score>string_quantiles[4]]<-4
  
  ## as levels
  return(df)
}

#####
##modeling and prediction 
require(e1071)
Naive_Bayes_Model <- naiveBayes(class ~., data=train_data) # train model
NB_Predictions <- predict(Naive_Bayes_Model,train_data,type="raw")  # prediction with post probabality
table(ifelse(NB_Predictions[,1]>=NB_Predictions[,2],"neg","pos"),train_data$class) # confuse matrix

require(ggplot2)
require(precrec)
eval1 <- evalmod(scores = NB_Predictions[,2],labels = train_data$class)
autoplot(eval1)
aucs <- auc(eval1) # auc 0.77

####  import integrated data line by line to prediction
##method1

con <- file("/home/wuzefeng/MyResearch/networks/2network_prediction/2result/prediction_set_continous.txt", "rt") # input file
outfile<- file('yourOutputFile', 'w')  # output line

while (length(input<- readLines(con, n=1000))>0){ 
  for (i in 2:length(input)){ 
       temp_line_df<-read.table(textConnection(input[[i]]))
    
       colnames(temp_line_df)<-c("X1","X2",colnames(train_data)[-1])
       temp_line_df<-binnings(temp_line_df) # bininig
       ##prediction
       prediction <- format(predict(Naive_Bayes_Model,temp_line_df,type="raw") , digits = 3)
       if (as.numeric(prediction[2])-as.numeric(prediction[1])){
          writeLines(paste(c(as.character(temp_line_df$X1),as.character(temp_line_df$X2),prediction),collapse = "\t"), con=outfile)
       }
      }
}
close(outfile)


#method2 (quick)
con <- file("/home/wuzefeng/MyResearch/networks/2network_prediction/2result/prediction_set_continous2.txt", "rt") # input file
input<-readLines(con,1)

while (length(input<- readLines(con, n=10000))>0){ 
    temp_line_df<-read.table(textConnection(input),stringsAsFactors = FALSE,header = FALSE)
    colnames(temp_line_df)<-c("X1","X2",colnames(train_data)[-1])
    temp_line_df<-binnings(temp_line_df) # bininig
    
    temp_line_df$pcc<-factor(temp_line_df$pcc,levels = c(1,2,3,4,NA))
    temp_line_df$tf<-factor(temp_line_df$tf,levels = c(1,2,3,4,NA))
    temp_line_df$profile_phy1<-factor(temp_line_df$profile_phy1,levels = c(1,2,3,4,NA))
    temp_line_df$profile_phy2<-factor(temp_line_df$profile_phy2,levels = c(1,2,3,4,NA))
    temp_line_df$profile_phy3<-factor(temp_line_df$profile_phy3,levels = c(1,2,3,4,NA))
    temp_line_df$profile_phy4<-factor(temp_line_df$profile_phy4,levels = c(1,2,3,4,NA))
    temp_line_df$profile_phy5<-factor(temp_line_df$profile_phy5,levels = c(1,2,3,4,NA))
    temp_line_df$domain<-factor(temp_line_df$domain,levels = c(1,2,3,4,NA))
    temp_line_df$score<-factor(temp_line_df$score,levels = c(1,2,3,4,NA))
    
    ##prediction
    prediction <-predict(Naive_Bayes_Model,temp_line_df,type="raw")
    predict_df<-cbind(data.frame(X1=temp_line_df$X1,X2=temp_line_df$X2),as.data.frame(prediction))
    write.table(format(predict_df[predict_df$pos/predict_df$neg>=9.141,],digit=4),file = "3result/tissue_naive_prediction2.txt",append = TRUE,row.names = FALSE,col.names=!file.exists("3result/tissue_naive_prediction.txt"),quote = FALSE)
}









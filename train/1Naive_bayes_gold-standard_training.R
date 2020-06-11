### multiple plot
score_list <-list()

##################################################################################
## method1 no impution with continous

train_data <- read.table("/home/wuzefeng/MyResearch/networks/2.5feature_plot/2018.4.16/train_data2.txt",sep="\t",stringsAsFactors = FALSE,header  = TRUE)
train_data <- train_data[,4:ncol(train_data)]
train_data$class<-as.factor(train_data$class)

require(e1071)
Naive_Bayes_Model <- naiveBayes(class ~., data=train_data) # train model
NB_Predictions <- predict(Naive_Bayes_Model,train_data,type="raw")  # prediction with post probabality
table(ifelse(NB_Predictions[,1]>=NB_Predictions[,2],"neg","pos"),train_data$class) # confuse matrix

  ### too slow
  #require(pROC)
  #modelroc <- roc(train_data[,target_variable],NB_Predictions)
  #plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
       #grid.col=c("green", "red"), max.auc.polygon=TRUE,
       #auc.polygon.col="skyblue", print.thres=TRUE)
require(ggplot2)
require(precrec)
eval1 <- evalmod(scores = NB_Predictions[,2],labels = train_data$class)
autoplot(eval1)
aucs <- auc(eval1) # auc 0.70

score_list[[1]]<-NB_Predictions[,2]


### continous with imputation
## continous wirth imputation
train_data <- read.table("/home/wuzefeng/MyResearch/networks/2.5feature_plot/2018.4.16/train_data2.txt",sep="\t",stringsAsFactors = FALSE,header  = TRUE)
train_data <- train_data[,4:ncol(train_data)]
train_data2<-train_data[,2:ncol(train_data)]
train_data2 <-as.data.frame(impute(train_data2)) #impute
train_data2$class<-train_data$class
train_data<-train_data2
train_data$class<-as.factor(train_data$class)    # as.factor for character variable

require(e1071)
Naive_Bayes_Model <- naiveBayes(class ~., data=train_data) # train model
NB_Predictions <- predict(Naive_Bayes_Model,train_data,type="raw")  # prediction with post probabality
table(ifelse(NB_Predictions[,1]>NB_Predictions[,2],"neg","pos"),train_data$class) # confuse matrix

require(ggplot2)
require(precrec)
eval1 <- evalmod(scores = NB_Predictions[,2],labels = train_data$class)
autoplot(eval1)
aucs <- auc(eval1) # auc 0.69

score_list[[2]]<-NB_Predictions[,2]

#### method2: binninig (best)
####  quantile cut into bin to prediction # before or after imputation (discreation and indispensiable is preferred before naive bayes) (best)
train_data <- read.table("/home/wuzefeng/MyResearch/networks/2.5feature_plot/2018.4.16/train_data2.txt",sep="\t",stringsAsFactors = FALSE,header  = TRUE)
train_data <- train_data[,4:ncol(train_data)]
train_data$class<-as.factor(train_data$class)
require(OneR)
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

#train_data$biogrid<-ifelse(is.na(train_data$biogrid),0,1)
#train_data$domain<-ifelse(is.na(train_data$domain),0,1)
#score_list[[4]]<-NB_Predictions[,2] #run at end! 

##predict
require(e1071)
Naive_Bayes_Model <- naiveBayes(class ~., data=train_data) # train model
NB_Predictions <- predict(Naive_Bayes_Model,train_data,type="raw")  # prediction with post probabality
table(ifelse(NB_Predictions[,1]>=NB_Predictions[,2],"neg","pos"),train_data$class) # confuse matrix


require(ggplot2)
require(precrec)
eval1 <- evalmod(scores = NB_Predictions[,2],labels = train_data$class)
autoplot(eval1)
aucs <- auc(eval1) # auc 0.77
score_list[[3]]<-NB_Predictions[,2]

##################### mutiple methods for same data
 
label<-train_data$class
modnames <-c("continous","continous imputation","binning continous","binning continous and ppi NA to 0")
eval_all <- mmdata(list(score_list),labels = train_data$class,modnames = modnames)
autoplot(evalmod(eval_all))



### methods3: binning and  cross_validation:
### method3: bining and cross validations
require(precrec)
train_data <- read.table("/home/wuzefeng/MyResearch/networks/2.5feature_plot/2018.4.16/train_data2.txt",sep="\t",stringsAsFactors = FALSE,header  = TRUE)
train_data <- train_data[,4:ncol(train_data)]
train_data$class<-as.factor(train_data$class)
require(OneR)
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


#cross validation
score_list<-list()  # record each cor
label_list<-list()
modenames<- paste("Cross validation",seq(1:10),sep="")
dsids <-seq(0,9)

nrFolds <- 10
folds <- rep_len(1:nrFolds, nrow(train_data))
folds <- sample(folds, nrow(train_data))  # https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation

par(mfrow=c(4,3),mar=c(4,5,2,1))
for(k in 1:nrFolds) {
  message(k)
  require(e1071)
  fold <- which(folds == k)
  data.train <- train_data[-fold,]
  data.test <- train_data[fold,]
  message(dim(data.train),"**",dim(data.test)) # 90% used to train, 10% used to validate
  model<- naiveBayes(class ~., data=data.train) 
  predictions<-predict(model,data.test,type="raw")
  
  ##
  require(ggplot2)
  require(precrec)
  #eval1 <- evalmod(scores = predictions[,2],labels = data.test$class)
  #autoplot(eval1)
  #print(auc(eval1))
  
  ## add socre list
  score_list[[k]]<-predictions[,2]
  label_list[[k]]<-data.test$class
}

mdata <- mmdata(score_list, labels = label_list,modnames = modenames,dsids = dsids)
mmcurves <- evalmod(mdata)
autoplot(mmcurves)

#### 






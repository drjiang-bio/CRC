setwd('E:\\Rwork\\machine_learning\\17-12-26')   #设置工作路径g
library(Biobase)
library(genefilter)
###---外部测试集预处理---
ceshi <- read.csv('./original/ceshi.csv',header = T,row.names=1,stringsAsFactors = F)
str(ceshi[,1:4])   #检查数据
#---
ceshi_m <- as.matrix(ceshi)
colnames(ceshi_m) <- c(1:length(colnames(ceshi_m)))
summary(as.vector(ceshi_m))   #描述测试集总体分布
Z#---
hist(as.vector(ceshi_m),breaks = 100,prob=T,xlab = 'Expression Levels',
     main = 'Histogram of Overall Expression Levels')  #直方图
abline(v = c(median(as.vector(ceshi_m),na.rm = T),shorth(as.vector(ceshi_m),na.rm = T),
             quantile(as.vector(ceshi_m),c(0.25,0.75),na.rm = T)),
       lty=2,col=c(2,3,4,4))  #加参考线
legend('topright',c('Median','Shorth','1stQ','3rdQ'),
       lty = 2,col = c(2,3,4,4))   #加图例
#---
save(ceshi,file = './cache/ceshi.Rdata')
write.csv(ceshi,file = './cache/ceshi.csv')  #保存中间数据
###---训练集预处理---<24,25,26>
eso_cer <- read.csv('./original/eso_cer(logFC6).csv',header = T,stringsAsFactors = F)
eso <- read.csv('./original/GSE45670_series_matrix.csv',header = T,row.names=1,stringsAsFactors = F)
cer <- read.csv('./original/GSE63514_series_matrix.csv',row.names=1,header = T,stringsAsFactors = F)
cer[1,grep('^N',cer[1,])] <- 30
cer[1,grep('^Ca',cer[1,])] <- 31
cer[,grep('^C',cer[1,])] <- NULL

eso[1,grep('^n',eso[1,])] <- 10
eso[1,grep('^e',eso[1,])] <- 11

eso_m <- apply(eso[-1,],1,function(x) as.numeric(x))
cer_m <- apply(cer[-1,],1,function(x) as.numeric(x))
summary(as.vector(eso_m))
summary(as.vector(cer_m))
t.test(colMeans(eso_m),colMeans(cer_m))

hist(as.vector(eso_m),breaks = 100,prob=T,xlab = 'Expression Levels',
     main = 'Histogram of Overall Expression Levels')  #直方图
abline(v = c(median(as.vector(eso_m),na.rm = T),shorth(as.vector(eso_m),na.rm = T),
             quantile(as.vector(eso_m),c(0.25,0.75),na.rm = T)),
       lty=2,col=c(2,3,4,4))  #加参考线
legend('topright',c('Median','Shorth','1stQ','3rdQ'),
       lty = 2,col = c(2,3,4,4))   #加图例

hist(as.vector(cer_m),breaks = 100,prob=T,xlab = 'Expression Levels',
     main = 'Histogram of Overall Expression Levels')  #直方图
abline(v = c(median(as.vector(cer_m),na.rm = T),shorth(as.vector(cer_m),na.rm = T),
             quantile(as.vector(cer_m),c(0.25,0.75),na.rm = T)),
       lty=2,col=c(2,3,4,4))  #加参考线
legend('topright',c('Median','Shorth','1stQ','3rdQ'),
       lty = 2,col = c(2,3,4,4))   #加图例
##---特征选择---
###1.纯两者差异基因
cer_f <- cer[c(1,which(row.names(cer) %in% eso_cer$ID)),]
eso_f <- eso[c(1,which(row.names(eso) %in% eso_cer$ID)),]
cer_t <- apply(cer_f,1,function(x) as.numeric(x))
eso_t <- apply(eso_f,1,function(x) as.numeric(x))
cer_t <- as.data.frame(cer_t)
eso_t <- as.data.frame(eso_t)
cer_t$Class <- factor(cer_t$Class, levels=c(30,31), 
                   labels=c("Normal_cer", "cer"))
eso_t$Class <- factor(eso_t$Class, levels=c(10,11), 
                      labels=c("Normal_eso", "eso"))
train <- rbind(cer_t[25:52,],eso_t[11:38,])
train$Class <- factor(train$Class)
train <- na.omit(train)
summary(as.vector(as.matrix((train[,-1]))))
test <- as.data.frame(t(ceshi))
##总结数据
summary(as.vector(eso_m))
summary(as.vector(cer_m))
summary(as.vector(ceshi_m))
summary(as.vector(as.matrix((train[,-1]))))
layout(matrix(c(1,2,3,4),2,2,byrow=T))
#4图
#内部测试
library(e1071)
set.seed(1234)
n <- sample(nrow(train), 0.7*nrow(train))
train_n <- train[n,]
test_n <- train[-n,]
table(train_n$Class)
table(test_n$Class)

tuned <- tune.svm(Class~., data=train,
                  gamma=10^(-6:1),cost=10^(-10:10))
tuned
fit.svm1 <- svm(Class~., data=na.omit(train_n))#, gamma=.01, cost=1
fit.svm1
svm.pred1 <- predict(fit.svm1, na.omit(test_n))
svm.pred1
svm.perf1 <- table(na.omit(test_n)$Class,
                  svm.pred1, dnn=c("Actual", "Predicted"))
svm.perf1
#外部测试
set.seed(1234)
tuned <- tune.svm(Class~., data=train,
                  gamma=10^(-6:1),cost=10^(-10:10))
tuned
fit.svm <- svm(Class~., data=na.omit(train))#, gamma=.01, cost=1
fit.svm
svm.pred <- predict(fit.svm, na.omit(test))
svm.pred
svm.perf <- table(na.omit(test)$Class,
                  svm.pred, dnn=c("Actual", "Predicted"))
svm.perf
###2.纯各自差异基因交集
###3.过滤
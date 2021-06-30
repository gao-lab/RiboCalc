library('caret')
library('glmnet')
load('../RiboCalc/RiboCalc.RData')

# read data
yeast <- read.table(file='raw_data/coding_features.yeast.data',header=T)
yeast_Li <- read.csv(file='raw_data/nar-00812-a-2017-File019.filt.csv',header=T)

# format yeast data for RiboCalc human model
yeast.tmp <- yeast
yeast.tmp$X5UTR_length <- log2(yeast.tmp$UTR5_length+1)
yeast.tmp$X3UTR_length <- log2(yeast.tmp$UTR3_length+1)
yeast.tmp$X3UTR_GC <- yeast.tmp$UTR3_GC
yeast.tmp$X5UTR_GC <- yeast.tmp$UTR5_GC
yeast.tmp$CDS_length <- log2(yeast.tmp$CDS_length+1)
yeast.tmp$Length <- log2(yeast.tmp$Length+1)
yeast.tmp$RNA_TPM <- log2(10^yeast.tmp$RNA)
yeast.tmp$Ribo_TPM <- yeast.tmp$RNA_TPM + log2(10^yeast.tmp$TR)

# intersect features
human_features <- colnames(training.filt)
yeast_features <- colnames(yeast.tmp)
features.tmp <- intersect(human_features,yeast_features)

# scaled as human training data
yeast.tmp.scaled <- yeast.tmp[,features.tmp]
for (i in features.tmp){
  yeast.tmp.scaled[,i] <- normalization2(x=yeast.tmp[,i],y=training.raw[,i])}

# fill abscent features in yeast with zeros
data <- training.filt[1:nrow(yeast.tmp.scaled),]
data[,] <- 0
data[,features.tmp] <- yeast.tmp.scaled[,features.tmp]
x.test.yeast <- model.matrix(Ribo_TPM ~., data)[,-1]
pred.human <- predict(lasso.model, newx = x.test.yeast)

# RiboCalc human model testing for Ribo-TPM
set.seed(1000)
inTrain <- sample(length(yeast$TR),1000)
cor(pred.human[-inTrain], yeast.tmp$Ribo_TPM[-inTrain]) #Table 3: RiboCalc human Ribo-TPM

# RiboCalc human model testing for TR
scaling_line <- lm(yeast.tmp$Ribo_TPM[inTrain]~pred.human[inTrain])
pred.human.tmp <- pred.human[-inTrain]*18.69025 - 13.35659 #scaling_line coefficient and intercept
pred.human.TR <- log10(2^(pred.human.tmp - yeast.tmp$RNA_TPM[-inTrain]))
cor(pred.human.TR, yeast.tmp$TR[-inTrain]) #Table 3: RiboCalc human TR


# Build RiboCalc yeast model
# format data
yeast$UTR5_length <- log10(yeast$UTR5_length)
yeast$UTR3_length <- log10(yeast$UTR3_length)
yeast$CDS_length <- log10(yeast$CDS_length)
yeast$Length <- log10(yeast$Length)
yeast$prot <- yeast$RNA + yeast$TR

# scale data
training.yeast <- as.data.frame(sapply(yeast[inTrain,-1],normalization1))
testing.yeast <- yeast[-inTrain,-1]
for (i in 1:ncol(yeast[-inTrain,-1])){
  testing.yeast[,i] <- normalization2(x=yeast[-inTrain,-1][i],y=yeast[inTrain,-1][i])}

# RiboCalc yeast: lasso regression model
set.seed(123)
# CDS length and MTDR were removed because of high correlation
# start codon PWMs were removed because they are human motifs
cv.yeast <- cv.glmnet(x = as.matrix(training.yeast[,-c(66:70,79:81)]), y = training.yeast$prot, alpha = 1,nfolds = 5)
yeast.lasso.model<-glmnet(x = as.matrix(training.yeast[,-c(66:70,79:81)]), y = training.yeast$prot, alpha = 1, lambda = cv.yeast$lambda.min)
#RiboCalc yeast model coefficients
coef(yeast.lasso.model) #Figure 2d/2c y axis, Table S6 

pred.yeast <- predict(yeast.lasso.model,newx = as.matrix(testing.yeast[,-c(66:70,79:81)]))
pred.yeast.tmp<- recover_scaling(pred.yeast, yeast$prot[inTrain])
pred.yeast.TR <- pred.yeast.tmp[,1] - yeast[-inTrain,"RNA"]

# Li's model
yeast_Li_model <- lm(data = yeast_Li[inTrain,-c(1:4,6,7,17)], log10.TR~.)
pred.Li.TR <- predict(yeast_Li_model,newdata = yeast_Li[-inTrain,])
pred.Li <- pred.Li.TR + yeast_Li[-inTrain,"log10.RNA"]

# Testing: correlation between predicted and observed values
cor(pred.yeast[,1],testing.yeast$prot) #Table 3: RiboCalc yeast Ribo-TPM
cor(pred.yeast.TR, yeast[-inTrain,"TR"]) #Table 3: RiboCalc yeast TR
cor(pred.Li, yeast_Li[-inTrain,"log10.TR"]+yeast_Li[-inTrain,"log10.RNA"]) #Table 3: Li's model RiboTPM
cor(pred.Li.TR, yeast_Li[-inTrain,"log10.TR"]) #Table 3: Li's model TR

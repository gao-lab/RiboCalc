library('caret')
library('glmnet')

# read data path
#args=commandArgs(T)
#sid <- args[1]  #SRR1803151, SRX870805, SRR970565, SRR627625, SRR3208870
sid <- "SRR1803151"
file_path <- paste("raw_data/", sid, ".data",sep = "")

# read data and format
data <- read.table(file=file_path, header=T)
format_data <- function(dt,type){
	data1 <- dt
	data1$X3UTR_length <- log2(dt$UTR3_length+1)
	data1$X3UTR_GC <- dt$UTR3_GC
	data1$X5UTR_length <- log2(dt$UTR5_length+1)
	data1$X5UTR_GC <- dt$UTR5_GC
	data1$CDS_length <- log2(dt$CDS_length+1)
	if (type==1){
		data1 <- data1[dt$Ribo_TPM>0,-c(72:75)]
	}else{
		data1 <- data1[dt$Ribo_TPM==0,-c(72:75)]}}
data <- format_data(data,1)

# split training and testing set
set.seed(1000)
inTrain <- sample(nrow(data),3000)
training.raw <- data[inTrain,-1] # remove IDs
testing.raw <- data[-inTrain,-1] #remove IDs
normalization1 <- function(x){(x-min(x))/(max(x)-min(x))}
normalization2 <- function(x,y){(x-min(y))/(max(y)-min(y))}
recover_scaling <- function(x,y){x*(max(y)-min(y))+min(y)}

training <- as.data.frame(sapply(training.raw,normalization1))
testing <- testing.raw
for (i in 1:ncol(testing)){
	testing[,i] <- normalization2(x=testing.raw[,i],y=training.raw[,colnames(testing)[i]])}


# build models
CorMatrix<-cor(training)
highlyCorrelated<-findCorrelation(CorMatrix,cutoff = .9)
set.seed(123)
cv <- cv.glmnet(x = as.matrix(training[,-78]), y = training$Ribo_TPM, alpha = 1,nfolds = 5)
lasso.model<-glmnet(x = as.matrix(training[,-78]), y = training$Ribo_TPM, alpha = 1, lambda = cv$lambda.min)
x.test <- model.matrix(Ribo_TPM ~., testing)[,-1]
x.pred <- predict(lasso.model,newx = x.test)
cor(x.pred,testing$Ribo_TPM) #Table 1
coef(lasso.model) #Table S6

#TE
x.pred.re <- recover_scaling(x.pred, training.raw$Ribo_TPM)
cor.test(x.pred.re-testing.raw$RNA_TPM, testing.raw$Ribo_TPM-testing.raw$RNA_TPM)

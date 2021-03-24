library('caret')
library('glmnet')
load('../RiboCalc/RiboCalc.RData')

# format data
data <- read.table(file='./feature_data/OCTOPOS_feature.testing.tab',header=T)
data$Length <- log2(data$Length+1)
data$Ribo_TPM <- log2(data$ProtAbundance + 1)
data.tmp <- format_data(data[,colnames(training.filt)])
for (i in colnames(training.filt)){
	data.tmp[,i] <- normalization2(x=data.tmp[,i], y=training.raw[,i])}

# testing directly using RiboCalc model
data.test <- model.matrix(Ribo_TPM ~., data.tmp[,colnames(training.filt)])[,-1]
data.pred.tmp <- predict(lasso.model, newx = data.test)
data.pred <- recover_scaling(data.pred.tmp, training.raw$Ribo_TPM)
cor.test(data.pred, data$Ribo_TPM)

set.seed(100)
inTest <- sample(length(data$Ribo_TPM), 461)
# choose 461 transcripts as testing set as OCTOPOS method (2136 * 30%) 

training <- data.tmp[-inTest, colnames(training.filt)[1:82]]
# remove translation regutor features because they are the same value in one sample
testing <- data.tmp[inTest, colnames(training.filt)[1:82]]
testing$Ribo_TPM <- 0
ctrl <- trainControl(method = "repeatedcv", number = 4, repeats = 5)
lm.refitted <- train(Ribo_TPM ~ ., data = training, method ="lm", trControl=ctrl)
test.pred.tmp  <- predict(lm.refitted, newdata = testing)
test.pred <- recover_scaling(test.pred.tmp, training.raw$Ribo_TPM)
testing$Ribo_TPM <- data$Ribo_TPM[inTest]

cor.test(test.pred, testing$Ribo_TPM)
cor.test(data.pred[inTest], testing$Ribo_TPM)

plot(testing$Ribo_TPM,data.pred[inTest],xlab="log2(Protein Abundance + 1)",ylab="RiboCalc predicted value",pch=20) #Figure S3 a
plot(testing$Ribo_TPM,test.pred,xlab="log2(Protein Abundance + 1)",ylab="RiboCalc re-fitted model predicted value",pch=20) #Figure S3 b

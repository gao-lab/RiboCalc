library('caret')
library('glmnet')

# read data
training <- read.table(file='raw_data/coding_features.train.data',header=T) #training set
testing <- read.table(file='raw_data/coding_features.test.data',header=T) #testing set
noncoding <- read.table(file='raw_data/noncoding_features.test.data',header=T) #noncoding testing set
both_co.all <- read.table(file='raw_data/CDCT_coding.test.data',header=T) #coding CDCT testing set
both_nc.all <- read.table(file='raw_data/CDCT_noncoding.test.data',header=T) #noncoding CDCT testing set

# log transformation of feature value
format_data <- function(data){
	data1 <- data
	data1$X3UTR_length <- log2(data$X3UTR_length+1)
	data1$X5UTR_length <- log2(data$X5UTR_length+1)
	data1$CDS_length <- log2(data$CDS_length+1)
	data1}

training.raw <- format_data(training)
testing.raw <- format_data(testing)

# remove information lines: ID/sample/label/isoform_type
removeLines <- c(1,83:85)

# scale data:
## normalization1 for training data
normalization1 <- function(x){(x-min(x))/(max(x)-min(x))}
## normalization2 for testing data (y is the training data)
normalization2 <- function(x,y){(x-min(y))/(max(y)-min(y))}
## recover to the unscales value
recover_scaling <- function(x,y){x*(max(y)-min(y))+min(y)}

training <- as.data.frame(sapply(training.raw[,-removeLines],normalization1))
testing <- testing.raw[,-removeLines]
for (i in 1:ncol(testing)){
	testing[,i] <- normalization2(x=testing[,i],y=training.raw[,colnames(testing)[i]])}


# remove highly correlated features
CorMatrix<-cor(training)
highlyCorrelated<-findCorrelation(CorMatrix,cutoff = .9)
training.filt<-training[,-highlyCorrelated]

# RiboCalc model: lasso regression
set.seed(123)
cv <- cv.glmnet(x = as.matrix(training.filt[,-81]), y = training.filt$Ribo_TPM, alpha = 1,nfolds = 5)
lasso.model<-glmnet(x = as.matrix(training.filt[,-81]), y = training.filt$Ribo_TPM, alpha = 1, lambda = cv$lambda.min)
x.test <- model.matrix(Ribo_TPM ~., testing[-highlyCorrelated])[,-1]
x.pred <- predict(lasso.model,newx = x.test)
cor(x.pred, testing$Ribo_TPM) #Figure 2a correlation result
plot(x.pred, testing$Ribo_TPM, pch=20, xlab="RiboCalc Predicted", ylab="Observed Ribo-TPM") #Figure 2a
coef(lasso.model) #Table 2, Figure 2d/2c x axis, Table S3

# TE Prediction
x.pred.re <- recover_scaling(x.pred, training.raw$Ribo_TPM)
cor.test(x.pred.re-testing.raw$RNA_TPM, testing.raw$Ribo_TPM-testing.raw$RNA_TPM)
save("format_data","normalization1","normalization2","recover_scaling",
     "training.raw","training.filt","lasso.model",file='RiboCalc.RData') #save RiboCalc model to RData file


# feature importance
all_features <- rownames(coef(lasso.model))
all_features[1] <- "Ribo_TPM"
init.all <- all_features[c(67:70,79)]
elongation.all <- all_features[c(2:65,76:78)]
express.all <- all_features[c(66,71,81)]
background.all <- all_features[83:206]
structure.all <- all_features[c(72:75,80,82)]

data.tmp <- testing[,all_features]
data.tmp[,] <- 0
# context model (trans-)
data.tmp[,express.all] <- testing[,express.all]
data.tmp[,background.all] <- testing[,background.all]
x.test <- model.matrix(Ribo_TPM ~., data.tmp)[,-1]
cor.test(predict(lasso.model,newx = x.test), testing$Ribo_TPM) #Figure 2b left bar
# intrinsic model (cis-)
data.tmp[,] <- 0
data.tmp[,init.all] <- testing[,init.all]
data.tmp[,elongation.all] <- testing[,elongation.all]
data.tmp[,structure.all] <- testing[,structure.all]
x.test <- model.matrix(Ribo_TPM ~., data.tmp)[,-1]
cor.test(predict(lasso.model,newx = x.test), testing$Ribo_TPM) #Figure 2b right bar

tmp <- rbind(cbind("Translation initiation",coef(lasso.model)[init.all,1]),
             cbind("Translation elongation",coef(lasso.model)[elongation.all,1]),
             cbind("Expression abundance",coef(lasso.model)[express.all,1]),
             cbind("Translation regulators",coef(lasso.model)[background.all,1]),
             cbind("Transcript structure",coef(lasso.model)[structure.all,1]))
features.plot <- as.data.frame(tmp)
colnames(features.plot) <- c("Class","Coefficient")
features.plot$Coefficient <- as.numeric(as.character(features.plot$Coefficient))
features.plot$Class <- factor(features.plot$Class,levels = c("Expression abundance",
				"Translation regulators","Translation initiation",
				"Translation elongation","Transcript structure"))
ggplot(data = features.plot, aes(x=Class,y=log10(abs(Coefficient))))+geom_boxplot() #Figure 2c


# preprocess: noncoding and CDCT data format and scaling
noncoding.raw <- format_data(noncoding)
both_co <- format_data(both_co.all[both_co.all$Ribo_TPM>0,colnames(training.filt)])
both_nc <- format_data(both_nc.all[both_nc.all$Ribo_TPM<=0,colnames(training.filt)])

noncoding <- noncoding.raw[,-removeLines]
for (i in 1:ncol(noncoding)){
	noncoding[,i] <- normalization2(x=noncoding[,i],y=training.raw[,colnames(noncoding)[i]])}
for (i in colnames(training.filt)){
	both_co[,i] <- normalization2(x=both_co[,i],y=training.raw[,i])}
for (i in colnames(training.filt)){
	both_nc[,i] <- normalization2(x=both_nc[,i],y=training.raw[,i])}

noncoding.tmp <- noncoding[,colnames(training.filt)]
x.test.nc <- model.matrix(Ribo_TPM ~., noncoding.tmp)[,-1]
both.test.nc <- model.matrix(Ribo_TPM ~., both_nc[,colnames(training.filt)])[,-1]
both.test.co <- model.matrix(Ribo_TPM ~., both_co[,colnames(training.filt)])[,-1]

# CDCT prediction
nc.fitted <- predict(lasso.model, newx = x.test.nc)
test.fitted <- x.pred
both_co.fitted <- predict(lasso.model, newx = both.test.co)
both_nc.fitted <- predict(lasso.model, newx = both.test.nc)
predicted <- c(nc.fitted, test.fitted, both_co.fitted, both_nc.fitted)
type <- c(rep("noncoding",length(nc.fitted)), rep("coding",length(test.fitted)),
		rep("coding CDCTs",length(both_co.fitted)), rep("noncoding CDCTs",length(both_nc.fitted)))
plotr <- as.data.frame(cbind(predicted,type))
plotr$predicted<-as.numeric(as.character(plotr$predicted))
ggplot(plotr, aes(y=predicted,x=type))+geom_boxplot()+ylab("RiboCalc predicted RiboTPM") #Figure 3a

# noncoding prediction
data <- rbind(training.filt[,all_features],testing[,all_features],noncoding[,all_features])
data_type <- c(rep("coding",nrow(training.filt)+nrow(testing)),rep("noncoding",nrow(noncoding)))
annotation <- as.factor(c(as.character(training.raw$isoform_type),
				as.character(testing.raw$isoform_type),
				as.character(noncoding.raw$isoform_type)))
coding.anno <- levels(annotation)[c(4,6,17,26:29)] #the biotypes of coding RNAs
#noncoding.anno <- levels(annotation)[c(1:3,5,7:16,18:25,30:36)]
# function for defining which RNA biotype is coding and which is noncoding
find_anno <- function(x,y=coding.anno){
  if (is.element(x,y)){
    x <- "mRNAs"
  }else{
    x <- "ncRNAs"}
  x
}
data_type <- paste(data_type,sapply(annotation, find_anno))
x.test <- model.matrix(Ribo_TPM ~., data)[,-1]
x.pred <- predict(lasso.model, newx = x.test)
ggplot(data, aes(y=x.pred,x=data_type))+geom_boxplot()+ylab("RiboCalc predicted RiboTPM") #Figure S9

# read ribo-lnc data
library('caret')
library('glmnet')

load('../RiboCalc/RiboCalc.RData')

BMC <- read.table(file='raw_data/BMC.test',header=T)
elife <- read.table(file='raw_data/elife.test',header=T)

BMC.tmp <- format_data(BMC[,colnames(training.filt)])
elife.tmp <- format_data(elife[,colnames(training.filt)])

for (i in colnames(training.filt)){
  BMC.tmp[,i] <- normalization2(x=BMC.tmp[,i],y=training.raw[,i])}
for (i in colnames(training.filt)){
  elife.tmp[,i] <- normalization2(x=elife.tmp[,i],y=training.raw[,i])}

BMC.test <- model.matrix(Ribo_TPM ~., BMC.tmp[,colnames(training.filt)])[,-1]
elife.test <- model.matrix(Ribo_TPM ~., elife.tmp[,colnames(training.filt)])[,-1]

BMC.pred <- predict(lasso.model,newx = BMC.test)
elife.pred <- predict(lasso.model,newx = elife.test)

BMC.plot <- BMC[,c("Ribo_TPM","Label")]
BMC.plot$Ribo_TPM <- BMC.pred[,1]
ggplot(BMC.plot, aes(y=Ribo_TPM, x=Label)) + geom_boxplot() #Figure 3C

elife.plot <- elife[,c("Ribo_TPM","Label")]
elife.plot$Ribo_TPM <- elife.pred[,1]
ggplot(elife.plot, aes(y=Ribo_TPM, x=Label)) + geom_boxplot() #Figure 3B

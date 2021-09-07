library('caret')
library('glmnet')

# Load data
load('../../RiboCalc/RiboCalc.RData')
data<-read.table(file='./raw_data/feature_data_LiJJ.tab',header=T)

# Format data
data$Length <- log2(data$Length+1)
data$X3UTR_length <- log2(data$X3UTR_length+1)
data$X5UTR_length <- log2(data$X5UTR_length+1)
data$CDS_length <- log2(data$CDS_length+1)
data$RNA_TPM <- log2(data$RNA_TPM)
data$Ribo_TPM <- log2(data$Ribo_TPM)

data.tmp <- format_data(data[,colnames(training.filt)])
for (i in colnames(training.filt)){
	data.tmp[,i] <- normalization2(x=data.tmp[,i], y=training.raw[,i])}

# Testing directly using RiboCalc model
data.test <- model.matrix(Ribo_TPM ~., data.tmp[,colnames(training.filt)])[,-1]
data.pred <- predict(lasso.model, newx = data.test)
cor(data.pred, data.tmp$Ribo_TPM)


# Remove feature lines that are not available in Li's data
data.tmp<-data[,-c(1,9:132,208)]

# Select training data
set.seed(100)
inTest <- sample(length(data$Ribo_TPM), 2000)

# Scale feature value to [0,1] interval
training <- as.data.frame(sapply(data.tmp[-inTest,], normalization1))
testing <- data.tmp[inTest,]

for (i in colnames(testing)){
	testing[,i] <- normalization2(x=testing[,i], y=data.tmp[-inTest,][i])}

# Re-trained data with RPKM
ctrl <- trainControl(method = "repeatedcv", number = 4, repeats = 5)
lm.refitted <- train(Ribo_TPM ~ ., data = training, method ="lm", trControl=ctrl)


# Ribo-RPKM testing (Table S6)
pred.tmp <- predict(lm.refitted, newdata=testing)
pred <- recover_scaling(pred.tmp, data.tmp$Ribo_TPM[-inTest])

cor(pred, data.tmp$Ribo_TPM[inTest]) #RiboCalc
cor(log2(10^data$pred_TR[inTest])+data.tmp$RNA_TPM[inTest], data.tmp$Ribo_TPM[inTest]) #Li human

# TR testing (Table S6)
cor(pred-data.tmp$RNA_TPM[inTest], data.tmp$Ribo_TPM[inTest]-data.tmp$RNA_TPM[inTest]) #RiboCalc
cor(log2(10^data$pred_TR[inTest]), data.tmp$Ribo_TPM[inTest]-data.tmp$RNA_TPM[inTest]) #Li human

# New model coefficients
all_features <- rownames(coef(lasso.model))
all_features[1] <- "Ribo_TPM"
init.all <- all_features[c(67:70,79)]
elongation.all <- all_features[c(2:65,76:78)]
express.all <- all_features[c(66,71,81)]
background.all <- all_features[83:206]
structure.all <- all_features[c(72:75,80,82)]
tmp <- rbind(cbind("Translation initiation",coef(lm.refitted$finalModel)[init.all]),
	cbind("Translation elongation",coef(lm.refitted$finalModel)[elongation.all]),
	cbind("Expression abundance",coef(lm.refitted$finalModel)[express.all]),
	cbind("Transcript structure",coef(lm.refitted$finalModel)[structure.all]))
tmp[is.na(tmp)] <- 0
features.plot <- as.data.frame(tmp)
colnames(features.plot) <- c("Class","Coefficient")
features.plot$Coefficient <- as.numeric(as.character(features.plot$Coefficient))
features.plot$Class <- factor(features.plot$Class,levels = c("Expression abundance",
	"Translation initiation","Translation elongation","Transcript structure"))
ggplot(data = features.plot, aes(x=Class,y=log10(abs(Coefficient))))+geom_boxplot() #Compared to Figure 2C

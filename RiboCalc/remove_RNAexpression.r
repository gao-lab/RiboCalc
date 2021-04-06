# Since Ribo-seq TPM is reported to be correlated RNA-seq TPM
# We test our model by removing the feature RNA-TPM

library('caret')
library('glmnet')
load('RiboCalc.RData')

# Load testing data
testing <- read.table(file='raw_data/coding_features.test.data',header=T)
testing.raw <- format_data(testing)
removeLines <- c(1,83:85)
testing <- testing.raw[,-removeLines]
for (i in 1:ncol(testing)){
	testing[,i] <- normalization2(x=testing[,i],y=training.raw[,colnames(testing)[i]])}
CorMatrix<-cor(training.raw[,-removeLines])
highlyCorrelated<-findCorrelation(CorMatrix,cutoff = .9)
testing.filt <- testing[-highlyCorrelated]

# Testing RiboCalc directly without RNA TPM feature
testing.filt$RNA_TPM <- 0
x.test <- model.matrix(Ribo_TPM ~., testing.filt)[,-1]
x.pred <- predict(lasso.model,newx = x.test)
cor(x.pred, testing$Ribo_TPM)


# Rebuild model wthout RNA TPM feature
set.seed(123)
cv2 <- cv.glmnet(x = as.matrix(training.filt[,-c(80,81)]), y = training.filt$Ribo_TPM, alpha = 1,nfolds = 5)
lasso.model.tmp<-glmnet(x = as.matrix(training.filt[,-c(80,81)]), y = training.filt$Ribo_TPM, alpha = 1, lambda = cv2$lambda.min)
x.test <- model.matrix(Ribo_TPM ~., testing.filt[,-80])[,-1]
x.pred <- predict(lasso.model.tmp,newx = x.test)
cor.test(testing.filt$Ribo_TPM, x.pred)


# Plot new model feature coefficients
all_features <- rownames(coef(lasso.model))
all_features[1] <- "Ribo_TPM"

init.all <- all_features[c(67:70,79)]
elongation.all <- all_features[c(2:65,76:78)]
express.all <- all_features[c(66,71)]
background.all <- all_features[83:206]
structure.all <- all_features[c(72:75,80,82)]
tmp <- rbind(cbind("Expression abundance",coef(lasso.model.tmp)[express.all,1]),
	cbind("Translation regulators",coef(lasso.model.tmp)[background.all,1]),
	cbind("Translation initiation",coef(lasso.model.tmp)[init.all,1]),
	cbind("Translation elongation",coef(lasso.model.tmp)[elongation.all,1]),
	cbind("Transcript structure",coef(lasso.model.tmp)[structure.all,1]))
features.plot <- as.data.frame(tmp)
colnames(features.plot) <- c("Class","Coefficient")
features.plot$Coefficient <- as.numeric(as.character(features.plot$Coefficient))
features.plot$Class <- factor(features.plot$Class,levels = c("Expression abundance",
				"Translation regulators","Translation initiation",
				"Translation elongation","Transcript structure"))
ggplot(data = features.plot, aes(x=Class,y=log10(abs(Coefficient))))+geom_boxplot()

### Load the packages and create necessary functions

library(randomForest)
library(Seurat)
library(dplyr)
library(reshape2)


plotConfusionMatrix = function(X, x.order=default.x.order, y.order=default.y.order, row.scale=TRUE, col.scale=FALSE, cols.use=gray.colors(10), max.size=5, ylab.use="Known", xlab.use="Predicted"){
  
  if (!col.scale & row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
  if (col.scale & !row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
  if(col.scale & row.scale){
    print("Only one of row.scale or col.scale should be true. performing row scaling by default")
    X = t(scale(t(X), center=FALSE, scale=rowSums(X)))
    X=X*100
  }
  default.x.order = colnames(X)
  default.y.order = rownames(X) # Default order is reverse (starts from axis 0)
  #X = X[rev(1:dim(X)[1]),]
  X = melt(X)
  colnames(X) = c("Known", "Predicted", "Percentage")
  X$Known = factor(X$Known, levels=rev(unique(X$Known)));
  X$Predicted = as.factor(X$Predicted)
  p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low ="#edf8b1",   high = "#2c7fb8", limits=c(0, 100 ))+scale_size(range = c(1, max.size))+   theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic"))
  p = p + scale_x_discrete(limit = x.order) + scale_y_discrete(limit = y.order)
  print(p)
}

### Load and look at previously generated clusters

load("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/Shafer_Hypo_66k.Robj") # named hypo
hypo.zeb <- hypo
hypo.zeb@meta.data$dataset <- "Zebrafish" 
load("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_64k.Robj") #named hypo
hypo.ast <- hypo
hypo.ast@meta.data$dataset <- "Astyanax"

TSNEPlot(hypo.zeb, label.size = 5, do.label = T, no.legend = T)
TSNEPlot(hypo.ast, label.size = 5, do.label = T, no.legend = T)

### Set up genes used for calling Random Forest assignments (intersection of both var.gene datasets vectors)

genes.tfs <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/GO_DNA_binding.csv", header = F)
genes.tfs <- unique(as.character(genes.tfs[,2]))
genes.tfs <- genes.tfs[genes.tfs %in% rownames(hypo.zeb@data)]
# genes.tfs <- genes.tfs[genes.tfs %in% rownames(hypo.ast@data)]

genes.nps <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/GO_neuropeptide.csv", header = F)
genes.nps <- unique(as.character(genes.nps[,2]))
genes.nps <- genes.nps[genes.nps %in% rownames(hypo.zeb@data)]
genes.nps <- genes.nps[genes.nps %in% rownames(hypo.ast@data)]

######################################################################################
############### Predict the TFs which regulate neuropeptide expression ###############
######################################################################################

hypo <- hypo.zeb
hypo.zeb <- SetAllIdent(hypo, id = "orig.ident")
hypo.zeb <- SubsetData(hypo.zeb, max.cells.per.ident = 2000)

set.seed(123)

index <- sample(1:ncol(hypo.zeb@data), size = 0.7*ncol(hypo.zeb@data))

train <- as.data.frame(as.matrix(t(hypo.zeb@data[c(genes.tfs, "hcrt"),index])))
test <- as.data.frame(as.matrix(t(hypo.zeb@data[c(genes.tfs, "hcrt"),-index])))

hypo.train <- SubsetData(hypo.zeb, cells.use = row.names(hypo.zeb@meta.data)[index])
hypo.test <- SubsetData(hypo.zeb, cells.use = row.names(hypo.zeb@meta.data)[-index])

# Decision tree (singular)
library(rpart)
library(rattle)

rt <- rpart(hcrt ~ ., data = as.data.frame(as.matrix(t(hypo.train@data[c(genes.tfs, "hcrt"),]))))
dev.new()
fancyRpartPlot(rt)

test.pred.rtree.hcrt <- predict(rt, test)

hypo.test@data <- rbind(hypo.test@data, test.pred.rtree.hcrt)

RMSE.rtree <- sqrt(mean((test.pred.rtree-test$crhb)^2))
MAE.rtree <- mean(abs(test.pred.rtree-test$crhb))

# printcp(rt)
# rt$variable.importance

dev.new()
DotPlot(object = hypo.test, genes.plot = c(names(rt$variable.importance), "hcrt", "test.pred.rtree.hcrt"), plot.legend = TRUE, x.lab.rot = TRUE, group.by = "SubclusterType") + theme(legend.position = "right")


# linear regression

best.guess <- mean(train[,"crhb"])
RMSE.baseline <- sqrt(mean((best.guess-test[,"crhb"])^2))
MAE.baseline <- mean(abs(best.guess-test[,"crhb"]))

lin.reg <- lm(log(crhb+1) ~ ., data = train)

exp(lin.reg$coefficients["hopx"])

test.pred.lin <- exp(predict(lin.reg,test))-1
RMSE.lin.reg <- sqrt(mean((test.pred.lin-test$crhb)^2))
MAE.lin.reg <- mean(abs(test.pred.lin-test$crhb))



# Random forest

train <- train[,1:1248]
test <- test[,1:1248]


hypo <- hypo.zeb
hypo.zeb <- SetAllIdent(hypo, id = "orig.ident")
hypo.zeb <- SubsetData(hypo.zeb, max.cells.per.ident = 500)

set.seed(123)

index <- sample(1:ncol(hypo.zeb@data), size = 0.7*ncol(hypo.zeb@data))

genes.tfs <- genes.tfs[-grep("-", genes.tfs)]
genes.tfs <- genes.tfs[-grep(":", genes.tfs)]
train <- as.data.frame(as.matrix(t(hypo.zeb@data[c(genes.tfs, "crhb"),index])))
test <- as.data.frame(as.matrix(t(hypo.zeb@data[c(genes.tfs, "crhb"),-index])))


rf <- randomForest(crhb ~ ., data = train, importance = T, ntree = 500)

which.min(rf$mse)
plot(rf)

imp <- as.data.frame(sort(randomForest::importance(rf)[,1],decreasing = TRUE),optional = T)

test.pred.forest <- predict(rf, test)
RMSE.forest <- sqrt(mean((test.pred.forest-test$crhb)^2))
MAE.forest <- mean(abs(test.pred.forest-test$crhb))

dev.new()
DotPlot(object = hypo.zeb, genes.plot = c(row.names(imp)[1:20], "crhb"), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE) + theme(legend.position = "right")


######################################################################################
############## Predict the adult clusters that juvenile clusters map to ##############
######################################################################################

############### Predict the adult clusters that adult clusters map to ################

# Set up training and test sets (mutually exclusive lists), by 1) Getting at least 75% of each cluster (up to 500 cells)
# Returns roughly ~62% as training set, leaving ~38% as test set of the adult data OR 2) Picking a random subset (62% of data) (see below)

training.set <- c(); test.set <- c()
training.label <- c(); test.label <- c();

for (i in 1:length((levels(hypo.zeb@ident)))) {
  cells.in.clust <- WhichCells(hypo.zeb, levels(hypo.zeb@ident)[i]); # Pull cells from cluster
  n <- min(500, round(length(cells.in.clust)*0.75)) # Take the minimum between either 500 or 75% of the cluster
  train.temp <- cells.in.clust[sample(length(cells.in.clust))][1:n] # Sample n of the cell names
  test.temp <- setdiff(cells.in.clust, train.temp) # Calculate the set difference between cluster cell names and training set
  training.set <- c(training.set,train.temp)
  test.set <- c(test.set,test.temp) # Sum the current for loop into the training and test sets
  training.label <- hypo.zeb@data["galn", training.set]
  test.label <- hypo.zeb@data["galn", test.set]
}

# Collect the data matrix from the data slot of the seurat object, scale, and remove NAs to 0s

predictor_Data <- t(scale(t(as.matrix(hypo.zeb@data[genes.tfs,])), center = TRUE, scale = TRUE))
predictor_Data[is.na(predictor_Data)] <- 0

# Run randomForest algorithms on the training set, then predict the test set

tmp <- as.vector(table(training.label)) # number of cells selected for each cluster
sampsizes <- rep(min(tmp),length(tmp)) # the cluster with the minimum number of cells selected

rf_output <- randomForest(x = t(predictor_Data[,training.set]), 
                          y = factor(training.label), 
                          importance = TRUE, 
                          ntree = 501, # Try 1001: No change in prediction
                          proximity = TRUE, 
                          sampsize = sampsizes, 
                          keep.inbag = TRUE, 
                          replace = FALSE) 

Conf_OOB0 <- rf_output$confusion # Confusion matrix generated by randomForest

# Assign cells to cluster based on randomForest prediction

test.predict.prob <- predict(rf_output,t(predictor_Data[,test.set]), type = "prob") # Matrix of test cells and the predicted probabilities for each adult cluster
thresh <- 0.1 # The class with the maximum probability needs to have at least this margin
test.predict <- apply(test.predict.prob,1,function(x) if (max(x) > thresh){ which.max(x) } else {38}) # Assign each cell to only 1 cluster
Conf_test_test <- table(test.label,test.predict) # New Confusion Matrix (not sure why better, also screws up cluster numbers (+1))

# Generate a graph of the Confusion Matrix

pdf("ConfusionMatrixAdultToAdult.pdf", width = 14, height = 12)
plotConfusionMatrix(Conf_test_test, row.scale=TRUE, max.size = 12, xlab.use="Adult Actual", ylab.use="Adult Predicted")
plotConfusionMatrix(Conf_OOB0, row.scale=TRUE, max.size = 12, xlab.use="Adult Actual", ylab.use="Adult Predicted")
dev.off()

############### Predict the adult clusters that Juvenile clusters map to ################

# Collect the data matrix from the data slot of the seurat object, scale to adult data, and remove NAs to 0s

Juvenile.rf <- as.matrix(raj@data[genes.use,])
Juvenile.rf <- t(scale(t(Juvenile.rf), center=rowMeans(as.matrix(hypo@data[genes.use,])), scale=TRUE))
Juvenile.rf[is.na(Juvenile.rf)] <- 0

#Adult.predict.prob = predict(rf_output,t(Adult.rf), type="prob")
#Adult.predict.vote = predict(rf_output,t(Adult.rf))
#max.prob.adult = apply(Adult.predict,1,max)
#confident_assignments = which(max.prob.adult > 0.2) # test
#Conf_test = table(Adult.ident[confident_assignments], Adult.predict.vote[confident_assignments])
#Output probabilities 

# Predict juvenile data as adult clusters

Juvenile.ident <- raj@ident

Juvenile.predict.prob <- predict(rf_output,t(Juvenile.rf), type = "prob")

# Assign cells to cluster based on randomForest prediction

thresh <- 0.1 # The class with the maximum probability needs to have at least this margin
Juvenile.predict <- apply(Juvenile.predict.prob,1,function(x) if (max(x) > thresh){ which.max(x) } else {38})
Conf_test_predict <- table(Juvenile.ident, Juvenile.predict) 

# Generate a graph of the Confusion Matrix

pdf("ConfusionMatrixJuvenileToAdulte.pdf", width = 14, height = 13)
plotConfusionMatrix(Conf_test_predict, row.scale=TRUE, max.size = 12, xlab.use="Adult", ylab.use="Juvenile")
dev.off()

######################################################################################
############## Predict the juvenile clusters that adult clusters map to ##############
######################################################################################


############ Predict the juvenile clusters that juvenile clusters map to #############

# Set up training and test sets (mutually exclusive lists), by 1) Getting at least 75% of each cluster (up to 500 cells)
# Returns roughly ~74% as training set, leaving ~38% as test set of the adult data OR 2) Picking a random subset (62% of data) (see below)

training.set <- c(); test.set <- c()
training.label <- c(); test.label <- c();
for (i in as.numeric(levels(raj@ident))){
  cells.in.clust <- WhichCells(raj,i); # Identify cell names from cluster
  n <- min(500, round(length(cells.in.clust)*0.75)) # Take the minimum between either 500 or 75% of the cluster
  train.temp <- cells.in.clust[sample(length(cells.in.clust))][1:n] # Sample 500 of the cell names
  test.temp <- setdiff(cells.in.clust, train.temp) # Calculate the set difference between cluster cell names and training set
  training.set <- c(training.set,train.temp); test.set=c(test.set,test.temp) # Sum the current for loop into the training and test sets
  training.label <- c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp))); # Same for labels
}

# Collect the data matrix from the data slot of the seurat object, scale, and remove NAs to 0s

predictor_Data <- t(scale(t(as.matrix(raj@data[genes.use,])), center = TRUE, scale = TRUE))
predictor_Data[is.na(predictor_Data)] <- 0

# Run randomForest algorithms on the training set, then predict the test set

tmp <- as.vector(table(training.label)) # number of cells selected for each cluster
sampsizes <- rep(min(tmp),length(tmp)) # the cluster with the minimum number of cells selected

rf_output <- randomForest(x = t(predictor_Data[,training.set]), 
                          y = factor(training.label), 
                          importance = TRUE, 
                          ntree = 1001, # Try 1001: No change in prediction
                          proximity = TRUE, 
                          sampsize = sampsizes, 
                          keep.inbag = TRUE, 
                          replace = FALSE) 

Conf_OOB0 <- rf_output$confusion # Confusion matrix generated by randomForest

# Assign cells to cluster based on randomForest prediction

test.predict.prob <- predict(rf_output,t(predictor_Data[,test.set]), type = "prob") # Matrix of test cells and the predicted probabilities for each adult cluster
thresh <- 0.16 # The class with the maximum probability needs to have at least this margin
test.predict <- apply(test.predict.prob,1,function(x) if (max(x) > thresh){ which.max(x) } else {13}) # Assign each cell to only 1 cluster
Conf_test_test <- table(test.label,test.predict) # New Confusion Matrix (not sure why better, also screws up cluster numbers (+1))

# Generate a graph of the Confusion Matrix

pdf("ConfusionMatrixJuvenileToJuvenile.pdf", width = 14, height = 12)
plotConfusionMatrix(Conf_test_test, row.scale=TRUE, max.size = 12, xlab.use="Juvenile Actual", ylab.use="Juvenile Predicted")
# plotConfusionMatrix(Conf_OOB0, row.scale=TRUE, max.size = 12, xlab.use="Adult Actual", ylab.use="Adult Predicted")
dev.off()

############ Predict the juvenile clusters that juvenile clusters map to #############

# Collect the data matrix from the data slot of the seurat object, scale to juvenile data, and remove NAs to 0s

Adult.rf <- as.matrix(hypo@data[genes.use,])
Adult.rf <- t(scale(t(Adult.rf), center=rowMeans(as.matrix(raj@data[genes.use,])), scale=TRUE))
Adult.rf[is.na(Adult.rf)] <- 0

#Adult.predict.prob = predict(rf_output,t(Adult.rf), type="prob")
#Adult.predict.vote = predict(rf_output,t(Adult.rf))
#max.prob.adult = apply(Adult.predict,1,max)
#confident_assignments = which(max.prob.adult > 0.2) # test
#Conf_test = table(Adult.ident[confident_assignments], Adult.predict.vote[confident_assignments])
#Output probabilities 

# Predict Adult data as adult clusters

Adult.ident <- hypo@ident

Adult.predict.prob <- predict(rf_output,t(Adult.rf), type = "prob")

# Assign cells to cluster based on randomForest prediction

thresh <- 0.16 # The class with the maximum probability needs to have at least this margin
Adult.predict <- apply(Adult.predict.prob,1,function(x) if (max(x) > thresh){ which.max(x) } else {13})
Conf_test_predict <- table(Adult.ident, Adult.predict) 

# Re-order columns to make graph look better

order <- c("2", "6", "1", "7", "4", "9", "5", "10", "3", "8", "12", "11", "13")

# Generate a graph of the Confusion Matrix

pdf("ConfusionMatrixAdultToAJuvenile.pdf", width = 14, height = 13)
plotConfusionMatrix(Conf_test_predict, x.order = order, row.scale=TRUE, max.size = 12, xlab.use="Juvenile", ylab.use="Adult") 
dev.off()


######################################################################################
######################################################################################
######################################################################################




# # 2) Picking a random subset (62% of data)
# 
# set.seed(624)
# samp <- sample(nrow(hypo@meta.data), 0.62 * nrow(hypo@meta.data))
# nrow <- c(0:nrow(hypo@meta.data))
# samp2 <- subset(nrow, !(nrow %in% samp))
# train <- hypo@data[,samp]
# test <- hypo@data[,samp2]

load("/Users/ShristiPandey/Dropbox/SingleCellAnalysis/ForPaper/10xadult/FinalAdultObj.RObj")


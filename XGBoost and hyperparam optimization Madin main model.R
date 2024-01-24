
rm(list=ls()) 

library(xgboost)
library(ParBayesianOptimization)
library(mlbench)
library(dplyr)
library(skimr)
library(recipes)
library(resample)
library(ggplot2)
library(tidyr)
library(rsample)
library(caret)

##LOAD DATAFRAME
genes.motility.counts.dbb = read.csv("C:/CU_Boulder/Chemotaxis/FinalAnalysis/Chemotaxis_motility_Madin.csv", header = T)
length(unique(genes.motility.counts.dbb$Gene))#21 GENES
length(unique(genes.motility.counts.dbb$Genome))#1225 GENOMES
genes.motility.counts.dbb = genes.motility.counts.dbb[,c(2,3,5,7)]
head(genes.motility.counts.dbb)

##ADD TAXONOMY TO DATAFRAME
meta.gtdb = read.csv("C:/CU_Boulder/Chemotaxis/bac120_metadata_r207.tsv", sep = "\t")
head(meta.gtdb)

##SUBSET METADATA BY GENOMES OF INTEREST AND EXTRACT TAXONOMIC INFORMATION
meta.gtdb.genomes = meta.gtdb[meta.gtdb$gtdb_genome_representative %in% genes.motility.counts.dbb$Genome,]
head(meta.gtdb.genomes)
length(unique(meta.gtdb.genomes$gtdb_genome_representative))

meta.gtdb.genomes.phyla = meta.gtdb.genomes[,colnames(meta.gtdb.genomes) %in% c("gtdb_genome_representative", "gtdb_taxonomy")]
head(meta.gtdb.genomes.phyla)

library(stringr)
meta.gtdb.genomes.phyla = cbind(meta.gtdb.genomes.phyla, str_split_fixed(meta.gtdb.genomes.phyla$gtdb_taxonomy, ";", 7))
meta.gtdb.genomes.phyla= meta.gtdb.genomes.phyla[,c(1,4,5)]
colnames(meta.gtdb.genomes.phyla) = c("Genome", "Phylum", "Class")
meta.gtdb.genomes.phyla$Phylum = gsub("p__", "", meta.gtdb.genomes.phyla$Phylum)
meta.gtdb.genomes.phyla$Class = gsub("c__", "", meta.gtdb.genomes.phyla$Class)
meta.gtdb.genomes.phyla$Phylum = gsub("Firmicutes_A", "Firmicutes", meta.gtdb.genomes.phyla$Phylum)
meta.gtdb.genomes.phyla$Phylum = gsub("Firmicutes_C", "Firmicutes", meta.gtdb.genomes.phyla$Phylum)
meta.gtdb.genomes.phyla$Phylum = gsub("Firmicutes_B", "Firmicutes", meta.gtdb.genomes.phyla$Phylum)
meta.gtdb.genomes.phyla$Phylum = gsub("Firmicutes_E", "Firmicutes", meta.gtdb.genomes.phyla$Phylum)
meta.gtdb.genomes.phyla$Phylum = gsub("Firmicutes_F", "Firmicutes", meta.gtdb.genomes.phyla$Phylum)
meta.gtdb.genomes.phyla$Phylum = gsub("Desulfobacterota_I", "Desulfobacterota", meta.gtdb.genomes.phyla$Phylum)
head(meta.gtdb.genomes.phyla)

genes.motility.counts.dbb.taxonomy = merge(genes.motility.counts.dbb, meta.gtdb.genomes.phyla, by = "Genome")
genes.motility.counts.dbb.taxonomy = unique(genes.motility.counts.dbb.taxonomy)
#genes.motility.counts.dbb.taxonomy = unique(genes.motility.counts.dbb.taxonomy[,c(1,2,5,6)])
head(genes.motility.counts.dbb.taxonomy)


##CREATE STRATIFIED SAMPLING TRAINING AND TEST SETS
dat1 = genes.motility.counts.dbb.taxonomy %>% 
  tidyr::pivot_wider(names_from = Gene, values_from = Presence)
model_df <- dat1
model_df$Strata = paste(model_df$Phylum, model_df$Motility)
model_df$Motility = as.numeric(as.factor(model_df$Motility))-1

set.seed(123)
splits <- initial_split(model_df, prop = 0.7)#, strata = c("Strata"))
train_df <- as.data.frame(as.matrix(training(splits)))
train_df = train_df[,c(2,5:25)]
train_labels <- as.numeric(train_df$Motility)
test_df <- as.data.frame(as.matrix(testing(splits)))
test_df = test_df[,c(2,5:25)]
test_labels <- as.numeric(test_df$Motility)



##START HYPERPARAM OPTIMIZATION##

#TRAINING MATRIX
X <- xgb.DMatrix(data = data.matrix(train_df[,!names(train_df) %in% c("Motility")])-1, label = data.matrix(train_labels)) 

#TEST MATRIX
test_df_xgb = xgb.DMatrix(data = data.matrix(test_df[,!names(test_df) %in% c("Motility")])-1, label = data.matrix(test_labels))

#TARGET VARIABLE
y <- train_df %>%
  pull(Motility)
#y = as.numeric(as.factor(y))-1

#CROSS VALIDATION FOLDS
folds <- list(fold1 = as.integer(seq(1, nrow(X), by = 3)),
              fold2 = as.integer(seq(2, nrow(X), by = 3)),
              fold3 = as.integer(seq(3, nrow(X), by = 3)))

#CREATE OBJECTIVE FUNCTION: Function must take the hyper-parameters as inputs
obj_func <- function(max_depth, min_child_weight, subsample) {
  
  Pars <- list(
    booster = "gbtree"
    , eta = 0.01
    , max_depth = max_depth
    , min_child_weight = min_child_weight
    , subsample = subsample
    , objective = "binary:logistic"
    , eval_metric = c("auc", "rmse")
  )
  xgbcv <- xgb.cv(
    params = Pars
    , data = X
    , nround = 100
    , folds = folds
    , prediction = TRUE
    , showsd = TRUE
    , early_stopping_rounds = 5
    , maximize = TRUE
    , verbose = 0
  )
  
  return(
    list(
      Score = max(xgbcv$evaluation_log$test_auc_mean)
      , nrounds = xgbcv$best_iteration
    )
  )
}

bounds <- list(max_depth = c(1L, 10L),
               min_child_weight = c(1, 100),
               subsample = c(0, 1))

#RUN OPTIMISER TO FIND SET OF OPTIMAL HYPERPARAMETERS
set.seed(15)
bayes_out <- bayesOpt(FUN = obj_func, bounds = bounds, initPoints = length(bounds) + 2, iters.n = 10, iters.k = 2, acq = "ei", gsPoints = 10, 
                      parallel = F, verbose = 0)
bayes_out


# Show relevant columns from the summary object 
bayes_out$scoreSummary[1:15, ]

# Get best parameters
data.frame(getBestPars(bayes_out))
best.pars = getBestPars(bayes_out)
best.pars

# EVALUATE IF THE BAYESIAN OPTIMIZATION CONVERGES INTO A STABLE "UTILITY" OPTIMUM
optObjSimp <- addIterations(bayes_out,10,verbose=FALSE)
plot(optObjSimp)



##FINALLY, FIT MODEL WITH THE OPTIMIZED PARAMETERS
# Combine best params with base params
opt_params <- append(list(booster = "gbtree", 
                          objective = "binary:logistic", 
                          eval_metric = c("auc", "rmse")), 
                     best.pars)

# Run cross validation 
xgbcv <- xgb.cv(
  params = opt_params
  , data = X
  , nround = 100
  , folds = folds
  , prediction = TRUE
  , showsd = TRUE
  , early_stopping_rounds = 5
  , maximize = TRUE
  , verbose = 0)

# Get optimal number of rounds
nrounds = xgbcv$best_iteration
nrounds

# Fit a xgb model
mdl <- xgboost(data = X, #label = y, 
               params = opt_params, 
               maximize = T, 
               early_stopping_rounds = 5, 
               nrounds = nrounds, 
               verbose = 0)
mdl

# Evaluate performance 
actuals <- test_df$Motility
xgpred <- predict(mdl, test_df_xgb)
xgpred.coded = ifelse(xgpred > 0.5, 1, 0)


library(Metrics)
perf.dbb = data.frame(Observed = actuals, Predicted = xgpred.coded)
perf.dbb$Errors = abs(as.numeric(perf.dbb$Observed) - as.numeric(perf.dbb$Predicted))
total.error = sum(perf.dbb$Errors)/nrow(perf.dbb)
total.error #ACCURACY = 1 - 0.038 = 0.962 ACCURACY STRATA = 0.949

mdl #TRAINING SET AUC = 0.98
xgbcv #TEST SET AUC = 0.98

perf.dbb$Consensus = ifelse(perf.dbb$Errors == 0, "Correct", "Incorrect")
perf.dbb$Motility = ifelse(perf.dbb$Observed == 1, "Yes", "No")
head(perf.dbb)

pred.dbb.plot = perf.dbb %>% dplyr::group_by(Consensus, Motility) %>% dplyr::count() 
pred.dbb.plot$Motility = factor(pred.dbb.plot$Motility, levels = c("Yes", "No"))
pred.dbb.plot

accuracy = ggplot(pred.dbb.plot, aes(x = Motility, y = n)) + geom_bar(aes(fill = Consensus), position="fill", stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(size = 16), axis.text.y=element_text(size=16), legend.position = "top", 
        legend.text = element_text(size = 14), legend.title = element_blank(), axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size =18, margin = margin(r = 10)), plot.title = element_text(size=18, face = "bold", 
  colour = "black", vjust = 2.5), plot.margin = margin(1, 1, 1, 1, "cm")) + xlab("Motility (observed)") + 
  scale_fill_manual(values=c("grey", "black")) + ylab("Proportion of correct predictions")


tiff("C:/CU_Boulder/Chemotaxis/FinalAnalysis/Model_performance.tiff", units="in", width=4, height=5, res=300)
accuracy
dev.off()

###FEATURE IMPORTANCE
# Compute feature importance matrix
importance_matrix = xgb.importance(colnames(X), model = mdl)
xgb.plot.importance(importance_matrix)
importance_matrix

xgb.ggplot.importance(importance_matrix) + theme_bw() +
  theme(axis.text.x = element_text(size = 16), axis.text.y=element_text(size=12), legend.position = "top", legend.text = element_text(size = 16),
        legend.title = element_text(size = 16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size =18, margin = margin(r = 10)), 
        plot.title = element_text(size=18, face = "bold", colour = "black", vjust = 2.5), 
        plot.margin = margin(1, 1, 1, 1, "cm"))  


#saveRDS(mdl, file = "C:/CU_Boulder/Chemotaxis/boosted_regression_model_flagellar_motility_binary_Madin.rda")



###PCA TO EVALUATE RELATIONSHIPS BETWEEN GENES AND MOTILITY
library(tidyr)
library(ggfortify)
library(factoextra)

df <- dat1[,c(5:25)]
pca_res <- prcomp(df, scale. = FALSE)
fviz_eig(pca_res)
fviz_pca_ind(pca_res)

dat1$Motility = factor(dat1$Motility, levels = c("Yes", "No"))
pca.stratified = ggplot(data = dat1, aes(x = as.data.frame(pca_res$x)$PC1, y = as.data.frame(pca_res$x)$PC2, color = Motility)) + theme_bw() + 
  theme(axis.text.x = element_text(size = 16), axis.text.y=element_text(size=16), legend.position = "top", legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 18), axis.title.y = element_text(size =18), 
        plot.title = element_text(size=18, face = "bold", colour = "black", vjust = 2.5), plot.margin = margin(1, 1, 1, 1, "cm")) + geom_jitter(width = 0.75, height = 0.75, size = 5, alpha = 0.75) +
  geom_segment(data = as.data.frame(pca_res$rotation), aes(x=0, y=0, xend=PC1*3.75, yend=PC2*1.75), arrow = arrow(length=unit(0.5, 'cm')), color = "black") +
  geom_text(data = as.data.frame(pca_res$rotation), aes(x = PC1*5.25, 
  y = PC2*2.5, label = rownames(pca_res$rotation)), colour = "black", size = 5) + xlab("PC1 (65.4%)") + ylab("PC2 (4.1%)") + scale_color_manual(values = c("#721422", "#B8B799"))

tiff("C:/CU_Boulder/Chemotaxis/FinalAnalysis/PCA_final.tiff", units="in", width=6, height=5.5, res=300)
pca.stratified
dev.off()


pca.stratified.nolabs = ggplot(data = dat1, aes(x = as.data.frame(pca_res$x)$PC1, y = as.data.frame(pca_res$x)$PC2, color = Motility)) + theme_bw() + 
  theme(axis.text.x = element_text(size = 16), axis.text.y=element_text(size=16), legend.position = "top", legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 18), axis.title.y = element_text(size =18), 
        plot.title = element_text(size=18, face = "bold", colour = "black", vjust = 2.5), plot.margin = margin(1, 1, 1, 1, "cm")) + geom_jitter(width = 0.75, height = 0.75, size = 5, alpha = 0.75) +
  geom_segment(data = as.data.frame(pca_res$rotation), aes(x=0, y=0, xend=PC1*3.75, yend=PC2*1.75), arrow = arrow(length=unit(0.5, 'cm')), color = "black") +
  xlab("PC1 (65.4%)") + ylab("PC2 (4.1%)") + scale_color_manual(values = c("#B8B799", "#721422"))

tiff("C:/CU_Boulder/Chemotaxis/FinalAnalysis/PCA_final_stratified_nolabs.tiff", units="in", width=6, height=5.5, res=300)
pca.stratified.nolabs
dev.off()

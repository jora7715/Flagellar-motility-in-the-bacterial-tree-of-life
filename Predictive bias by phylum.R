####TESTING BIAS IN THE PREDICTIVE ACCURACY ACROSS PHYLA
rm(list=ls()) 

####1 MADIN ET AL 2020
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

madin.motility.phyla = genes.motility.counts.dbb.taxonomy
unique(madin.motility.phyla$Phylum)
head(madin.motility.phyla)

##RESTRICT TO THE GENOMES IN THE TEST SET (LOAD TEST_DF FROM MAIN GENOMIC MODEL)
test_df

madin.motility.phyla.test = madin.motility.phyla[madin.motility.phyla$Genome %in% test_df$Genome,]
head(madin.motility.phyla.test)


##SUBSET BY PHYLUM OF INTEREST
madin.motility.proteos = subset(madin.motility.phyla.test, Phylum == "Proteobacteria")
head(madin.motility.proteos)
madin.motility.actinos = subset(madin.motility.phyla.test, Phylum == "Actinobacteriota")
head(madin.motility.actinos)
madin.motility.firm = subset(madin.motility.phyla.test, Phylum == "Firmicutes")
head(madin.motility.firm)
madin.motility.bact = subset(madin.motility.phyla.test, Phylum == "Bacteroidota")
head(madin.motility.bact)
madin.motility.acido = subset(madin.motility.phyla.test, Phylum == "Acidobacteriota")
head(madin.motility.acido)


##MAKE PREDICTIONS FOR EACH PHYLUM
df2 <- madin.motility.proteos
head(df2)
length(unique(df2$Genome))#444 unique genomes

library(tidyr)
dat1 = df2 %>% 
  pivot_wider(names_from = "Gene", values_from = "Presence")
head(dat1)


#GET BOOSTED REGRESSION MODEL FROM "C:/CU_Boulder/Chemotaxis/XGBoost and hyperparam optimization.r"
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

mdl = readRDS(file = "C:/CU_Boulder/Chemotaxis/FinalAnalysis/boosted_regression_model_flagellar_motility_binary_Madin.rda")
mdl

#LOAD ORDER OF FLAGELLAR GENES - THEY NEED TO MATCH!
order.genes = mdl$feature_names

query_dbb_xgb = dat1[,-c(1:3)]
#query_dbb_xgb$Flagellin_IN = rep(0, nrow(query_dbb_xgb))

head(query_dbb_xgb)

query_dbb_xgb = query_dbb_xgb[order.genes]

order.genes == colnames(query_dbb_xgb)#ALL TRUE, WE CAN GO AHEAD

query_dbb_xgb = as.matrix(as.data.frame(query_dbb_xgb))
head(query_dbb_xgb)

xgpred <- predict(mdl, query_dbb_xgb)
xgpred

dat1.complete = as.data.frame(cbind(xgpred, dat1))
dat1.complete$Predicted.motility = ifelse(dat1.complete$xgpred > 0.5,"Yes", "No")
head(dat1.complete)

#write.csv(dat1.complete, "C:/CU_Boulder/Chemotaxis/FinalAnalysis/Predicted_motility_Acidobacteriota.csv")


###PLOT RESULT
dat1.complete$Consensus = ifelse(dat1.complete$Motility == dat1.complete$Predicted.motility, "Yes", "No") 
colnames(dat1.complete)
head(dat1.complete)

consensus.summary = dat1.complete %>% group_by(Motility, Consensus) %>% count()
consensus.summary

consensus.summary$Motility = factor(consensus.summary$Motility, levels = c("Yes", "No"))
consensus.summary$Consensus = ifelse(consensus.summary$Consensus == "Yes", "Correct", "Incorrect")
consensus.summary

p = ggplot(consensus.summary, aes(x = Motility, y = n)) + geom_bar(aes(fill = Consensus), position="fill", stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(size = 16), axis.text.y=element_text(size=16), legend.position = "top", 
        legend.text = element_text(size = 14), legend.title = element_blank(), axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size =18, margin = margin(r = 10)), plot.title = element_text(size=18, face = "bold", 
                                                                                                  colour = "black", vjust = 2.5), plot.margin = margin(1, 1, 1, 1, "cm")) + xlab("Motility (observed)") + 
  scale_fill_manual(values=c("grey", "black")) + ylab("Proportion of correct predictions")
p

tiff("C:/CU_Boulder/Chemotaxis/FinalAnalysis/Model_performance_proteobacteria.tiff", units="in", width=4, height=5, res=300)
p
dev.off()

# setwd('/home/labs/elinav/uria/rafa_code/final_notebooks')
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)

# util func ####
run.lmer.feature <- function(f, table){
  form <- paste0(f, " ~ rDay * rGroup  + (1  | Participant)")
  m <- lmer(form, data = table)
  aov.res <- as.data.frame(anova(m, type ='II'))
  list(Group.pvalue = aov.res['rGroup','Pr(>F)'], interaction.pvalue = aov.res['rDay:rGroup','Pr(>F)'], conv.mess = m@optinfo$conv$lme4$messages)
}

# read data ####
# data preview on the next cell
raw.data <- read.table("./Suez2018.txt",header = T,  sep = "\t")
raw.data$rGroup <- factor(raw.data$rGroup, levels = c('CTR','PBX','FMT'))
raw.data$Participant <- factor(raw.data$Participant)
raw.data$rDay <- as.double(raw.data$rDay)
meta.cols <- c(1:6)
data.cols <- unlist(colnames(raw.data)[-c(1:6)])




# log scale the relative abundance data 
raw.data.logscaled <- data.frame(raw.data)
raw.data.logscaled[,data.cols] <- log( 100 * raw.data.logscaled[,data.cols])
# raw.data.logscaled[,data.cols] <- scale(raw.data.logscaled[,data.cols], center = F)

# apply univariate p-values (see run.lmer.feature function)
ra.lmer <- sapply(data.cols, FUN = function(f){
  run.lmer.feature(f,raw.data.logscaled)
}) %>% t
ra.lmer <- as.data.frame(ra.lmer)
ra.lmer$Group.padj <- p.adjust(ra.lmer$Group.pvalue,method = 'BH') # FDR for group factor
ra.lmer$interaction.padj <- p.adjust(ra.lmer$interaction.pvalue,method = 'BH') # FDR for interaction



p.list2 <- read.csv("./tcam_post_abx_features_Sue2018.csv", header = T) # read the pruned features list
ra.lmer.pruned2 <- ra.lmer[rownames(ra.lmer) %in% as.character(unlist(p.list2)), ] # filter the pvalue table 
ra.lmer.pruned2$interaction.padj <- p.adjust(ra.lmer.pruned2$interaction.pvalue
                                             ,method = 'BH') # recompute FDR for interaction
ra.lmer.pruned2$Group.padj <- p.adjust(ra.lmer.pruned2$Group.pvalue
                                       ,method = 'BH') # recompute FDR for group factor
ra.lmer.pruned2 <- data.frame(lapply(ra.lmer.pruned2[,-c(3)], FUN = unlist)) # make a data frame


write.table(x=as.matrix(ra.lmer.pruned2), file = 'Suez2018_pruned_ra_lmer_results.txt', sep ="\t", quote = F, row.names = T)
write.table(x=as.matrix(ra.lmer), file = 'Suez2018_ra_lmer_results.txt', sep ="\t", quote = F, row.names = T)
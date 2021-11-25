library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)

# util func ####
run.lmer.feature <- function(f, table){
  form <- paste0(f, " ~ Week * Group  + (1  | Participant)")
  m <- lmer(form, data = table)
  aov.res <- as.data.frame(anova(m, type ='II'))
  
  list(Group.pvalue = aov.res['Group','Pr(>F)'], interaction.pvalue = aov.res['Week:Group','Pr(>F)'], conv.mess = m@optinfo$conv$lme4$messages)
}

# read data ####
# data preview on the next cell

raw.data <- read.table("Deehan2020.tsv",header = T,  sep = "\t", comment.char = "")
meta <- read.table("metadata_Deehan2020.tsv", header=TRUE, sep="\t")
colnames(meta) <- c("Participant", "Group")
raw.data <- merge(meta, raw.data, by="Participant")
meta.cols <- c(1:3)
data.cols <- unlist(colnames(raw.data)[-c(1:3)])


# log scale the relative abundance data 
raw.data.logscaled <- data.frame(raw.data)
raw.data.logscaled[,data.cols] <- log(raw.data.logscaled[,data.cols])


# apply univariate p-values (see run.lmer.feature function)
ra.lmer <- sapply(data.cols, FUN = function(f){
  run.lmer.feature(f,raw.data.logscaled)
}) %>% t
ra.lmer <- as.data.frame(ra.lmer)
ra.lmer$Group.padj <- p.adjust(ra.lmer$Group.pvalue,method = 'fdr') # FDR for group factor
ra.lmer$interaction.padj <- p.adjust(ra.lmer$interaction.pvalue,method = 'fdr') # FDR for interaction



p.list2 <- read.csv("Deehan2020_pruned_features.csv", header = T) # read the pruned features list
p.list2$X0 <- gsub(";", ".", p.list2$X0)
p.list2$X0 <- gsub("\\[", ".", p.list2$X0)
p.list2$X0 <- gsub("\\]", ".", p.list2$X0)
p.list2$X0 <- gsub("-", ".", p.list2$X0)
ra.lmer.pruned2 <- ra.lmer[rownames(ra.lmer) %in% as.character(unlist(p.list2)), ] # filter the pvalue table 
ra.lmer.pruned2$interaction.padj <- p.adjust(ra.lmer.pruned2$interaction.pvalue
                                             ,method = 'fdr') # recompute FDR for interaction
ra.lmer.pruned2$Group.padj <- p.adjust(ra.lmer.pruned2$Group.pvalue
                                       ,method = 'fdr') # recompute FDR for group factor
ra.lmer.pruned2 <- data.frame(lapply(ra.lmer.pruned2[,-c(3)], FUN = unlist)) # make a data frame
write.table(ra.lmer.pruned2, "Deehan2020_ra.lmer.pruned.tsv", row.names=TRUE, sep="\t", quote=FALSE)


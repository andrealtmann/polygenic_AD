#load libraries
require(ADNIMERGE)
require(data.table)

#require(lme4)
#require(survival)

#source helper functions
source("./helper_functions.R")

polygenic_scores <- fread("../data/polygenic_scores.csv", data.table=F)
apoe.info2 <- merge(polygenic_scores, desikanlab, by="RID", all=F)

#selected the name of the score to be used
#(options: PHS, PRScs_auto, PRS1, PRS2)
if (!exists("mkroc"))
  score.use <- "PRS2"

#list the variables that are used as confounding variables
#confound <- paste(c("PTGENDER", "PTEDUCAT"),collapse=" + ")
confound <- paste(c("PTGENDER", "PTEDUCAT", paste("PCA",1:5,sep="")),collapse=" + ")

#this removes subjects who contributed to ADGC
#(set to true if running AD-PRS based scores)
exclude.adgc <- T

#this removes subjects with Ashkenazi Jewish ancestry
#they form an outlier group with PCA1 > 0.06, verified by SNPweights
exclude.aj <- F

if(exclude.adgc)
  apoe.info2 <- subset(apoe.info2, EXCLUDE==0)
if (exclude.aj)
  apoe.info2 <- subset(apoe.info2, PCA1<=0.06)

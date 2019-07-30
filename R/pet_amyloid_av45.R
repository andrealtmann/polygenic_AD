source("./initialize.R")

### AV45 for all cortical regions ###

#ROI used for normalization of PET intensities
norm <- "WHOLECEREBELLUM_SUVR"

#select a subset of the data with amyloid image at any visit
amybase <- subset(adnimerge, subset=!is.na(AV45))
tmp <- unique(amybase$RID)

#identify first visit of subject for AV45 scan
vis.sel <- t(sapply(tmp, function(rid){
  tmp <- subset(amybase, subset=RID==rid)
  idx <- order(tmp$Years.bl)[1]
  return(tmp[idx,c("RID","VISCODE")])
}))
vis.sel <- data.frame(vis.sel)
amybase2 <- merge(amybase, vis.sel, by=c("RID","VISCODE"), all=F)

#merge with detailed UCBerkeley data:
amydata     <- merge(amybase2, ucberkeleyav45, by=c("RID", "VISCODE"), all=F)
#add polygenic info
amydata2    <- merge(amydata, apoe.info2, by="RID", all=F)

#restrict to white non hispanic subjects
amydata2.wn <- subset(amydata2, subset=PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino")

###average data from left and right hemisphere

##Note: these have the _SUVR postfix, but they are not SUVRs
target.rois.rh <- colnames(amydata2.wn)[grep("CTX_RH_.*_SUVR",colnames(amydata2.wn))]
target.rois.lh <- colnames(amydata2.wn)[grep("CTX_LH_.*_SUVR",colnames(amydata2.wn))]

lbase <- gsub("CTX_LH_","",target.rois.lh)
rbase <- gsub("CTX_RH_","",target.rois.rh)

discord <- sum(lbase!=rbase)
stopifnot(discord==0)

summed <- amydata2.wn[,target.rois.lh] + amydata2.wn[,target.rois.rh]
colnames(summed) <- lbase
amydata2.wn <- data.frame(amydata2.wn, summed)

target.rois <- c("COMPOSITE_SUVR", setdiff(lbase,"UNKNOWN_SUVR"))

#run the model for every ROI (and the composite)

av45.full <- t(sapply(target.rois, function(roi){
  #basic model only adjusting for APOE4 carrier status
  f0 <- paste( "eval(", roi,"/", norm, ")","~ eval(APOE4>0)*1 + eval(AGE+Years.bl) +", confound)
  #adding the scaled polygenic score
  f1 <- paste(f0, "+ scale(", score.use, ")")
  #adding APOE4 burden
  f2 <- paste(f0, "+ APOE4")

  m0 <- lm(as.formula(f0), data=amydata2.wn)
  m1 <- lm(as.formula(f1), data=amydata2.wn)
  m2 <- lm(as.formula(f2), data=amydata2.wn)
  base_score   <- anova(m0, m1)[2,6]
  simple_score <- anova(m0, m2)[2,6]

  #like f0, but adjusting for APOE4 locus
  g0 <- paste( "eval(", roi,"/", norm, ")","~ APOE4 + rs7412 + eval(AGE+Years.bl) + ", confound)
  g1 <- paste(g0, "+ scale(", score.use, ")")

  n0 <- lm(as.formula(g0), data=amydata2.wn)
  n1 <- lm(as.formula(g1), data=amydata2.wn)

  full_score <- anova(n0, n1)[2,6]
  full_meta <- metaAnalyze(amydata2.wn, g1, c("CN","MCI","Dementia"), paste("scale(", score.use,")",sep=""))

  return(c(base_score, full_score, simple_score, full_meta[4,4]))
}))

rownames(av45.full) <- gsub("_SUVR","",rownames(av45.full))
colnames(av45.full) <- c(paste(score.use, " (",c("Base","APOE"), ")", sep=""), "APOEe44", paste(score.use,"(APOE) meta"))

print(av45.full)

### END AV45  ###

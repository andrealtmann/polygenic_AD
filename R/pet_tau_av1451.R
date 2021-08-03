source("./initialize.R")

### Tau imaging AV1451 ###

#ROI to be used for normalization
norm <- "INFERIOR_CEREBGM_SUVR"

#add tau data to adnimerge
taubase <- merge(adnimerge, ucberkeleyav1451_pvc, by=c("RID","VISCODE"),all=F)

#add polygenic data to adnimerge
taubase2 <- merge(taubase, apoe.info2, by="RID", all=F)
tmp <- unique(taubase2$RID)

#find first visit with tau imaging
vis.sel <- t(sapply(tmp, function(rid){
  tmp <- subset(taubase2, subset=RID==rid)
  idx <- order(tmp$Years.bl)[1]
  return(tmp[idx,c("RID","VISCODE")])
}))
vis.sel <- data.frame(vis.sel)

#select a subset of the data with tau image at any visit
taubase3 <- merge(taubase2, vis.sel, by=c("RID","VISCODE"), all=F)

#limit to white non-hisp/latino participants
taudata.wn <- subset(taubase3, subset=PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino")
#use this to additionally limit to people w/o Dementia diagnosis
#taudata.wn <- subset(taubase3, subset=PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino" & DX != "Dementia")


##avg left and right hemisphere (just like in the case with amyloid)
target.rois.lh <- colnames(taudata.wn)[grep("CTX_LH_.*_SUVR",colnames(taudata.wn))]
target.rois.rh <- colnames(taudata.wn)[grep("CTX_RH_.*_SUVR",colnames(taudata.wn))]

lbase <- gsub("CTX_LH_","",target.rois.lh)
rbase <- gsub("CTX_RH_","",target.rois.rh)
discord <- sum(lbase!=rbase)
stopifnot(discord==0)
summed <- taudata.wn[,target.rois.lh] + taudata.wn[,target.rois.rh]
colnames(summed) <- lbase
taudata.wn <- data.frame(taudata.wn, summed)


#target regions of interest
target.rois <- c("BRAAK12_SUVR", "BRAAK34_SUVR", "BRAAK56_SUVR", setdiff(lbase,"UNKNOWN_SUVR"))
#target.rois <- c("BRAAK12_SUVR", "BRAAK34_SUVR", "BRAAK56_SUVR")


tau1451.full <- t(sapply(target.rois, function(roi){

  #only adjust for APOE4-status
  f0 <- paste( "eval(", roi,"/", norm, ")","~ eval(APOE4>0)*1 +", confound, "+  eval(AGE+Years.bl)")
  #add scaled score
  f1 <- paste(f0, "+ scale(", score.use, ")")
  #add APOE4 count (to test effect of APOE4 only)
  f2 <- paste(f0, "+ APOE4")

  #run models
  m0 <- lm(as.formula(f0), data=taudata.wn)
  m1 <- lm(as.formula(f1), data=taudata.wn)
  m2 <- lm(as.formula(f2), data=taudata.wn)
  base_score   <- anova(m0, m1)[2,6]
  simple_score <- anova(m0, m2)[2,6]

  #like f0 but using adjustment for the APOE locus
  g0 <- paste( "eval(", roi,"/", norm, ")","~ APOE4 + rs7412 +", confound, "+ eval(AGE+Years.bl)")
  g1 <- paste(g0, "+ scale(", score.use, ")")
  n0 <- lm(as.formula(g0), data=taudata.wn)
  n1 <- lm(as.formula(g1), data=taudata.wn)

  full_score <- anova(n0, n1)[2,6]
  return(c(base_score, full_score, simple_score))
}))

colnames(tau1451.full) <- c(paste(score.use, " (", c("Base","APOE"),")",sep=""),"APOEe44")
print(tau1451.full)

### End Tau imaging ###

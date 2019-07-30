source("./initialize.R")

require(lme4)



##settings for the linear mixed effects model
#Only random intercept
mc <- "(1 | RID)"

#Random intercept and slope
#mc <- "(1 + Years.bl | RID)"

### begin BSI ###


#use foxlab data
fuscans <- subset(foxlabbsi, LONIUID_BASE != "" & MRSEQUENCE != "Acc")
rids <- unique(fuscans$RID)
fuscan.idx <- unlist(sapply(rids, function(rid){
  tmp <- subset(fuscans, RID==rid)
  if (length(table(tmp$MRFIELD)) > 1){
    #if there are 1.5T and 3T scans
    tmp <- subset(tmp, MRFIELD==3.0)
  }
  return(unlist(tmp[,"LONIUID"]))
}))

fuscan.idx2 <- setdiff(fuscan.idx, c("I729643","I384779"))

fuscans2 <- subset(fuscans, is.element(LONIUID, fuscan.idx2))
IMAGEUID <- gsub("I","",paste(fuscans2$LONIUID))
fuscans2 <- data.frame(fuscans2, IMAGEUID)

#manual fixes in the data
fuscans2[fuscans2$LONIUID == "I598363","VISCODE"] <- "m48"

##################################################################################

##NOTE: latest version of ADNIMERGE broke this: VISCODE not available for all...
#bsidata <- merge(adnimerge, fuscans2, by=c("RID","VISCODE"), all=F)
#try IMAGEUID instead, but getting incorrect sample sizes wrt earlier version...
bsidata <- merge(adnimerge, fuscans2, by=c("RID","IMAGEUID"), all=F)

##################################################################################

#restrict to white non-hisp/latino subjects
bsidata.wn <- subset(bsidata, subset=PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino")

#merge with polygenic data
bsidata2.wn <- merge(bsidata.wn, apoe.info2, by="RID", all.y=F)

bsidata2.wn$RID <- as.factor(bsidata2.wn$RID)

#select the target BSI's to work with
targetbsi <- colnames(bsidata2.wn)[grep("BSI",colnames(bsidata2.wn))]

bsi.full <- t(sapply(targetbsi, function(tar){
  #base model including correction for APOE4 carrier status, and all interactions with time, and also scaled score
  f0 <- paste( tar, "~ scale(", score.use,") + ( AGE +", confound, "+ eval(APOE4>0)*1 ) * Years.bl")
  #adding score-by-time interaction
  f1 <- paste(f0, " + scale(", score.use,"):Years.bl",sep="")
  m0 <- lmer(as.formula(paste(f0, "+", mc)), data=bsidata2.wn)
  m1 <- lmer(as.formula(paste(f1, "+", mc)), data=bsidata2.wn)
  base_score <- anova(m0, m1)[2,8]

  #like f0, but replacing the score with APOE4 count (simple score)
  g0 <- paste( tar, "~ APOE4 + ( AGE +", confound, "+ eval(APOE4>0)*1 ) * Years.bl")
  g1 <- paste(g0, " + APOE4:Years.bl",sep="")
  n0 <- lmer(as.formula(paste(g0, "+", mc)), data=bsidata2.wn)
  n1 <- lmer(as.formula(paste(g1, "+", mc)), data=bsidata2.wn)
  simple_score <- anova(n0, n1)[2,8]

  #like f0, but full adjustment for APOE4 locus
  h0 <- paste( tar, "~ scale(", score.use,") + ( AGE +", confound, "+ APOE4 + rs7412 ) * Years.bl")
  h1 <- paste(h0, " + scale(", score.use,"):Years.bl",sep="")
  o0 <- lmer(as.formula(paste(h0, "+", mc)), data=bsidata2.wn)
  o1 <- lmer(as.formula(paste(h1, "+", mc)), data=bsidata2.wn)
  full_score <- anova(o0, o1)[2,8]

  return(c(base_score, full_score, simple_score))
}))

colnames(bsi.full) <- c(paste(score.use, " (",c("Base","APOE"),")",sep=""), "APOEe44")

print(bsi.full)

### END BSI ###

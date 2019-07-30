source("./initialize.R")
require(lme4)

#Only random intercept
mc <- "(1 | RID)"
#Random intercept and slope
#mc <- "(1 + Years.bl | RID)"

#set target variables
target.vars <- c("CDRSB","ADNI_MEM","ADNI_EF")

### begin cognitive decline analysis ###

#require amyloid pet, non-dementia diagnosis and MRI
cogbase <- subset(adnimerge, subset=!is.na(AV45.bl) & DX.bl!= "Dementia" & !is.na(ICV.bl))

#add composite scores for MEM and EF
cogbase2 <- merge(cogbase, uwnpsychsum, by=c("RID","VISCODE"),all.x=T)
dummy <- subset(uwnpsychsum, subset=VISCODE=="bl")[,c("RID","ADNI_MEM","ADNI_EF")]
colnames(dummy)[2:3] <- c("ADNI_MEM.bl","ADNI_EF.bl")
cogbase3 <- merge(cogbase2, dummy, by="RID")

#restict to white non-hisp/latino subjects
cogbase.wn <- subset(cogbase3, subset=PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino")
cogbase.wnnd <- subset(cogbase.wn, subset=DX.bl!="AD")

#add av45 data
cogdata <- merge(cogbase.wnnd, subset(ucberkeleyav45, VISCODE=="bl"), by="RID", all.y=F)

#add polygenic info
cogdata2 <- merge(cogdata, apoe.info2, by="RID", all=F)


cog.nd <- sapply(target.vars, function(tar){
  tbl <- paste(tar, ".bl", sep="")

  #base model APOE4 status only and only interactions with time
  e0 <- paste( "eval(", tar, "-", tbl, ")", "~ ( scale(AGE) +", confound, " + scale(eval(Entorhinal.bl)) + scale(eval(FRONTAL_SUVR/WHOLECEREBELLUM_SUVR)) + eval(APOE4>0)*1 ):Years.bl")
  e1 <- paste(e0, " + scale(", score.use, "):Years.bl", sep="")

  m0 <- lmer(as.formula(paste(e0, "+", mc)), data=cogdata2)
  m1 <- lmer(as.formula(paste(e1, "+", mc)), data=cogdata2)
  base_score_1 <- anova(m0, m1)[2,8]

  #base model APOE4 status only and interactions with time and variables alone
  f0 <- paste( "eval(", tar, "-", tbl, ")", "~ ", score.use," + ( AGE +",confound,"+ scale(eval(Entorhinal.bl)) + scale(eval(FRONTAL_SUVR/WHOLECEREBELLUM_SUVR)) + eval(APOE4>0)*1 ) * Years.bl")
  f1 <- paste(f0, " + ", score.use,":Years.bl",sep="")

  m0 <- lmer(as.formula(paste(f0, "+", mc)), data=cogdata2)
  m1 <- lmer(as.formula(paste(f1, "+", mc)), data=cogdata2)
  base_score_2 <- anova(m0, m1)[2,8]

  #use the simple score instead of polygenic score
  s0 <- paste( "eval(", tar, "-", tbl, ")", "~ APOE4 + ( AGE +", confound," + eval(Entorhinal.bl) + eval(FRONTAL_SUVR/WHOLECEREBELLUM_SUVR) + eval(APOE4>0)*1 ) * Years.bl")
  s1 <- paste(s0, " + APOE4:Years.bl")
  m0 <- lmer(as.formula(paste(s0, "+", mc)), data=cogdata2)
  m1 <- lmer(as.formula(paste(s1, "+", mc)), data=cogdata2)
  simple_score <- anova(m0, m1)[2,8]

  #use the adjustment for full APOE4 locus
  p0 <- paste( "eval(", tar, "-", tbl, ")", "~ ", score.use," + ( AGE +", confound," + eval(Entorhinal.bl) + eval(FRONTAL_SUVR/WHOLECEREBELLUM_SUVR) + APOE4 + rs7412) * Years.bl")
  p1 <- paste(p0, " + ", score.use,":Years.bl", sep="")
  m0 <- lmer(as.formula(paste(p0, "+", mc)), data=cogdata2)
  m1 <- lmer(as.formula(paste(p1, "+", mc)), data=cogdata2)
  full_score <- anova(m0, m1)[2,8]

  return(c(base_score_1, base_score_2, full_score, simple_score))
})

rownames(cog.nd) <- c( paste(score.use, " (", c("Base1","Base2","APOE"),")",sep=""),"APOEe44")
print(cog.nd)

### end cognitive decline analysis ###

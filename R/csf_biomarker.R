source("./initialize.R")

### CSF biomarker ###

#TAU and PTAU are more normal distributed after log logtransform
#ABETA still a bit biomodal after log transform though
#However, p-values were largely unaffected by this.
do.logtransform <- T


#define target biomarker
target.csf <- c("TAU.bl","PTAU.bl","ABETA.bl")

#get subjects with complete CSF biomarker at baseline visit
csfbase <- subset(adnimerge, subset= !is.na(TAU.bl) & !is.na(ABETA.bl) & !is.na(PTAU.bl) & Years.bl==0)

#merge with polygenic info
csfbase2 <- merge(csfbase, apoe.info2, by="RID", all=F)

#restrict to self-reported white non-hisp/latino subjects
csfdata.wn <- subset(csfbase2, subset=PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino")

#replace values for ABETA
csfdata.wn$ABETA.bl <- gsub(">1700","1700",csfdata.wn$ABETA.bl)
#csfdata.wn$ABETA.bl <- gsub(">1700",NA,csfdata.wn$ABETA.bl)
#csfdata.wn <- subset(csfdata.wn, subset=!is.na(ABETA.bl))


##for the machine learning bit:
n_fold <- 10
set.seed(31415)
cvset.wn  <- getCVfold(csfdata.wn, n_fold)

pred.mod <- c()

#loop through target biomaker
csf.full <- t(sapply(target.csf, function(csf){

  #remove smaller than and larger than entries
  bmark <- as.double(gsub("<","",gsub(">","",csfdata.wn[,csf])))
  if (do.logtransform)
    bmark <- log(bmark)

  #only adjust for APOE4 carrier status
  f0 <- paste( "bmark ~ eval(APOE4>0)*1 + AGE +", confound)
  #add scaled polygenic score to f0
  f1 <- paste(f0, "+ scale(", score.use,")")
  #add APOE4 count to f0 (test effect of e4/e4 subjects)
  f2 <- paste(f0, "+ APOE4")

  m0 <- lm(as.formula(f0), data=csfdata.wn)
  m1 <- lm(as.formula(f1), data=csfdata.wn)
  m2 <- lm(as.formula(f2), data=csfdata.wn)
  base_score <- anova(m0, m1)[2,6]
  simple_score <- anova(m0, m2)[2,6]

  #stratify by diagnosis
  base_meta <- metaAnalyze(data.frame(csfdata.wn,bmark), f1, c("CN","MCI","Dementia"), paste("scale(", score.use,")",sep=""))
  simple_meta <- metaAnalyze(data.frame(csfdata.wn,bmark), f2, c("CN","MCI","Dementia"), "APOE4")

  #like f0 but with full adjustment for APOE locus
  g0 <- paste( "bmark ~ APOE4 + rs7412 + AGE +", confound)
  g1 <- paste(g0, "+ scale(", score.use, ")")


  n0 <- lm(as.formula(g0), data=csfdata.wn)
  n1 <- lm(as.formula(g1), data=csfdata.wn)

  full_score <- anova(n0, n1)[2,6]
  full_meta <- metaAnalyze(data.frame(csfdata.wn,bmark), g1, c("CN","MCI","Dementia"), paste("scale(", score.use,")",sep=""))


  #predictive regression analysis
  my.forms <- c(f0, f1, f2, g0, g1)
  names(my.forms) <- c("f0", "f1", "f2", "g0", "g1")

  dummy.data <- data.frame(csfdata.wn, bmark )

  #do for each setting
  mycv <- lapply(my.forms, function(x){
    ttt <- formula2cv(paste(x,"+ DX"), dummy.data, "bmark", cvset.wn, my.fam=gaussian(), measure="cor")
  })
  perf.xxx <- cvsummary(mycv)
  colnames(perf.xxx) <- paste(csf, c("mean","sd","p"))
  pred.mod <<- cbind(pred.mod, perf.xxx)

  return(c(base_score, full_score, simple_score, base_meta[4,4], full_meta[4,4], simple_meta[4,4]))
 }))

colnames(csf.full) <- paste(c(paste(score.use, " (",c("Base","APOE"), ")", sep=""), "APOEe44"),rep(c(""," meta"),each=3),sep="")

print(csf.full)

rownames(pred.mod) <- c("Base", paste(score.use, "(Base)"), "APOEe44", "(APOE)", paste(score.use, "(APOE)"))
print(pred.mod)

### end CSF biomarker ###

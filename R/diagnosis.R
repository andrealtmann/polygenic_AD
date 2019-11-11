source("./initialize.R")

### DIAGNOSIS ###

#obtain lastest diagnosis
rids <- unique(adnimerge$RID)

n_subjects <- length(rids)
n_update   <- floor(n_subjects/10)

n_fold <- 10
set.seed(31415)

if (!exists("ROC_curves"))
  ROC_curves <- list()

if (!exists("admerge.last.idx")){
  cnt <- 0
  message("extracting diagnostic information")
  admerge.last.idx <- t(sapply(rids, function(rid){
    if (cnt %% n_update == 0)
      message(".", appendLF=F)
    cnt <<- cnt + 1
    tmp <- subset(adnimerge, RID==rid & !is.na(DX))
    idx <- order(tmp$Years.bl, decreasing=T)
    last <- idx[1]
    tmp.dx <- tmp[last,"DX"]

    #get first age or last age
    if (is.na(tmp.dx)){
      return( c(rid,"bl",NA,NA))
    }
    if (tmp.dx == "CN"){
      lage <- tmp[last,"AGE"] + tmp[last, "Years.bl"]
      fage <- tmp[last,"AGE"]
    }
    if (tmp.dx == "MCI"){
      lage <- tmp[last,"AGE"] + tmp[last, "Years.bl"]
      mix  <- max(which(tmp[idx,"DX"] == "MCI"))
  	  fage <- tmp[last,"AGE"] + tmp[mix, "Years.bl"]
    }
    if (tmp.dx == "Dementia"){
      lage <- tmp[last,"AGE"] + tmp[last, "Years.bl"]
  	  mix  <- max(which(tmp[idx,"DX"] == "Dementia"))
      fage <- tmp[last,"AGE"] + tmp[mix, "Years.bl"]
    }
    return(c(tmp[last,c("RID","VISCODE")],fage, lage))
  }))
  message(" done!")
  admerge.last.idx <- data.frame(admerge.last.idx)
  colnames(admerge.last.idx)[3:4] <- c("firstAGE","lastAGE")
}

admerge.last <- merge(adnimerge, admerge.last.idx, all=F)
admerge.last <- subset(admerge.last, !is.na(DX))

#restrict analysis to self-reported white non hisp/latino
admerge.last.wn <- subset(admerge.last, subset=PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino")
dxdata <- merge(admerge.last.wn, apoe.info2, all=F)

#f0 is the baseline mode with only adjustment for APOE4 status
f0 <- paste("DX ~ eval(APOE4>0) + testAGE +", confound)

#in f1 we add the scaled polygenic score to model f0
f1 <- paste(f0, " + scale(", score.use,")")

#in f2 we add APOE4 burden to model f0
f2 <- paste(f0, " + ", "APOE4")

##correct modeling of APOE locus
#like f0, but fully accounting for APOE4 locus
g0 <- paste("DX ~ APOE4 + rs7412 + testAGE + ", confound)

#like g0, but adding scaled polygenic score
g1 <- paste(g0, " + scale(", score.use,")")

dxcn  <- dxdata$DX=="CN"
dxmci <- dxdata$DX=="MCI"
dxdem <- dxdata$DX=="Dementia"

my.forms <- c(f0, f1, f2, g0, g1)
names(my.forms) <- c("f0","f1","f2","g0","g1")

testAGE <- unlist(dxdata$firstAGE)
testAGE[dxcn] <- unlist(dxdata[dxcn,"lastAGE"])

#HC vs AD
m0A <- glm(as.formula(f0), family=binomial, data=dxdata, subset=DX!="MCI")
m1A <- glm(as.formula(f1), family=binomial, data=dxdata, subset=DX!="MCI")
m2A <- glm(as.formula(f2), family=binomial, data=dxdata, subset=DX!="MCI")

#prepare df for predictive modeling
dummy <- data.frame(dxdata, testAGE, dxdem)
dxdata.hcad <- subset(dummy, DX!="MCI")
cvset.hcad  <- getCVfold(dxdata.hcad, n_fold)

#do for each setting
mycv <- lapply(my.forms, function(x){
  ttt <- formula2cv(x, dxdata.hcad, "dxdem", cvset.hcad, my.fam=binomial())
})
perf.hcad <- cvsummary(mycv)

#get the ROC curves
for (rmod in names(my.forms)){
  ROC_curves[[paste("HCAD",score.use, rmod, sep="_")]] <- getROC(mycv, rmod)
}

#HC vs MCI
m0B <- glm(as.formula(f0), family=binomial, data=dxdata, subset=DX!="Dementia")
m1B <- glm(as.formula(f1), family=binomial, data=dxdata, subset=DX!="Dementia")
m2B <- glm(as.formula(f2), family=binomial, data=dxdata, subset=DX!="Dementia")

#prepare df for predictive modeling
dummy <- data.frame(dxdata, testAGE, dxmci)
dxdata.hcmci <- subset(dummy, DX!="Dementia")
cvset.hcmci  <- getCVfold(dxdata.hcmci, n_fold)

#do for each setting
mycv <- lapply(my.forms, function(x){
  ttt <- formula2cv(x, dxdata.hcmci, "dxmci", cvset.hcmci, my.fam=binomial())
})
perf.hcmci <- cvsummary(mycv)

#get the ROC curves
for (rmod in names(my.forms)){
  ROC_curves[[paste("HCMCI",score.use, rmod, sep="_")]] <- getROC(mycv, rmod)
}


#HC vs MCI and AD
m0D <- glm(as.formula( gsub("DX","dxcn",f0) ), family=binomial, data=dxdata)
m1D <- glm(as.formula( gsub("DX","dxcn",f1) ), family=binomial, data=dxdata)
m2D <- glm(as.formula( gsub("DX","dxcn",f2) ), family=binomial, data=dxdata)

#prepare df for predictive modeling
dummy        <- data.frame(dxdata, testAGE, dxcn)
dxdata.hcall <- subset(dummy)
cvset.hcall  <- getCVfold(dxdata.hcall, n_fold)

#do for each setting
mycv <- lapply(my.forms, function(x){
  ttt <- formula2cv(x, dxdata.hcall, "dxcn", cvset.hcall, my.fam=binomial())
})
perf.hcall <- cvsummary(mycv)

#get the ROC curves
for (rmod in names(my.forms)){
  ROC_curves[[paste("HCALL",score.use, rmod, sep="_")]] <- getROC(mycv, rmod)
}


#HC vs AD
p1 <- anova(m0A, m1A, test="Chisq")[2,5]
p2 <- anova(m0A, m2A, test="Chisq")[2,5]

#HC vs MCI
q1 <- anova(m0B, m1B, test="Chisq")[2,5]
q2 <- anova(m0B, m2B, test="Chisq")[2,5]

#HC vs MCI and AD
s1 <- anova(m0D, m1D, test="Chisq")[2,5]
s2 <- anova(m0D, m2D, test="Chisq")[2,5]


n0A <- glm(as.formula(g0), family=binomial, data=dxdata, subset=DX!="MCI")
n1A <- glm(as.formula(g1), family=binomial, data=dxdata, subset=DX!="MCI")


n0B <- glm(as.formula(g0), family=binomial, data=dxdata, subset=DX!="Dementia")
n1B <- glm(as.formula(g1), family=binomial, data=dxdata, subset=DX!="Dementia")

n0D <- glm(as.formula( gsub("DX","dxcn",g0)), family=binomial, data=dxdata)
n1D <- glm(as.formula( gsub("DX","dxcn",g1)), family=binomial, data=dxdata)

p3 <- anova(n0A, n1A, test="Chisq")[2,5]
q3 <- anova(n0B, n1B, test="Chisq")[2,5]
s3 <- anova(n0D, n1D, test="Chisq")[2,5]

##MCI vs AD
testAGE <- unlist(dxdata$firstAGE)
testAGE[dxmci] <- unlist(dxdata[dxmci,"lastAGE"])
m0C <- glm(as.formula(f0), family=binomial, data=dxdata, subset=DX!="CN")
m1C <- glm(as.formula(f1), family=binomial, data=dxdata, subset=DX!="CN")
m2C <- glm(as.formula(f2), family=binomial, data=dxdata, subset=DX!="CN")

#prepare df for predictive modeling
dummy        <- data.frame(dxdata, testAGE, dxdem)
dxdata.mciad <- subset(dummy)
cvset.mciad  <- getCVfold(dxdata.mciad, n_fold)

#do for each setting
mycv <- lapply(my.forms, function(x){
  ttt <- formula2cv(x, dxdata.mciad, "dxdem", cvset.mciad, my.fam=binomial())
})
perf.mciad <- cvsummary(mycv)

#get the ROC curves
for (rmod in names(my.forms)){
  ROC_curves[[paste("MCIAD",score.use, rmod, sep="_")]] <- getROC(mycv, rmod)
}


n0C <- glm(as.formula(g0), family=binomial, data=dxdata, subset=DX!="CN")
n1C <- glm(as.formula(g1), family=binomial, data=dxdata, subset=DX!="CN")

r1 <- anova(m0C, m1C, test="Chisq")[2,5]
r2 <- anova(m0C, m2C, test="Chisq")[2,5]
r3 <- anova(n0C, n1C, test="Chisq")[2,5]


#put all in a matrix
dx.results <- matrix(NA, nrow=3, ncol=4)
rownames(dx.results) <- c(paste(score.use, " (", c("Base","APOE"), ")",sep=""), "APOEe44")
colnames(dx.results) <- c("HC|AD","HC|MCI","MCI|AD","HC|MCI_AD")

dx.results[1,] <- c(p1,q1,r1,s1)
dx.results[2,] <- c(p3,q3,r3,s3)
dx.results[3,] <- c(p2,q2,r2,s2)

print(dx.results)

pred.perf <- cbind(perf.hcad, perf.hcmci, perf.mciad, perf.hcall)
colnames(pred.perf) <- rep(c("mean","sd","pv"),4)
rownames(pred.perf) <- c("Base", paste(score.use, "(Base)"), "APOEe44", "(APOE)", paste(score.use, "(APOE)"))

print(pred.perf)

### END DIAGNOSES###

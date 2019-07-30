require(survival)

source("./initialize.R")

##settings



### begin conversion analysis ###
rids <- unique(adnimerge$RID)
target.dx <- "Dementia"

n_subjects <- length(rids)
n_update   <- floor(n_subjects/10)

#extract time to dementia diagnosis
if (!exists("timetoad")){
  cnt <- 0
  message("extracting time to dementia")
  timetoad <- t(sapply(rids, function(x){
    if (cnt %% n_update == 0)
      message(".", appendLF=F)
    cnt <<- cnt + 1

    tmp <- subset(adnimerge, RID==x)
    N <- nrow(tmp)
    idx <- order(tmp$Years.bl)
    bldx <- paste(tmp[idx[1],"DX"])
    blage <- tmp[idx[1],"AGE"]

    conv <- which(tmp[idx,"DX"] == target.dx)
    event <- 0
    if (length(conv)== 0){
      fu <- max(which(tmp[idx,"DX"] != target.dx))
    } else {
      event <- 1
      fu <- min(which(tmp[idx,"DX"] == target.dx))
    }
    endage <- blage + tmp[idx[fu],"Years.bl"]
    enddx  <- paste(tmp[idx[fu],"DX"])

    return(c(x, blage, endage, event))
  }))
  message(" done!")

  colnames(timetoad) <- c("RID","AGE","AGE2","DEM")
  timetoad <- data.frame(timetoad)
}

#baseline info
ambl <- subset(adnimerge, subset=Years.bl==0)
#add conversion info
surv1 <- merge(ambl, timetoad, by=c("RID","AGE"),all=F)
#add polygenic info
surv2 <- merge(surv1, apoe.info2, by="RID", all=F)

#restrict to white non-hisp/latino
surv2.wn <- subset(surv2, subset=PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino")
#exclude people with AD
surv2.wnnd <- subset(surv2.wn, subset=DX.bl != "AD" & AGE2-AGE>0)

#typle of model:
#left trunkated, right censored
sv0 <- Surv(surv2.wnnd$AGE, surv2.wnnd$AGE2, surv2.wnnd$DEM)

#classic (right censored, but work on time since inclusion)
#sv0 <- Surv(surv2.wnnd$AGE2 - surv2.wnnd$AGE, surv2.wnnd$DEM)

#adjusted only for APOE4 status
f0 <- paste("sv0 ~ eval(APOE4>0)*1 + eval(AGE-mean(AGE,na.rm=T)) + strata(DX) + ", confound)
#add scaled score
f1 <- paste(f0, "+ scale(", score.use, ")")
#add APOE4 count
f2 <- paste(f0, "+ APOE4")

#run the models
m0 <- coxph(as.formula(f0), data=surv2.wnnd)
m1 <- coxph(as.formula(f1), data=surv2.wnnd)
m2 <- coxph(as.formula(f2), data=surv2.wnnd)
base_score   <- anova(m0,m1)
simple_score <- anova(m0,m2)

#adjusted for full APOE4 locus
g0 <- paste("sv0 ~  APOE4 + rs7412 + eval(AGE-mean(AGE,na.rm=T)) + strata(DX) +", confound)
g1 <- paste(g0, "+ scale(", score.use, ")")
n0 <- coxph(as.formula(g0), data=surv2.wnnd)
n1 <- coxph(as.formula(g1), data=surv2.wnnd)
full_score <- anova(n0,n1)

clin_conv <- c(base_score[2,4], full_score[2,4], simple_score[2,4])
names(clin_conv) <- c(paste(score.use, " (", c("Base","APOE"),")",sep=""),"APOEe44")

print(clin_conv)

### end conversion analysis ###

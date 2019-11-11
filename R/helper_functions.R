
#ds is the dataset
#form is the formula (as string)
#grps are the group names
#targetvar is the variable name to be extracted from results

metaAnalyze <- function(ds, form, grps, targetvar){


  #message(form)

  #run models for each group separately
  mods <- sapply(grps,function(gg){
    tmp <- lm(as.formula(form), data=ds, subset=DX==gg)
    res <- summary(tmp)$coeff[targetvar,]
    return(res)
  })

  #combine effect sizes using inv variance weighted
  wi <- 1 / mods["Std. Error",]^2
  SE <- sqrt(1/sum(wi))
  BE <- sum(mods["Estimate",] * wi)/sum(wi)
  Z <- BE/SE
  P <- pnorm(abs(Z), lower.tail=F)*2
  META <- c(BE, SE, Z, P)

  return(cbind(mods, META))
}

#add correlation to pairs plot
panel.cor <- function(x, y, ...){
  par(usr = c(0,1,0,1))
  txt <- as.character(format(cor(x,y)^2, digit=2))
  text(0.5, 0.5, txt, cex=6 * abs(cor(x,y)))
}


### cross-validation and predictive modeling

## create a CV set ###
getCVfold <- function(x, nfold){
    avail <- 1:nrow(x)
    nsize <- ceiling(nrow(x)/nfold)
    cv.idx <- rep(nfold, nrow(x))
    for (i in 1:(nfold-1)){
        tmp <- sample(avail, nsize)
        cv.idx[tmp] <- i
        avail <- setdiff(avail, tmp)
    }
    return(cv.idx)
}

### run CV with a glm ###
### y is the label/outcome and x are the features
### set fam to binomial for log.reg, gaussian for lin.reg

require(ROCR)
runLinearCV <- doubleCV <- function(y, x, nfold, measure="auc", fam=gaussian(), myfold=c()){

    ##if you don't provide a CV fold, one will be created
    if (length(myfold) == 0){
      cv.idx <- getCVfold(x, nfold)
    } else {
        cv.idx <- myfold
        nfold  <- max(myfold)
    }

    my2CV.predict <- list()
    my2CV.labels  <- list()

    #run CV
    for(i in 1:nfold){
        message("outer fold: ",i)
        test.y <- y[cv.idx==i]
        test.x <- x[cv.idx==i,]

        train.y <- y[cv.idx!=i]
        train.x <- x[cv.idx!=i,]

        mymod <- glm( train.y ~ ., family=fam, data=as.data.frame(train.x))
        #print(mymod)

        my2CV.predict[[i]] <- predict(mymod, newdata=as.data.frame(test.x), type="response")
        my2CV.labels[[i]] <- test.y
    }

    res <- rep(0, nfold)

    if (measure=="auc"){
      res <- unlist(performance(prediction(my2CV.predict, my2CV.labels), "auc")@y.values)
    } else {
      res <- t(sapply(1:nfold, function(i){
        cr   <- cor(my2CV.predict[[i]], my2CV.labels[[i]])
        rmse <- sqrt(mean((my2CV.predict[[i]] - my2CV.labels[[i]])^2))
        #return(c(cr, rmse))
        return(cr)
      }))
    }
    metrics <-  res
    toreturn <- list()
    toreturn[["metrics"]] <- metrics
    toreturn[["prediction"]] <- my2CV.predict
    toreturn[["labels"]]     <- my2CV.labels
    return(toreturn)
}



formula2cv <- function(my.formula, dat, lab, cv.set, my.fam, measure="auc", ...){

  xmat   <- model.matrix(as.formula(my.formula), dat)
  xmat   <- xmat[,2:ncol(xmat)]

  if (measure == "surv"){
    mycv.form <- runCoxPhCV(dat[,lab], xmat, myfold=cv.set)
    return(mycv.form)
  }
  mycv.form  <- runLinearCV( dat[,lab], xmat, fam=my.fam, myfold=cv.set, measure=measure)


}

cvres2entry <- function(res1, res2=NA){

  mmm <- mean(res1$metrics)
  sss <- sd(res1$metrics)
  pv <- NA
  if ( sum(is.na(res2)) < 1 ){
    pv <- wilcox.test( res2$metrics, res1$metrics, paired=T, alternative="less")$p.value
  }
  return(c(round(mmm,3), round(sss,3), pv))
}

cvsummary <- function(mycvs){

  tmp.perf <- c()
  #base model
  tmp.perf <- rbind(tmp.perf, cvres2entry(mycvs[["f0"]]))
  tmp.perf <- rbind(tmp.perf, cvres2entry(mycvs[["f1"]], mycvs[["f0"]]))
  tmp.perf <- rbind(tmp.perf, cvres2entry(mycvs[["f2"]], mycvs[["f0"]]))
  tmp.perf <- rbind(tmp.perf, cvres2entry(mycvs[["g0"]], mycvs[["f0"]]))
  tmp.perf <- rbind(tmp.perf, cvres2entry(mycvs[["g1"]], mycvs[["g0"]]))

  return(tmp.perf)

}

getROC <- function(mycvs, model){

  lab  <- mycvs[[model]][["labels"]]
  pred <- mycvs[[model]][["prediction"]]
  return(performance(prediction(pred, lab), 'tpr','fpr'))

}

require(survival)
require(survcomp)
runCoxPhCV <- doubleCV <- function(y, x, nfold, myfold=c()){

    ##if you don't provide a CV fold, one will be created
    if (length(myfold) == 0){
      cv.idx <- getCVfold(x, nfold)
    } else {
        cv.idx <- myfold
        nfold  <- max(myfold)
    }

    my2CV.predict <- list()
    my2CV.conc <- list()
    my2CV.labels  <- list()

    #run CV
    for(i in 1:nfold){
        message("outer fold: ",i)
        test.y <- y[cv.idx==i]
        test.x <- x[cv.idx==i,]

        train.y <- y[cv.idx!=i]
        train.x <- x[cv.idx!=i,]

        Ttrain <- data.frame(sv=train.y, train.x)
        Ttest  <- data.frame(sv=test.y, test.x)

        mymod <- coxph( sv ~ ., data=Ttrain, x=T)

        my2CV.predict[[i]] <- predict(mymod, newdata=Ttest, type="lp")
        my2CV.labels[[i]] <- test.y

        my2CV.conc[[i]] <- concordance(mymod)


    }

    toreturn <- list()
    toreturn[["metrics"]] <- NA
    toreturn[["prediction"]]  <- my2CV.predict
    toreturn[["labels"]]      <- my2CV.labels
    toreturn[["concordance"]] <- my2CV.conc

    return(toreturn)
}

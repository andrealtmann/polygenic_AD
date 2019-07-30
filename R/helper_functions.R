
#ds is the dataset
#form is the formula (as string)
#grps are the group names
#targetvar is the variable name to be extracted from results

metaAnalyze <- function(ds, form, grps, targetvar){


  message(form)

  mods <- sapply(grps,function(gg){
    tmp <- lm(as.formula(form), data=ds, subset=DX==gg)
    res <- summary(tmp)$coeff[targetvar,]
    return(res)
  })

  #inv variance weighted
  wi <- 1 / mods["Std. Error",]^2
  SE <- sqrt(1/sum(wi))
  BE <- sum(mods["Estimate",] * wi)/sum(wi)
  Z <- BE/SE
  P <- pnorm(abs(Z), lower.tail=F)*2
  META <- c(BE, SE, Z, P)

  return(cbind(mods, META))
}


panel.cor <- function(x, y, ...){
  par(usr = c(0,1,0,1))
  txt <- as.character(format(cor(x,y)^2, digit=2))
  text(0.5, 0.5, txt, cex=6 * abs(cor(x,y)))
}

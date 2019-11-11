##this is a script to just generate all the ROC curves needed
##and then plot a subset of them

mkroc <- T

pr.scores <- c("PHS", "PRS1", "PRS2","PRScs_auto")
prstr <- c("PHS","PRS1","PRS2","PRS-cs")
names(prstr) <- pr.scores

grps <- c("HCAD","HCMCI","MCIAD","HCALL")
grptr <- c("CN|AD", "CN|MCI", "MCI|AD", "CN|MCI + AD")
names(grptr) <- grps

#mycolors <- c(rgb(0,0,0), rgb(178,152,255), rgb(255,204,110), rgb(0,31,121), rgb(198,5,76)) #, rgb(15,66,0))
#mycolors <- c("#000000","#b298ff","#ffcc6e","#0f4200","#c6054c") # "#001f79")
mycolors <- c("#000000","#697cd4","#ba495c","#b0913b","#56ae6c")
names(mycolors) <- c("null", pr.scores)

for (score.use in pr.scores){
  source("diagnosis.R")
}


#grp <- "HCMCI"


### for the reduced dataset
#tiff("/tmp/ROC_reducedds.tiff", height=3, width=9, res=600, units='in')
par(mfrow=c(1,3))
for(grp in setdiff(grps, "HCALL")){
  plot(ROC_curves[[paste(grp, "PHS","g0",sep="_")]], avg="vertical", lwd=1, col=mycolors["null"], main=grptr[grp])
  #plot(ROC_curves[[paste(grp, "PHS","f0",sep="_")]], avg="vertical", lwd=1, lty=2, col=mycolors["null"], add=T)
  for (s in pr.scores){
    plot(ROC_curves[[paste(grp,s,"g1", sep="_")]], avg="vertical", lwd=1, col=mycolors[s], add=T)
  }
  abline(a=0,b=1, col="grey", lwd=1, lty=2)
  legend(0.45, 0.5, c("burden", paste("+",prstr[pr.scores])), lwd=1, lty=1, col=mycolors[c("null",pr.scores)], bty="n")
}

### for the full dataset
#tiff("/tmp/ROC_fullds.tiff", height=3, width=9, res=600, units='in')
par(mfrow=c(1,3))
for(grp in setdiff(grps, "HCALL")){
  plot(ROC_curves[[paste(grp, "PHS","g0",sep="_")]], avg="vertical", lwd=1, col=mycolors["null"], main=grptr[grp])
  plot(ROC_curves[[paste(grp, "PHS","f0",sep="_")]], avg="vertical", lwd=1, lty=2, col=mycolors["null"], add=T)
  s <- "PHS"
  plot(ROC_curves[[paste(grp,s,"g1", sep="_")]], avg="vertical", lwd=1, col=mycolors[s], add=T)
  #plot(ROC_curves[[paste(grp,s,"f1", sep="_")]], avg="vertical", lwd=1, lty=2, col=mycolors[s], add=T)
  abline(a=0,b=1, col="grey", lwd=1, lty=2)
  legend(0.55, 0.35, c("status","burden", paste("+",s)), lwd=1, lty=c(2,1,1), col=mycolors[c("null","null",s)], bty="n")
}

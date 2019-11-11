source("./initialize.R")


##for scatter plot
dummy <- merge(subset(adnimerge, Years.bl==0), apoe.info2, by="RID", all=F)
dummy2 <- subset(dummy, PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino")
pairs(dummy2[,c("PHS","PRS1","PRS2","PRScs_auto")], lower.panel=panel.cor, labels=c("PHS","PRS1","PRS2","PRS-cs"), pch=20, cex=0.5, col=dummy2[,"rs429358"]+1)



##for demographics
data_demog <- merge(subset(adnimerge, Years.bl==0 & PTRACCAT=="White" & PTETHCAT=="Not Hisp/Latino"), apoe.info2, by="RID")

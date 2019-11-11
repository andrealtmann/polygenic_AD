
mk_vtk_table <- function(pvals, lbase, maxint=10){
  fs.map <- read.csv("FS_mapping.csv")
  fs.map.left <- subset(fs.map, Hemisphere=="Left")

  pv <- -log10(p.adjust(pvals[lbase],'fdr'))
  lbase2 <- gsub("_SUVR","", lbase)

  fs_abbrev2 <- toupper(gsub("L_","",fs.map.left$FS_abbrev))
  vtkcode <- fs.map.left$code - 1000
  fs.map.left <- data.frame(fs.map.left, fs_abbrev2, vtkcode)

  fs.map.out <- data.frame(fs.map.left, PV=pv[fs.map.left$fs_abbrev2])

  omap <- matrix(NA,nrow=40, ncol=4)
  omap[,1] <- 1:40
  rownames(omap) <- 1:40
  
  cp <- colorRampPalette(c("red","yellow"))(10 * maxint)
  dummy <- floor(fs.map.out$PV*10)
  dummy[dummy == 0] <- 1
  dummy[dummy >10 * maxint] <- 10 * maxint
  hex.col <- cp[dummy]
  rgb.col <- data.frame(t(col2rgb(hex.col))) 
  rownames(rgb.col) <- fs.map.out$vtkcode

  omap[rownames(rgb.col),2] <- rgb.col[,1]
  omap[rownames(rgb.col),3] <- rgb.col[,2]
  omap[rownames(rgb.col),4] <- rgb.col[,3]
  colnames(omap) <- c("label","red","green","blue")
  omap["40",2:4] <- rep(192,3)
  omap[is.na(omap)] <- 0

  #gray out non sign results
  idx <- paste(fs.map.out[fs.map.out$PV < -log10(0.05),"vtkcode"])
  omap[idx,2] <- 192
  omap[idx,3] <- 192
  omap[idx,4] <- 192

  return(omap)
}
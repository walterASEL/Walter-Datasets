library(ks)
library(rgdal)
library(maptools)
library(PBSmapping)
library(rgeos)

# get input file
indata <- read.csv("White10pelicans.csv")
innames <- unique(indata$ID)
innames <- innames[1:5]
outnames <- innames

# set up output table
output <- as.data.frame(matrix(0,nrow=length(innames),ncol=9))
colnames(output) <- c("ID","noFixes","h11","h12","h21","h22",
                      "iso50areaKm","iso95areaKm","iso99areaKm")

# set up levels for home range
levels <- c(50,95,99)
# set up output directory for shp files
dirout <- "output"
# begin loop to calculate home ranges
for (i in 1:length(innames)){
  data <- indata[which(indata$ID==innames[i]),]
  if(dim(data)[1] != 0){
    # export the point data into a shp file
    data.xy = data[c("Longitude", "Latitude")]
    coordinates(data.xy) <- ~Longitude+Latitude
    sppt <- SpatialPointsDataFrame(coordinates(data.xy),data)
    proj4string(sppt) <- CRS("+proj=longlat +ellps=WGS84")
    writePointsShape(sppt,fn=paste(outnames[i],sep="/"),factor2char=TRUE)
    # start populating output table
    output$ID[i] <- as.character(data$ID[1])
    output$noFixes[i] <- dim(data)[1]
    locs <- cbind(data$Longitude,data$Latitude)
    try(HpiOut <- Hpi(locs,pilot="samse",binned=TRUE))
    if(is.null(HpiOut)=="FALSE"){
      output$h11[i] <- HpiOut[1,1]
      output$h12[i] <- HpiOut[1,2]
      output$h21[i] <- HpiOut[2,1]
      output$h22[i] <- HpiOut[2,2]
      fhatOut <- kde(x=locs,H=HpiOut)
    }
    if(is.null(fhatOut)=="FALSE"){
      for (j in 1:length(levels)){
        fhat.contlev <- contourLevels(fhatOut, cont=c(levels[j]))
        fhat.contlines <- contourLines(x=fhatOut$eval.points[[1]],y=fhatOut$eval.points[[2]], z=fhatOut$estimate, level=fhat.contlev)
        # convert contour lines into spatial objects to export as polygon shp file
        sldf <- ContourLines2SLDF(fhat.contlines)
        proj4string(sldf) <- CRS("+proj=longlat +ellps=WGS84")
        ps <- SpatialLines2PolySet(sldf)
        attr(ps,"projection") <- "LL"
        sp <- PolySet2SpatialPolygons(ps)
        dataframe <- as.data.frame(matrix(as.character(1,nrow=1,ncol=1)))
        spdf <- SpatialPolygonsDataFrame(sp,dataframe,match.ID=TRUE)
        # get area and export shp files
        if (j == 1){
        pls <- slot(spdf, "polygons")[[1]]
        gpclibPermit()
        xx <- checkPolygonsHoles(pls)
        a <- sapply(slot(xx, "Polygons"), slot, "area")
        h <- sapply(slot(xx, "Polygons"), slot, "hole")
        output$iso50areaKm[i] <- sum(ifelse(h, -a, a))/1000000
        writeOGR(spdf,dirout,paste(outnames[i],"KUD50",sep=""),"ESRI Shapefile")}
        if (j == 2){
        pls <- slot(spdf, "polygons")[[1]]
        gpclibPermit()
        xx <- checkPolygonsHoles(pls)
        a <- sapply(slot(xx, "Polygons"), slot, "area")
        h <- sapply(slot(xx, "Polygons"), slot, "hole")
        output$iso95areaKm[i] <- sum(ifelse(h, -a, a))/1000000
        writeOGR(spdf,dirout,paste(outnames[i],"KUD95",sep=""),"ESRI Shapefile")}
        if (j == 3){
        pls <- slot(spdf, "polygons")[[1]]
        gpclibPermit()
        xx <- checkPolygonsHoles(pls)
        a <- sapply(slot(xx, "Polygons"), slot, "area")
        h <- sapply(slot(xx, "Polygons"), slot, "hole")
        output$iso99areaKm[i] <- sum(ifelse(h, -a, a))/1000000
        writeOGR(spdf,dirout,paste(outnames[i],"KUD99",sep=""),"ESRI Shapefile")}
      }
    }
  }
  rm(data,data.xy,sppt,locs,HpiOut,fhatOut,fhat.contlev,
     fhat.contlines,sldf,ps,sp,dataframe,spdf)
    
}
# write output
write.csv(output,paste("output.csv",sep=""))



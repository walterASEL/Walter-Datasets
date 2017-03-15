library(adehabitatHR)

#Creates a Spatial Points Data Frame for 2 animals by ID
allcats <-read.csv("pantherjitter.csv", header=T)
#allcats2 <- subset(allcats, allcats$id=="FP048" | allcats$id=="FP131" | allcats$id=="FP135" | allcats$id=="FP143") 
#allcats <- read.csv("C:\\Walter\\WalterSpatialEcologyLab\\SpatialEcologyCourse\\Chapter5\\Overlap\\subsetHRlocs.csv", header=T)

data.xy = allcats[c("X","Y")]

#Creates class Spatial Points for all locations
xysp <- SpatialPoints(data.xy)
proj4string(xysp) <- CRS("+proj=utm +zone=17N +ellps=WGS84")

#Creates a Spatial Data Frame from 
sppt<-data.frame(xysp)

#Creates a spatial data frame of ID
idsp<-data.frame(allcats[1])

#Merges ID data frame with GPS locations data frame
#Data frame is called "idsp" comparable to the "relocs" from puechabon dataset
coordinates(idsp)<-sppt

#First we need to create utilization distributions for each panther
ud <- kernelUD(idsp[,1], h = "href", grid = 200, same4all = TRUE, hlim = c(0.1, 1.5), kern = c("bivnorm"), extent = 0.5)
#ud <- kernelUD(idsp[,1]), same4all=TRUE)

#output of UDs for each panther
image(ud)

#Then we can compute volume of intersection scores
kerneloverlaphr(ud, meth="HR", conditional=TRUE)
kerneloverlaphr(ud, meth="PHR", conditional=TRUE)
kerneloverlaphr(ud, meth="BA", conditional=TRUE)
kerneloverlaphr(ud, meth="UDOI", conditional=TRUE)
kerneloverlaphr(ud, meth="HD", conditional=TRUE)
kerneloverlaphr(ud, meth="VI", conditional=TRUE)


plot(idsp, col="yellow")
uds <- getverticeshr(ud)
plot(uds)

plot(idsp, col=idsp@data[,1], add=T)
ud1 <- getverticeshr(ud[[1]])
plot(ud1, lty=5)
ud2 <- getverticeshr(ud[[2]])
plot(ud2, lwd=5)


#An alternative way without creating UD in a separate step
kerneloverlap(idsp[,1], grid=200, method="HR", percent=95, conditional=TRUE)
kerneloverlap(idsp[,1], grid=200, method="PHR", percent=95, conditional=TRUE)
kerneloverlap(idsp[,1], grid=200, method="BA", percent=95, conditional=TRUE)
kerneloverlap(idsp[,1], grid=200, method="UDOI", percent=95, conditional=TRUE)
kerneloverlap(idsp[,1], grid=200, method="HD", percent=95, conditional=TRUE)
kerneloverlap(idsp[,1], grid=200, method="VI", percent=95, conditional=TRUE)


install.packages("gpclib")

library(adehabitatHR)
library(maptools)

#Loading shapefile is easier than text or CSV that require only x and y and does not seem to import as easily
#Load the LoCoH home range R code
source("NNCH.R")


#Load the LoCoH home range GUI interface
source("locoh_gui.R")

#Opens the GUI
locoh()

panther <- read.csv("pantherjitter2.csv")
str(panther)
panther$CatID <- as.factor(panther$CatID)

#Or explore with one panther with 381 relocations
cat159 <- subset(panther, CatID=="159")
str(cat159)
cat159$CatID <- factor(cat159$CatID)

#Get the relocation data from the source file
data.xy = cat159[c("x","y")]


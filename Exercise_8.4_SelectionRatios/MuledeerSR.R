library(adehabitatHS)

MDsr <- read.csv("MD_winter12.csv",header=T)
#Remove deer that cause errors in plot function later
MDsr <- subset(MDsr,MDsr$animal_id !="647579A")
MDsr$animal_id <- factor(MDsr$animal_id)
used <- subset(MDsr, MDsr$use == 1)
used <- used[c(-1,-3:-6,-8:-15)]
used <- xtabs(~used$animal_id + used$crop, used)
used <- as.data.frame.matrix(used[1:13, 1:5])

rand <- subset(MDsr, MDsr$use == 0)
rand <- rand[c(-1,-3:-6,-8:-16)]
rand <- xtabs(~rand$animal_id + rand$crop, rand)
rand <- as.data.frame.matrix(rand[1:13, 1:5])

# PVT Code for VegRSF #
pvt.W <- widesIII(used,rand,avknown = FALSE, alpha = 0.1)
pvt.W
par(mfrow=c(1,2))
plot(pvt.W)

#Five categories
#1 = Sunflower,summer crops, random crops, grassland
#2 = Winter crops
#3 = Alfalfa
#4 = Forest
#5 = Shrubland

#Now run on distance to roads binned into 10 categories
MDsr <- read.csv("MD_winter12.csv",header=T)
#Delete deer that have limited data
MDsr <- subset(MDsr,MDsr$animal_id !="647582A" & MDsr$animal_id !="647584A")
MDsr$animal_id <- factor(MDsr$animal_id)

#Bin roads into 4 categories instead of 10
MDsr$NewRoad <- MDsr$BinRoad
levels(MDsr$NewRoad)<-list(class1=c("0-200","200-400"), class2=c("400-600","600-800"),class3=c("800-1000","1000-12000","1200-1400"),class4=c("1400-1600","1600-1800","1800-2000"))

used <- subset(MDsr, MDsr$use == 1)
used <- used[c(-1:-6,-8:-15)]
used <- xtabs(~used$animal_id + used$NewRoad, used)
used <- as.data.frame.matrix(used[1:12, 1:4])

rand <- subset(MDsr, MDsr$use == 0)
rand <- rand[c(-1:-6,-8:-15)]
rand <- xtabs(~rand$animal_id + rand$NewRoad, rand)
rand <- as.data.frame.matrix(rand[1:12, 1:4])

pvt.road <- widesIII(use,rand,avknown = FALSE, alpha = 0.1)
pvt.road
par(mfrow=c(1,2))
plot(pvt.road)




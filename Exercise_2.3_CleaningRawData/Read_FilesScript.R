# TO RUN: highlight code, hit CTRL R

# NOTE: Be sure that the only .txt files in your working directory are 
# the flow files before  running this code

# Vector of files names in working directory
files <- list.files(pattern = ".txt")
files
# Total number of files in working directory (for loop below)
n.files <- length(files)
n.files
# Container to hold text files
files.list <- list()
#files.list <- as.factor(files$ELEV)

# (populate the container files.list with climate data sets
files.list <- lapply(files, read.table, header =T, sep="\t") 

# Set up matrix for weather station summary data
m1 <- matrix(NA,ncol=8,nrow=n.files)

# Loop for running through all weather station files
for(i in 1:n.files){
      
	# Assign elevation
        m1[i,1] <- files.list[[i]][1,10]

	#Assign Lat
        m1[i,2] <- files.list[[i]][1,11]

	#Assign Long
        m1[i,3] <- files.list[[i]][1,12]

	#Calculate mean snow depth
        SNWD_mm <- mean(files.list[[i]][,7],na.rm=T)

	#Convert snow depth mean to inches
	SNWD_in <- SNWD_mm/25.4

	#Assign snow depth
	m1[i,4] <- SNWD_in

	#Calculate mean maximum temp
        TMAX_C <- mean(files.list[[i]][,8],na.rm=T)

	#Convert max temp to F
	TMAX_F <- TMAX_C*0.18 + 32
	
	#Assign max temp
	m1[i,5] <- TMAX_F

	#Calculate mean minimum temp
	TMIN_C <- mean(files.list[[i]][,9],na.rm=T)

	#Convert min temp to F
	TMIN_F <- TMIN_C*0.18 + 32

	#Assign min temp
	m1[i,6] <- TMIN_F

	#Reassign GHCN number
	GHCN <- toString(files.list[[i]][1,1])

	#Assign Station Name
	m1[i,7] <- GHCN

	#Reassign Station Name
	SN <- toString(files.list[[i]][1,2])

	#Assign Station Name
	m1[i,8] <- SN
}

colnames(m1) <- c("Elevation","Lat","Long","SNWD","TMAX","TMIN","GHCN","Station")
write.csv(m1,paste(".","\\output.csv",sep=""))

#Removes quotation marks in output table
m1 <-noquote(m1)
m1[1:5,]

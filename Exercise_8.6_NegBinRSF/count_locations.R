#
# USED IN "desert_elk.R"
#
# CREATED BY RYAN NIELSON (WEST, INC.) 2012.
# MODIFIED BY FAWN HORNSBY (WEST, INC.) 2014.


#=====================================================================================

F.count.relocations <- structure(
function # Counts the number of animal locations within a circular sampling unit

# Description: Counts the number of animal locations within a circular sampling unit.
# The input for habitat.units.df should contain a systematic sample of points or a 
# random sample of points to be used as center points of each circular sampling unit.
# The input for habitat.unit.size is the radius of each circular sampling unit.
#
# Returns: A dataframe containing all the center points and covariates provided in 
# the habitat.units.df argument as well as the number of locations found within each
# sampling unit (column is called "n.locations") and the total animal locations (column
# is called "total")

(locations.df, 
### A dataframe containing animal locations.  There should be 3 columns: "ID", "x", and "y".
### The "ID" column should contain a unique ID for each animal.
### The "x" column should contain the x coordinates and the "y" column should contain
### the y coordinates of each animal's location(UTMs in meters).
habitat.units.df, 
### A data frame containing the center points for each sampling units
### There should be at least 2 columns: "x" and "y".
### The "x" column should contain the x coordinates and the "y" column should contain
### the y coordinates of the sampling unit center points (UTMs in meters).
### Additional covariate columns corresponding to each sampling unit center point should be included.
habitat.unit.size
### The radius in meters of each sampling unit.
){

  # read in animal locations
	 locs = locations.df
	 
  # read in sampled habitat units
	 hab = habitat.units.df
	 
  # remove habitat units with missing observations
  missing.data = rep(0, nrow(hab))
  for(i in 1:nrow(hab)){
      for(j in 1:ncol(hab)){
          if(is.na(hab[i,j]) == T) missing.data[i] = missing.data[i] + 1
      }
  }
     
  hab = hab[missing.data == 0,]
  
  # Prepare output file
  UDsample = NULL
  ids = unique(locs$ID)    
  for(i in 1:length(ids)){
  
      UDsample = rbind(UDsample, data.frame(hab, animal.ID=ids[i], n.locations=0, total=0))
  
  }
  
  # loop through habitat units and each animal to count the number of locations
  #   in each habitat unit.
  for(i in 1:nrow(UDsample)) {
  
      # temp1 contains only locations from 1 animal
      temp1 = locs[locs$ID == UDsample$animal.ID[i],]
  
      # calculate distance of each animal location to the center of habitat unit
      distance = sqrt((UDsample$x[i] - temp1$x)^2 + (UDsample$y[i] - temp1$y)^2)
  
      # temp2 contains locations within desired distance of center of habitat unit
      temp2 = temp1[distance <= habitat.unit.size,]
  
      # count number of locations within desired distance
      if(is.na(nrow(temp2)) == F) UDsample$n.locations[i] = nrow(temp2)
  
  }
  
  
  # count total number of locations in study area for each animal
  for(i in 1:length(ids)){
      
      total = nrow(locs[locs$ID == ids[i],])
      
      UDsample$total[UDsample$animal.ID == ids[i]] = total
  
  }
  
  return(UDsample)
  
})

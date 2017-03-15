############################################################
############################################################
###
### k-Nearest Neighbor Convex Hull (k-NNCH) Implementation
### 
### Locates homeranges of animals using the powerful k-NNCH
### algorithm.
###
### There is also a GUI tcltk application that can control
### this script. It can be invoked with "homerange<-locoh()".
###
### An academic paper on the k-NNCH (aka Fixed k LoCoH) algorithm is located at
### http://www.cnr.berkeley.edu/%7Egetz/Reprints04/Getz&WilmersEcoG_SF_04.pdf
###
### The NNCH function can be used to trigger Fixed k LoCoH,
### Fixed r LoCoH, or Adaptive LoCoH, depending on whether
### you pass it a k, r, or d parameter respectively.
###
### Sample usage:
###   library(adehabitat)
###   data(puechabon)
###   xys<-puechabon$locs[,c("X","Y")]
###   ids<-puechabon$locs[,c("Name")]
###   homerange<-NNCH(xys,id=ids,k=c(5,10,15))
###   plot(homerange)
###   plot(NNCH.k.area(homerange))
###
############################################################
############################################################


############################################################
## NNCH: Analyzes data and returns a homerange object.
## Useful Parameters:
##   xy: The x-y pairings of subject sightings
##   id: A tag of each point placing it into some form of
##       group. Can be the animal name, time of observation
##       etc. Seperate analysis will be run for each id.
##   k: The number of neighboring points out of which to
##      create convex hulls. A good value is the square root
##      of the data points. Can be a vector, in which case
##      we will repeatedly analize the data using each k.
##   r: The radius for use in the Fixed r LoCoH method.
##   a: The distance for use in the Adaptive LoCoH method.
##   min.k: If the generated k using Fixed r or Adaptive\
##   methods is less than min.k, set k equal to min.k. 
##   status: If true, prints out status messages
##   duplicates: Either a number for the distance to
##      displace duplicated points or 'delete' to delete
##		duplicated points, or 'ignore' to not worry about it
##	 hog.limit: If the number of data points in the data set
##		are above hog.limit, slow, but memory effecient,
##      methods will be used to analyze the homerange.
############################################################


NNCH<-function(xy, id=NULL, k=c(10), unin = c("m", "km"), unout = c("m2", "ha", "km2"), status=FALSE, duplicates=1, hog.limit=500, r=NULL, a=NULL, e=c(1), min.k=NULL, max.k=NULL)
{
	if (status)
		cat('Beginning homerange analysis on ',date(),'.\n',sep='')
		
	## figure out what method we're using
	
	if(is.null(a) && is.null(r)){
		mode<-'Fixed k'
		a<-1
		r<-1
		e<-1
	}else if (is.null(a)){
		mode<-'Fixed r'
		k<-1
		a<-1
	}else if(is.null(r)){
		mode<-'Adaptive'
		r<-1
		k<-1
	}else{
		stop('Could not determine what algorithm to use. Make sure to pass in only one of the following parameters: k (use Fixed k LoCoH), r (use Fixed r LoCoH), a (use Adaptive LoCoH).')
	}
	
	if (status)
		cat('Using ',mode,' LoCoH mode.\n',sep='')
	
	## make sure the data is valid
	
	k<-vectorize.parameter(k)
	a<-vectorize.parameter(a)
	r<-vectorize.parameter(r)
	e<-vectorize.parameter(e)
	
	if (ncol(xy) != 2) 
        stop("xy should have two columns")
    #set up the ids
    if (is.null(id)) 
        id <- rep(1, nrow(xy))
        
    if (nrow(xy) != length(id)) 
        stop("id should have the same length as xy")
    
        
    id <- id[!is.na(xy[, 1])]
    xy <- xy[!is.na(xy[, 1]), ]
    id <- id[!is.na(xy[, 2])]
    xy <- xy[!is.na(xy[, 2]), ]
    
    #remove duplicates in the data if duplicate set to unique
   
	if (duplicates=='delete'){
		xyid<-unique(cbind(xy,id))
		xy<-xyid[,1:2]
		id<-xyid[,3]
	}else if (is.numeric(duplicates)){
		dups<-which(duplicated(cbind(xy,id))==TRUE)
		unit<-duplicates
		for(p in dups){
			theta<-runif(1)*2*3.14
			xy[p,1]<-xy[p,1]+cos(theta)*unit
			xy[p,2]<-xy[p,2]+sin(theta)*unit
		}
	}else if (duplicates!='ignore'){
		stop("duplicates must either be a number for point displacement, or \'ignore\', or \'delete\'")
	}
    id <- factor(id)
  
	if(mode=='Fixed k'){
    	min.points<-min(sapply(levels(factor(id)),function(x) length(which(x==id)))) # finds the smallest number of points of anything uniquely deliminated by id
	    if (TRUE %in% (k > min.points)) 
	        stop(paste("too large a number for k, you can't have more k's than you have points in the smallest grouping (",min.points,")",sep=''))
	    if (TRUE %in% ( k < 3)) 
	        stop("too small a number for k, you can't have a smaller k than are needed for a triangle (3)")
    }else{
    	if(TRUE %in% (a<=0))
    		stop('a must be a positive number')
		if(TRUE %in% (e<=0))
	    		stop('e must be a positive number')
    	if(TRUE %in% (r<=0))
    		stop('r must be a positive number')
    }
    
	if (TRUE %in% (nrow(xy) < 5)) 
        stop("at least 5 locations are required to fit a homerange")
        
	#make sure data is of the correct type
	class(xy[,1]) <- "double"
	class(xy[,2]) <- "double"
	unin <- match.arg(unin)
	unout <- match.arg(unout)
	
	memory.hog<-(nrow(xy)<=hog.limit) #be quick and memory grubbing if we don't have many data points
	
	## end data setup     	
	
	#make sure we have the neccesary libraries
	if (!require(gpclib))
		stop("package gpclib required")

	## begin analysis
	
	res<-list()

	if (status)
		cat('Parameters created. Beginning data analysis...\n')
		
	for (kk in 1:nlevels(id)){						#iterate through each category
		
		if (status)
			cat('Analying data for id=',levels(id)[kk],'...\n',sep='')
			
		xyt<-xy[id==levels(id)[kk],]			#select data points for current individual
		npoints<-nrow(xyt)
		
		idt<-1:npoints
		if (memory.hog){
			dij<-as.matrix(dist(xyt))			#calculates the distance between the points
		}else{
			xproj<-xyt[,1]						#sort the points by x locations 
			xprojOrder=order(unlist(xproj))
			xyt<-xyt[xprojOrder,]
			xproj<-xproj[xprojOrder]
		}
		
		for(eVal in e){							#iterate though each e (Adaptive LoCoH)
		for(aVal in a){							#iterate though each d (Adaptive LoCoH)
		for(rVal in r){							#iterate through each r (Fixed r LoCoH)
		for(kVal in k){							#iterate through each k (Fixed k LoCoH)
			
			if (status){
				if(mode=='Fixed k'){
					cat('Analying data for k=',kVal,'...\n',sep='')
				}else if(mode=='Fixed r'){
					cat('Analying data for r=',rVal,'...\n',sep='')
				}else{
					cat('Analying data for a=',aVal,' e=',eVal,'...\n',sep='')
				}
			}
			
			#set up some empty variables
			li<-list()
			li2<-list()
			lin<-list()
			lin2<-list()
			ar<-0
			polyI<-0
			
			if(memory.hog){
				for (i in 1:npoints) {	
					if(mode=='Fixed k'){
						maxi<-kVal
					}else{
						if(mode=='Fixed r'){
							sdij<-sort(dij[i,])
							maxi<-rev(which(sdij<=rVal))[1]
						}else{
							sdij<-cumsum(sort(dij[i,]**eVal))
							maxi<-rev(which(sdij<=aVal))[1]
						}
						maxi<-min(max(maxi,min.k),max.k)
					}
					if (maxi>2){
						polyI<-polyI+1
						iid<-idt[order(dij[i,])][1:maxi]	#returns the indexes of the k nearest neighbors to point i
						xytmp<-xyt[iid,]					#converts the indexes into actual points
						ch<-chull(xytmp[,1], xytmp[,2])		#returns the indices of the boundary points
						li[[polyI]]<-as(xytmp[ch,], "gpc.poly")	#converts these indices into points and classes them as a polygon
						lin[[polyI]]<-iid						#saves the indexes of the nearest neighbors
					
					}
				}
			}else{
				for (i in 1:npoints) {
					xdist<-abs(xproj-xproj[i])
					
					curr.point<-xyt[i,]
					
					if((mode=='Fixed k') || (!is.null(min.k))){
						
						if(!is.null(min.k)){	#if are in adaptive mode and we have set a minum k
							kVal<-min.k
						}
						
						first.indices<-max(1,i-(kVal-1)):min(npoints,i+(kVal-1))
						test.points<-xyt[first.indices,]
					
						ntpoints<-nrow(test.points)
						dists<-array(dim=ntpoints)
						for (pointID in 1:ntpoints) {
							dists[pointID]<-((test.points[pointID,1]-curr.point[1])**2+(test.points[pointID,2]-curr.point[2])**2)**.5
						}
					
						first.dists<-dists
					
						maxdist<-max(dists)
						
					}else{
						first.indices<-c(i)
						first.dists<-c(0)
						maxdist<-0
					}
					
					if(mode=='Adaptive'){
						maxdist<-max(aVal**(1/eVal),maxdist)
					}else if(mode=='Fixed r'){
						maxdist<-max(rVal,maxdist)
					}
				
					test.indices<-which(xdist<=maxdist)
					
					test.indices<-setdiff(test.indices,first.indices)
					test.points<-xyt[test.indices,]
					
					ntpoints<-max(nrow(test.points),0)
					dists<-array(dim=ntpoints)
					if(ntpoints>0){
						for (pointID in 1:ntpoints) {
							dists[pointID]<-((test.points[pointID,1]-curr.point[1])**2+(test.points[pointID,2]-curr.point[2])**2)**.5
						}
					}
					test.indices<-c(first.indices,test.indices)
					
					orderDists<-order(c(first.dists,dists))
					dists<-sort(c(first.dists,dists))
					
					if(mode=='Fixed k'){
						maxi<-kVal
					}else{
						if(mode=='Fixed r'){
							maxi<-rev(which(dists<=rVal))[1]
						}else{
							sdij<-cumsum(dists**eVal)
							maxi<-rev(which(sdij<=aVal))[1]
						}
						maxi<-min(max(maxi,min.k),max.k)
					}
					if (maxi>2){
						polyI<-polyI+1
						iid<-test.indices[orderDists][1:maxi]	#returns the indexes of the k nearest neighbors to point i
						xytmp<-xyt[iid,]						#converts the indexes into actuall points
						ch<-chull(xytmp[,1], xytmp[,2])			#returns the indices of the boundary points
						li[[polyI]]<-as(xytmp[ch,], "gpc.poly")	#converts these indices into points and classes them as a polygon
						lin[[polyI]]<-iid						#saves the indexes of the nearest neighbors
					}
				}
			}
			
			if(!length(li)>0)
				stop('No hulls were created. If you are using Fixed r or Adaptive LoCoH, try entering a bigger value for r or a, or set min.k to a value greater than 2. If you are using Fixed LoCoH, you should never see this message.')
			
			aa<-unlist(lapply(li, area.poly))		#saves the area of the convex hulls
			
			
			ord<-order(aa)
			li<-li[ord]						#order hullls by size
			lin<-lin[ord]
			idbis<-idt[ord]
				
			if(mode!='Fixed k'){			#order hulls by k
				nump<-sapply(lin,length)
				ord<-rev(order(nump))
				li<-li[ord]			
				lin<-lin[ord]
				idbis<-idt[ord]
			}
			
			li2[[1]]<-li[[1]]						#this code creates a list of polygons from the smallest polygon with the union of consecutive polygons so the last item is the union of all polygons
			lin2[[1]]<-lin[[1]]						#lin2 is the total number of nearest neighbors encapsulated in that level of the convex hull creation
			for (i in 2:length(li)) {
				li2[[i]]<-union(li2[[i-1]], li[[i]])
				lin2[[i]]<-unique(c(lin2[[i-1]], lin[[i]]))
			}
			rm(li, lin)
			
			n<-unlist(lapply(lin2, length))/nrow(xyt)	#n is the percentage of points used at the  level
	
			ar<-unlist(lapply(li2, area.poly))	 		#ar contains the areas of all the polygons
			rm(lin2)

			#handle unit conversion
			if (unin == "m") {
				if (unout == "ha") 
					ar <- ar/10000
				if (unout == "km2") 
					ar <- ar/1e+06
			}
			if (unin == "km") {
				if (unout == "ha") 
					ar <- ar * 100
				if (unout == "m2") 
					ar <- ar * 1e+06
			}
			
			names(li2) <- round(n * 100)
			area <- data.frame(levels = round(n * 100, 2), area = ar)
			
			rm(ar)
			
        	dup <- !duplicated(area)
        	area = area[dup, ]
        	row.names(area) <- 1:nrow(area)

			if(mode=='Fixed k')
				name<-paste(levels(id)[kk],'.k',kVal,sep='')
			else if (mode=='Fixed r')
				name<-paste(levels(id)[kk],'.r',rVal,sep='')
			else
				name<-paste(levels(id)[kk],'.a',aVal,'.e',eVal,sep='')
				
        	res[[name]] <- list(area = area, polygons = li2[dup], xy = xyt, xy=xyt,k=kVal,id=levels(id)[kk], r=rVal, a=aVal, e=eVal)
			}
		}
		}
		}
		}
	
	if (status)
		cat('Analysis finished. Saving data...\n')
		
	#save the data and return it
	attr(res,'mode')<-mode
	if(mode=='Fixed k'){
		attr(res,'k')<-k
	}else if(mode=='Fixed r'){
		attr(res,'r')<-r
	}else{
		attr(res,'a')<-a
		attr(res,'e')<-e
	}
	attr(res, "units") <- unout
	attr(res,'id')<-id
	attr(res,'min.k')<-min.k
	attr(res,'max.k')<-max.k
	class(res)<-"NNCH"
	
	if (status)
		cat('Homerange analysis completed succesfully on ',date(),'.\n',sep='')
		
	return(res)
}

############################################################
## NNCH.iso.index: Finds a given isopleths index, given a
## percent.
############################################################

NNCH.iso.index<-function(x,percent){
	i<-max(which(x$area$levels<=percent))
	
	if(!is.integer(i)){
		warning(paste(percent,'% isopleth could not be created. More data points are probably needed in order to generate this isopleth.',sep=''))
	}
	
	return(i)
}

############################################################
## noholes.poly: Takes a gpclib polygon and removes the
## holes without changing its shape. It does this by
## connecting each hole to the outside of the polygon.
## Thus the perimeter of the polygon is changed, while area
## remains constant. Calling 'union' on a 'noholed'
## polygon will undo the 'noholeing'. Useful for plotting +
## exporting to polygon formats that do not support holes.
############################################################

noholes.poly<-function(x){

	#set up our variables
	#split apart the polygon class into its components
	
	pts<-get.pts(x)
	xvals<-lapply(pts, function(x) x$x)
	yvals<-lapply(pts, function(x) x$y)
	holes<-unlist(lapply(pts, function(x) x$hole))
	
	#identify which polygon a hole is in

	container<-array(dim=length(holes))
	hp<-which(holes==TRUE) #holes
	pp<-which(holes==FALSE)	#non-holes
	for (k in hp){
		for (i in pp){
			pa<-as(data.frame(x=xvals[[k]],y=yvals[[k]]),'gpc.poly')
			pb<-as(data.frame(x=xvals[[i]],y=yvals[[i]]),'gpc.poly')
	
			if (area.poly(union(pa,pb)) - (area.poly(pa) + area.poly(pb)) < 0) { #then pa and pb are on each other
				container[k]<-i
				break
			}
		}
	}
	
	newPoly<-new('gpc.poly') #create our new polygon
	
	#iterate through all our non-hole polygons
	#and connect its holes to its outside so it is
	#now one big continues contour
	
	for(v in 1:length(pp)){
		i=pp[v]
		
		#the points in the non hole polygon
		pxs<-unlist(xvals[[i]])
		pys<-unlist(yvals[[i]])
		
		#hp contains the indexs of all holes in this polygon
		hp<-which(container==i)
		
		#iterate through each of the holes
		while(length(hp)>0){
			sizeP<-length(pxs)
		
			xs<-xvals[hp]
			ys<-yvals[hp]
			
			#these variables make it easier to yank back the index of the hole
			#when we only have the index of the points in a list of all the points of all the holes
			sizeH<-sapply(xs, function(x) length(x))
			sumSizeH<-c(0,cumsum(sizeH))
			holeIDs<-array(dim=0)
			for(k in 1:length(sizeH)){
				holeIDs<-c(holeIDs,rep(k,sizeH[k]))
			}
			
			#calculate the distance between everypoint in the polygon and everypoitns in the holes
			#The polygons points are indexed by the x axis of the matrix
			#the whole points indexed by the y axis of the matrix
			#then we find the mininum distance between the polygon and a whole point
			
			allxs<-unlist(xs)
			allys<-unlist(ys)
			dists<-matrix(ncol=sizeP, nrow=length(allxs))
			mins<-array(dim=sizeP)
			minsIndex<-array(dim=sizeP)
			for(k in 1:sizeP){
				for(j in 1:length(allxs)){
					dists[j,k]<-((pxs[k]-allxs[j])**2+(pys[k]-allys[j])**2)**.5
				}
				minsIndex[k]<-which.min(unlist(dists[,k]))
				mins[k]<-unlist(dists[,k])[minsIndex[k]]
			}

			minP<-which.min(mins)
			minH<-minsIndex[minP]
			
			#lets find our hole index now that we have the point index
			holeID<-holeIDs[minH]
			minH<-minH-sumSizeH[holeID]
			
			holeID<-hp[holeID]

			#lets insert the hole's vertices into the non-hole polygon's vertices
			
			xpts<-unlist(xvals[holeID])
			ypts<-unlist(yvals[holeID])
			pxs<-c(pxs[1:minP],xpts[minH:length(xpts)],xpts[1:minH],pxs[minP:length(pxs)])
			pys<-c(pys[1:minP],ypts[minH:length(xpts)],ypts[1:minH],pys[minP:length(pys)])
			
			#lets remove the hole and go back and do the process again if we need to
			hp<-setdiff(hp,holeID)
			
		}
		
		#merge the polygons together
		#don't use 'union' it will undo the deholing
		newPoly<-append.poly(newPoly,as(data.frame(x=pxs,y=pys),"gpc.poly"))
	}
	
	return(newPoly)
}

############################################################
## print: Prints out some general information about the 
## homerange object.
############################################################

print.NNCH<-function(x)
{
	cat("***********************************************\n")
	cat("***   LoCoH: Nearest-Neighbor Convex Hull\n")
	cat("***********************************************\n")
	cat(paste("Homerange generated using",attr(x,'mode'),"mode available for", length(x), "groupings:\n"))
	print(names(x))
	cat("\n\nEach grouping is a component of the object. For each grouping,")
	cat("\nthe following information is available:\n")
	cat("\n$area:		home-range size estimated at various levels")
	cat("\n$polygons:	objects of class \"gpc.poly\" storing the home-range limits")
	cat("\n$hulls:		hobjects of class \"gpc.poly\" storing hulls")
	cat("\n$xy:		the relocations")
	cat("\n$id:		the id of the factor")
	if(attr(x,'mode')=='Fixed k'){
		cat("\n$k:		the value of k used to generate the homeranges")
	}else if(attr(x,'mode')=='Fixed r'){
		cat("\n$r:		the value of r used to generate the homeranges")
	}else{
		cat("\n$a:		the value of a used to generate the homeranges")
		cat("\n$e:		the value of e used to generate the homeranges")
	}
	cat('\n\nFurthermore, the following attributes are available about this object:')
	cat("\nmode:		the mode (Fixed, Adaptive FSI, or Adaptive ASI) of the LoCoH algorithm")
	cat("\nunits:		the output units, only valid if input units directly specified on input")
	cat("\nid:		the id of the factors")
	if(attr(x,'mode')=='Fixed k'){
		cat("\nk:		the values of k used to generate the homeranges")
	}else if(attr(x,'mode')=='Fixed r'){
		cat("\nr:		the values of r used to generate the homeranges")
		cat("\nmin.k:		the value of min.k used to generate the homerange")
		cat("\nmax.k:		the values of max.k used to generate the homeranges")
	}else{
		cat("\na:		the values of a used to generate the homeranges")
		cat("\ne:		the values of e used to generate the homeranges")
		cat("\nmin.k:		the value of min.k used to generate the homerange")
		cat("\nmax.k:		the values of max.k used to generate the homeranges")
	}
	cat("\n")
}



############################################################
## NNCH.area: Returns a data frame containt information
## about how much are the home range covers vs the isopleth
## level (percent). Ie, the number of points contained in
## that isopleth.
############################################################

NNCH.area<-function(x, percent=c(100,90,80,70,60,50,40,30,20,10),id=NULL,k=NULL,r=NULL,a=NULL,e=NULL)
{
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
		
	percent<-rev(vectorize.parameter(percent))

	homerange<-NNCH.select(x,k=k,id=id,a=a,r=r,e=e)

	res<-matrix(0,nrow=length(percent), ncol=length(homerange))

	for (kk in 1:length(homerange)) {
		for (i in 1:length(percent))
			res[i,kk]<-homerange[[kk]]$ar[NNCH.iso.index(homerange[[kk]],percent[i]),2]
	}
	res<-as.data.frame(res)
	row.names(res)<-percent
	names(res)<-names(homerange)
	class(res) <- c("isoarea", "data.frame")
	attr(res, "units") <- attr(homerange, "units")
	return(res)
}


############################################################
## plot.NNCH: Creates a graphical representation of the
## homerange object.
############################################################

plot.NNCH<-function(x, add.points=TRUE, pch=21, bgpts="white", colpts="black", cex=0.7, add=FALSE, same4all=TRUE, border = NA, percent=c(100,90,80,70,60,50,40,30,20,10), gr=grey(vectorize.parameter(percent)/100), id=NULL, k=NULL,r=NULL,a=NULL,e=NULL,...)
{
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
		
	percent<-rev(vectorize.parameter(percent))
	
	homerange<-NNCH.select(x,id=id,k=k,a=a,r=r,e=e)
	
	if (length(x) > 1) {
       opar <- par(mfrow = n2mfrow(length(x)))
       on.exit(par(opar))
    }
	if (same4all) {
		xxx<-do.call("rbind", lapply(homerange, function(x) x$xy))
		rx<-range(xxx[,1])
		ry<-range(xxx[,2])
	}
	for (kk in names(homerange)) {
		if (!same4all){ 
			rx<-range(homerange[[kk]]$xy[,1])
			ry<-range(homerange[[kk]]$xy[,2])
		}
		if(!add){
			if (length(homerange) > 1) 
                plot(homerange[[kk]]$xy, ty = "n", asp = 1, main = kk,  xlim = rx, ylim = ry, ...)
            if (length(homerange) == 1) 
                plot(homerange[[kk]]$xy, ty = "n", asp = 1, xlim = rx, ylim = ry, ...)
		}
		
		li2<-homerange[[kk]]$polygons
		for (i in 1:length(percent)){
			isoIndex<-NNCH.iso.index(homerange[[kk]],percent[i])
			if(is.integer(isoIndex)){
	 			plot(noholes.poly(li2[[isoIndex]]), poly.args = list(col = gr[i], border = border), add = TRUE)
	 		}
	 	}
		if (add.points)
			points(homerange[[kk]]$xy, pch=pch, bg=bgpts, col=colpts, cex=cex)
	}
}

############################################################
## NNCH.select: Selects a subset of the homerange object.
## Can be used to isolate one (or more) individual, or
## specific k, d, or r value.
############################################################

NNCH.select<-function(x, id=NULL, k=NULL, r=NULL, a=NULL, e=NULL){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	
	k<-vectorize.parameter(k)
	a<-vectorize.parameter(a)
	r<-vectorize.parameter(r)
	e<-vectorize.parameter(e)
	

	new.range=list()
	class(new.range)<-"NNCH"
	
	for (kk in names(x)) {
		if((is.null(k)) || (x[[kk]]$k %in% k)){
			if((is.null(id)) || (x[[kk]]$id %in% id)){
				if((is.null(a)) || (x[[kk]]$a %in% a)){
					if((is.null(r)) || (x[[kk]]$r %in% r)){
						if((is.null(e)) || (x[[kk]]$e %in% e)){
							new.range[[kk]]=x[[kk]]
						}
					}
				}
			}
		}
	}

	#save the data and return it
	
	attr(new.range, "units") <- attr(x, "units")
	attr(new.range, "min.k") <- attr(x, "min.k")
	attr(new.range, "max.k") <- attr(x, "max.k")
	attr(new.range,'k')<-attr(x,'k')
	attr(new.range,'a')<-attr(x,'a')
	attr(new.range,'r')<-attr(x,'r')
	attr(new.range,'e')<-attr(x,'e')
	attr(new.range,'mode')<-attr(x,'mode')
	attr(new.range,'id')<-attr(x,'id')	
	if(!is.null(k))
		attr(new.range,'k')<-attr(x,'k')[match(k, attr(x,'k'))]
	if(!is.null(a))
		attr(new.range,'a')<-attr(x,'a')[match(a, attr(x,'a'))]
	if(!is.null(r))
		attr(new.range,'r')<-attr(x,'r')[match(r, attr(x,'r'))]
	if(!is.null(e))
		attr(new.range,'e')<-attr(x,'e')[match(e, attr(x,'e'))]
	if(!is.null(id))
		attr(new.range,'id')<-attr(x,'id')[match(id, attr(x,'id'))]
	
	return(new.range)
}

############################################################
## NNCH.n.area: Creates a data.frame holding the area 
## occupied by the homerange vs the value of n used to
## generate it. n is the independent variable 'k','a','r'
############################################################

NNCH.n.area<-function(x,n, percent=c(100,50),id=NULL,k=NULL,r=NULL,a=NULL,e=NULL){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")

	percent=rev(vectorize.parameter(percent))
	
	homerange<-NNCH.select(x,k=k,id=id,a=a,r=r,e=e)

	ns<-unique(unlist(lapply(homerange,function(x) x[[n]])))
	ids<-unique(unlist(lapply(homerange,function(x) x$id)))
	
			
	result<-list()
	for(i in 1:length(ids)){
		
		cols<-length(percent)
		rows<-length(ns)
	
		res<-matrix(0,nrow=rows, ncol=cols)
	
		for (c in 1:cols) {
			for (r in 1:rows){
				shomerange<-NNCH.select(homerange,id=ids[i])
				IsoIndex<-NNCH.iso.index(shomerange[[r]],percent[c])
				res[r,c]<-shomerange[[r]]$ar[IsoIndex,2]
			}
		}
		
		res<-as.data.frame(res)
		row.names(res)<-ns
		names(res)<-percent
		class(res) <- c("data.frame")
		attr(res, "units") <- attr(homerange, "units")
		
		result[[ids[i]]]<-res
	}
	attr(result, "variable") <- n
	attr(result, "units") <- attr(homerange, "units")
	class(result)<-c("narea","list")
	
	return(result)
}

NNCH.k.area<-function(x,...){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	if (attr(x,"mode")!='Fixed k')
		stop("homerange must have been generated using the Fixed k algorithm")
	return(NNCH.n.area(x,n='k',...))
}

NNCH.a.area<-function(x,...){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	if (attr(x,"mode")!='Adaptive')
		stop("homerange must have been generated using the Adaptive algorithm")
		
	return(NNCH.n.area(x,n='a',...))
}

NNCH.e.area<-function(x,...){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	if (attr(x,"mode")!='Adaptive')
		stop("homerange must have been generated using the Adaptive algorithm")
		
	return(NNCH.n.area(x,n='e',...))
}

NNCH.r.area<-function(x,...){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	if (attr(x,"mode")!='Fixed r')
		stop("homerange must have been generated using the Fixed r algorithm")
	
	return(NNCH.n.area(x,n='r',...))
}

NNCH.ae.area<-function(x, id=NULL,percent=100,a=NULL,e=NULL){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	if (attr(x,"mode")!='Adaptive')
		stop("homerange must have been generated using an Adaptive algorithm")
	if(length(percent)>1)
		stop("percent for this method can only be a scalar")

	homerange<-NNCH.select(x,id=id,a=,e=e)

	ds<-unique(unlist(lapply(homerange,function(x) x[['a']])))
	es<-unique(unlist(lapply(homerange,function(x) x[['e']])))
	ids<-unique(unlist(lapply(homerange,function(x) x$id)))
			
	result<-list()
	for(i in 1:length(ids)){
		
		cols<-length(es)
		rows<-length(ds)
	
		res<-matrix(0,nrow=rows, ncol=cols)
	
		for (c in 1:cols) {
			for (r in 1:rows){
				shomerange<-NNCH.select(homerange,id=ids[i],e=es[c])
				
				res[r,c]<-shomerange[[r]]$ar[NNCH.iso.index(shomerange[[r]],percent),2]
			}
		}
		
		res<-as.data.frame(res)
		row.names(res)<-ds
		names(res)<-es
		class(res) <- c("data.frame")
		attr(res, "units") <- attr(homerange, "units")
		
		result[[ids[i]]]<-res
	}
	attr(result, "percent") <- percent
	attr(result, "units") <- attr(homerange, "units")
	class(result)<-c("aearea","list")
	
	return(result)
}

############################################################
## NNCH.der.n.area: Basically, the derivative of the
## NNCH.n.area graph.
############################################################

NNCH.der.n.area<-function(x, ...){
	areas<-NNCH.n.area(x,...)

	origns<-row.names(areas[[1]])
	if(length(origns)<2){
		stop('You must have at least 2 values, if you want to generate the derivative of n vs area data.')
	}
	
	result<-list()
	for(i in 1:length(areas)){
		ns<-origns[2:length(origns)]
		singlearea<-areas[[i]]
		percents<-names(singlearea)
		res<-matrix(0,nrow=length(ns), ncol=length(percents))
		for(kk in 1:length(percents)){
			res[,kk]<-diff(singlearea[[kk]])/diff(as.numeric(origns))
		}
		res<-as.data.frame(res)
		row.names(res)<-ns
		names(res)<-percents
		
		class(res) <- "data.frame"
		attr(res, "units") <- attr(x, "units")
		result[[names(areas)[i]]]<-res
	}

	class(result) <- c("narea", "list")
	attr(result, "units") <- attr(areas, "units")

	return(result)
}

NNCH.der.k.area<-function(x, ...){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	if (attr(x,"mode")!='Fixed k')
		stop("homerange must have been generated using the Fixed k algorithm")
		
	return(NNCH.der.n.area(x,n='k',...))
}

NNCH.der.r.area<-function(x, ...){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	if (attr(x,"mode")!='Fixed r')
		stop("homerange must have been generated using the Fixed r algorithm")
		
	return(NNCH.der.n.area(x,n='r',...))
}

NNCH.der.a.area<-function(x, ...){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	if (attr(x,"mode")!='Adaptive')
		stop("homerange must have been generated using the Adaptive algorithm")
		
	return(NNCH.der.n.area(x,n='a',...))
}

NNCH.der.e.area<-function(x, ...){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")
	if (attr(x,"mode")!='Adaptive')
		stop("homerange must have been generated using the Adaptive algorithm")
		
	return(NNCH.der.n.area(x,n='e',...))
}

############################################################
## plot.narea: Useful graph formatter.
############################################################

plot.narea<-function(x,...){
	if (!inherits(x, "narea")) 
		stop("should be of class narea")
	
	opar <- par(mfrow = n2mfrow(length(x)),mar=c(3.2,3.2,1.6,0.5),mgp=c(2.1, 1, 0))
	on.exit(par(opar))
	
	for(i in 1:length(x)){
    	matplot(row.names(x[[i]]),x[[i]],type="l", xlab = paste("Value of",attr(x,'variable')), ylab = "Homerange Area", main=names(x)[i], ...)  
    	matpoints(row.names(x[[i]]),x[[i]],type="p",pch = 16, cex = 0.6)
	}
}

############################################################
## plot.aearea: Useful graph formatter.
############################################################

plot.aearea<-function(x,...){
	if (!inherits(x, "aearea")) 
		stop("should be of class aearea")
	if(!require(ade4))
		stop("ade4 library required")
	

	opar <- par(mfrow = n2mfrow(length(x)),mar=c(3.2,3.2,1.6,0.5),mgp=c(2.1, 1, 0))
	on.exit(par(opar))

	for(i in 1:length(x)){
		curr<-x[[i]]
    	xv<-array()
    	yv<-array()
    	a<-array()
    	for(r in 1:nrow(curr)){
    		for(c in 1:ncol(curr)){
    			xv<-c(xv,names(curr)[c])
    			yv<-c(yv,row.names(curr)[r])
    			a<-c(a,curr[r,c])
    		}
    	}
    	a<-as.numeric(a[2:length(a)])
    	a<-formatC(a, format="g", digits=3)
    	xv<-as.numeric(xv[2:length(xv)])
    	yv<-as.numeric(yv[2:length(yv)])

    	rx<-range(xv)
		ry<-range(yv)
		
		plot(xv,yv,xlab ="Value of e", ylab = "Value of a", main=names(x)[i],xlim=rx,ylim=ry,...)

		scatterutil.eti(xv, yv, a, .7)
			
	}
}

############################################################
## plot.isoarea: Useful graph formatter.
############################################################

plot.isoarea<-function(x,...){
	if (!inherits(x, "isoarea")) 
		stop("should be of class isoarea")
	opar <- par(mfrow = n2mfrow(ncol(x)),mar=c(3.2,3.2,1.6,0.5),mgp=c(2.1, 1, 0))
	on.exit(par(opar))
	for (i in 1:ncol(x)) {
		plot(as.numeric(row.names(x)), x[, i], main = names(x)[i], 
		pch = 16, cex = 0.5, xlab = "Isopleth Level (%)", ylab = "Isopleth Area")
	}

}

############################################################
## summary.NNCH: Prints a summary of the homerange in text
## file format. Contains useful statistics about the range.
############################################################

summary.NNCH<-function(x,file='',id=NULL,k=NULL,r=NULL,a=NULL,e=NULL){

	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")

	homerange<-NNCH.select(x,k=k,id=id,a=a,r=r,e=e)
	
	if (file!='')
		sink(file=file)

	cat('#########################################\n','Homerange Analysis Summary\n','Printed on ',date(),'\n#########################################\n\n',sep='')
	
	cat('***LoCoH Mode:\n')
	print(attr(homerange,'mode'))
	
	cat('\n\n***ID\'s Analyzed:\n')
	print(levels(attr(homerange,'id')))
	
	if(attr(homerange,'mode')=='Fixed k'){
		cat('\n\n***k\'s Analyzed:\n')
		print(attr(homerange,'k'))
	}else if(attr(homerange,'mode')=='Fixed r'){
		cat('\n\n***r\'s Analyzed:\n')
		print(attr(homerange,'r'))
	}else{
		cat('\n\n***a\'s Analyzed:\n')
		print(attr(homerange,'a'))
		
		cat('\n\n***e\'s Analyzed:\n')
		print(attr(homerange,'e'))
	}

	cat('\n\n***Units if Specified (if not specified, ignore this):\n')
	print(attr(homerange,'units'))
	
	cat('\n\n***Areas of Isopleths:\n')
	print(try(NNCH.area(homerange)))

	if(attr(homerange,'mode')=='Fixed k'){
		cat('\n\n***k vs Isopleth Area Data:\n')
		print(try(NNCH.k.area(homerange)))
	
		cat('\n\n***Derivative of k vs Isopleth Area Data:\n')
		print(try(NNCH.der.k.area(homerange)))
	}else if(attr(homerange,'mode')=='Fixed r'){
		cat('\n\n***r vs Isopleth Area Data:\n')
		print(try(NNCH.r.area(homerange)))
	
		cat('\n\n***Derivative of r vs Isopleth Area Data:\n')
		print(try(NNCH.der.r.area(homerange)))
	}else{
		cat('\n\n***a vs Isopleth Area Data:\n')
		print(try(NNCH.a.area(homerange)))
	
		cat('\n\n***Derivative of a vs Isopleth Area Data:\n')
		print(try(NNCH.der.a.area(homerange)))
	}

	cat('\n\nEnd Summary\n########################################\n')
	
	if (file!='')
		sink()
}

############################################################
## asciigrid.NNCH: Turns homerange data to
## ESRI's ascii grid file format. 
############################################################

NNCH.asciigrid<-function(x,cellsize=1,percent=c(100,90,80,70,60,50,40,30,20,10),id=NULL,k=NULL,r=NULL,a=NULL,e=NULL){
    
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")

	if (!require(adehabitat))
		stop("package adehabitat required")
    
	percent<-rev(vectorize.parameter(percent))

	x<-NNCH.select(homerange,k=k,id=id,r=r,a=a,e=e)
	
	if (length(homerange)!=1)
		stop("you can only export one homerange to asciigrid at a time")
	
	isoIndex<-NNCH.iso.index(homerange[[1]],percent[1])
	if(is.integer(isoIndex)){
		baseIso<-homerange[[1]]$polygons[[isoIndex]]
	}else{
		stop('your biggest isopleth, contains no data')
	}
	
	bbox<-get.bbox(baseIso)
	
	mat<-matrix(data=0,nrow=ceiling((bbox$x[2]-bbox$x[1])/cellsize)+2,ncol=ceiling((bbox$y[2]-bbox$y[1])/cellsize)+2)
	
	asc<-as.asc(mat,xll=bbox$x[1]-cellsize,yll=bbox$y[1]-cellsize,cellsize=cellsize)
	
	rasters<-list()
	for (i in 1:length(percent)) {
		isoIndex<-NNCH.iso.index(homerange[[1]],percent[i])
		if(is.integer(isoIndex)){
			xys<-get.pts(noholes.poly(homerange[[1]]$polygons[[isoIndex]]))
			subPolys<-length(xys)
			temp.asc<-list()
			for (subPoly in 1:subPolys){
				xy<-xys[[subPoly]]
				temp.asc[[subPoly]]<-mcp.rast(cbind(xy$x,xy$y),asc)
			}
			
			for(k in 1:length(temp.asc)){
				temp.asc[[k]][is.na(temp.asc[[k]])]<-0
				temp.asc[[k]][which(temp.asc[[k]]==1)]<-percent[i]
			}			
			rasters[[i]]<-temp.asc[[1]]
			if(length(temp.asc)>1){
				for(k in 2:length(temp.asc)){
					rasters[[i]]<-rasters[[i]]+temp.asc[[k]]
				}
			}
		}
	}
	
	rasters[[1]][0==rasters[[1]]]<-101 #a big number, bigger than possible
	if(length(rasters)>1){
		for (i in 2:length(rasters)){
			rasters[[i]][0==rasters[[i]]]<-101
			rasters[[1]]<-pmin(rasters[[1]],rasters[[i]])
		}
	}
	rasters[[1]][101==rasters[[1]]]<-NA
	return(rasters[[1]])
}

############################################################
## asciigrid.export.NNCH: Turns homerange data to
## ESRI's ascii grid file format. And then saves it to file.
## File should have a .asc extension 
############################################################

NNCH.export.asciigrid<-function(x,file,...){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")

	if (!require(adehabitat))
		stop("package adehabitat required")
		
    export.asc(NNCH.asciigrid(x,...),file=file)
}

############################################################
## NNCH.export.shapefile: Export of homerange data to
## ESRI's shapefile. Saves the isopleths epecified by
## percent. 'file' should have no extension, those will be
## added by the function. If the extension '.shp' exists
## it will be removed.
############################################################

NNCH.export.shapefile<-function(x,file,...){
	shpobject<-NNCH.make.shapefile(x,...)
	
	if(tolower(substr(file,nchar(file)-3,nchar(file)))=='.shp')
		file<-substr(file,0,nchar(file)-4)
		
	write.shapefile(shpobject, file, arcgis=TRUE)
}

############################################################
## NNCH.make.shapefile: Creates a shapefile using the
## shapefile package
############################################################

NNCH.make.shapefile<-function(x,percent=c(100,90,80,70,60,50,40,30,20,10),k=NULL,id=NULL,r=NULL,a=NULL,e=NULL){
	if (!inherits(x, "NNCH"))
		stop("x should be of class \"NNCH\"")

	if (!require(shapefiles))
		stop("package shapefiles required")
	
	percent=rev(vectorize.parameter(percent))
	
	homerange<-NNCH.select(x,k=k,id=id,a=a,r=r,e=e)
	
	polyID<-0
	shp.data<-matrix(ncol=3,nrow=0)
	if(attr(homerange,'mode')!='Adaptive'){
		att.data<-matrix(ncol=4,nrow=0)
	}else{
		att.data<-matrix(ncol=5,nrow=0)
	}
	
	for(kk in names(homerange)) {
		for (i in length(percent):1){
			isoIndex=NNCH.iso.index(homerange[[kk]],percent[i])
			if(is.integer(isoIndex)){
				currPoly<-noholes.poly(homerange[[kk]]$polygons[[isoIndex]])
				subPolys<-length(get.pts(currPoly))
				for (subPoly in 1:subPolys){
					polyID<-polyID+1
			 		polypoints<-data.frame(get.pts(currPoly)[[subPoly]][1:2])
			 		polypoints<-rbind(polypoints,polypoints[1,])
			 		polypoints<-cbind(id=rep(polyID,nrow(polypoints)),polypoints)
					shp.data<-rbind(polypoints,shp.data)
					if(attr(homerange,'mode')=='Fixed k'){
						att.data<-rbind(c(polyID,homerange[[kk]]$id,homerange[[kk]]$k,percent[i]),att.data)
					}else if(attr(homerange,'mode')=='Fixed r'){
						att.data<-rbind(c(polyID,homerange[[kk]]$id,homerange[[kk]]$r,percent[i]),att.data)
					}else{
						att.data<-rbind(c(polyID,homerange[[kk]]$id,homerange[[kk]]$a,homerange[[kk]]$e,percent[i]),att.data)
					}
				}
			}
		}
	}
	
	if(attr(homerange,'mode')=='Fixed k'){
		att.data<-data.frame(Id=as.integer(att.data[,1]),HR_id=att.data[,2],k=as.integer(att.data[,3]),Isopleth=as.numeric(att.data[,4]))
	}else if(attr(homerange,'mode')=='Fixed r'){
		att.data<-data.frame(Id=as.integer(att.data[,1]),HR_id=att.data[,2],r=as.integer(att.data[,3]),Isopleth=as.numeric(att.data[,4]))
	}else{
		att.data<-data.frame(Id=as.integer(att.data[,1]),HR_id=att.data[,2],a=as.numeric(att.data[,3]),e=as.numeric(att.data[,4]),Isopleth=as.numeric(att.data[,5]))
	}
	
	shp.file=convert.to.shapefile(shp.data, att.data, 'Id', 5)

	return(shp.file)
}

############################################################
##  shapefile.points: Takes the root name of a shapefile
##  and returns a list of points in the shapefile
############################################################

shapefile.points<-function(file){
	if (!require(shapefiles))
		stop("package shapefiles required")
	if(tolower(substr(file,nchar(file)-3,nchar(file)))=='.shp')
		file<-substr(file,0,nchar(file)-4)	
		
	shp.file <- read.shapefile(file)
	homerange.points <- convert.to.simple(shp.file$shp)
	rm(shp.file)
	homerange.points<-homerange.points[,c('X','Y')]
	return(homerange.points)
}

############################################################
## NNCH.shapefile: Runs the homerange analysis on an ESRI
## shapefile of points. 'file' should have no extension.
## If 'file' ends with '.shp', it will be removed.
############################################################

NNCH.shapefile<-function(file,...){
	homerange.points<-shapefile.points(file)
	
	res<-NNCH(homerange.points,...)
	class(res)<-"NNCH"

	return(res)
}

############################################################
## NNCH.textfile: Runs the homerange analysis on a
## textfile. First line defines the data and is ignored.
## Every line after that is x-coord [SPACE OR TAB] y-coord.
############################################################

NNCH.textfile<-function(file,...){
	homerange.points <- read.delim(file)
	
	homerange.points<-homerange.points[,1:2]

	res<-NNCH(homerange.points,...)
	class(res)<-"NNCH"
	return(res)
}


############################################################
## vectorize.parameter: Takes a string or vector and turns it 
## into a vector of numbers. Sorts the parameter in ascending
## order.
############################################################

vectorize.parameter<-function(x){
	if(is.null(x))
		return(x)
		
	if(is.character(x))
		x<-as.numeric(unlist(strsplit(x,' *, *')))
	return(sort(x))
}

############################################################
## NNCH.test: Unit tests to ensure the correct function of
## the script.
############################################################

NNCH.test<-function(){
	cat('Beginning tests.\nLoading data...\n')
	load("NNCHtest.Rdata")

	cat('Data loaded succesfully.\n')
	
	cat('Beginning Fixed k LoCoH homerange test...\n')
	thomerange<-NNCH(txy,id=tid,k=c(5,10,15),duplicates='delete')
	cat('Homerange test succesful.\n')
	
	cat('Beginning memory hog test...\n')
	if(NNCH.area(NNCH(txy,k=4,duplicates='delete'))==NNCH.area(NNCH(txy,k=4,duplicates='delete',hog.limit=0))){
		cat('Memory hog test successful.\n')
	}else{
		cat('***Memory hog test failed.\n')
	}
	
	cat('Beginning area test...\n')
	tarea<-NNCH.area(thomerange,percent=c(100,95,90,80,70,60,50,40,30,20,10))
	if(tarea==f.area){
		cat('Area test successful.\n')
	}else{
		cat('***Area test failed.\n')
		cat("*****Wanted:")
		print(area)
		cat("*****Got:")
		print(tarea)
	}

	cat('Beginning area vs k test...\n')
	tdata<-NNCH.k.area(thomerange,percent=100)
	if(tdata[[1]]==k.area[[1]]){
		cat('k vs area test successful.\n')
	}else{
		cat('***k vs area test failed.\n')
		cat("*****Wanted:\n")
		print(k.area)
		cat("*****Got:\n")
		print(tdata)
	}

	cat('Beginning area vs k derivative test...\n')
	tdata<-NNCH.der.k.area(thomerange, percent=100)
	if(tdata[[1]]==der.k.area[[1]]){
		cat('k vs area derivative test successful.\n')
	}else{
		cat('***k vs area derivative test failed.\n')
		cat("*****Wanted:\n")
		print(der.k.area)
		cat("*****Got:\n")
		print(tdata)
	}
	
	cat('Beginning Adaptive LoCoH homerange test...\n')
	if(NNCH.area(NNCH(txy,a=1500,duplicates='delete',min.k=3))==NNCH.area(NNCH(txy,a=1500,duplicates='delete',hog.limit=0,min.k=3))){
		cat('Adaptive test successful.\n')
	}else{
		cat('***Adaptive test failed.\n')
	}
	
	cat('Beginning Fixed r LoCoH homerange test...\n')
	if(NNCH.area(NNCH(txy,r=200,duplicates='delete',min.k=3))==NNCH.area(NNCH(txy,r=200,duplicates='delete',hog.limit=0,min.k=3))){
		cat('Fixed r test successful.\n')
	}else{
		cat('***Fixed r test failed.\n')
	}

	
	cat('Unit tests completed.\n')
}
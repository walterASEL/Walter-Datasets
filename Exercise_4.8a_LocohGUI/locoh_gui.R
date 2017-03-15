##################################################
##
## LoCoH: A tcl/tk Interface to the k-NNCH Method
## Written by Scott Fortmann-Roe; Summer 2006
##
##################################################

locoh <- function(){
    if (data.class(result<-try(locoh.run(),TRUE))=="try-error")
    {
        tkmessageBox(title="An error has occured!",message=as.character(result),icon="error",type="ok")
        return(NULL)
    }
    else
        return (result)
}

locoh.run<-function(){
	if (!require(gpclib)) 
	  stop("Package gpclib required for running of this program.")
	if (!require(tcltk))
	  stop("Package tcltk required for running of this program.")
	
	# Create a new toplevel window
	tt <- tktoplevel()
	tkwm.title(tt,"LoCoH Homerange Generator")
	tkwm.resizable(tt, 0, 0) 
	
	Top.frm <- tkframe(tt)
	
	exitStatus <- tclVar(0)
	
	#create input frame
	Input.frm <- tkframe(Top.frm,relief="groove",borderwidth=2)
	
	TextFile.txt <- tklabel(Input.frm,text='Load a tab deliminated text file of points:')
	tkgrid(TextFile.txt,columnspan=2,sticky='w')
	TextFile.frm <- tkframe(Input.frm,relief="groove",borderwidth=2)
	TextFile.val <- tclVar('')
	TextFile.txt <- tklabel(TextFile.frm,text=tclvalue(TextFile.val))
	TextFile.but <- tkbutton(TextFile.frm,text="Select Textfile",command=function() tclvalue(TextFile.val)<-tclvalue(tkgetOpenFile(filetypes="{{Textfiles} {.txt}} {{All files} *}")))
	TextFile.hlp <- tkbutton(TextFile.frm,text="?",command=function() tkmessageBox(title="Help",message="You can select an tab deliminated textfile of points for analysis. The first column are the x coordinates of the points, the second column contains the y coordinates. The firs row labels the columns('x','y'; 'easting','norhting';etc)",icon="info",type="ok"))
	tkconfigure(TextFile.txt,textvariable=TextFile.val)
	tkgrid(TextFile.txt,TextFile.but,TextFile.hlp)
	tkgrid(TextFile.frm,columnspan=2,sticky='nesw')
	
	ShapeFile.txt <- tklabel(Input.frm,text='Load an ESRI shapefile of points:')
	tkgrid(ShapeFile.txt,columnspan=2,sticky='w')
	ShapeFile.frm <- tkframe(Input.frm,relief="groove",borderwidth=2)
	ShapeFile.val <- tclVar('')
	ShapeFile.txt <- tklabel(ShapeFile.frm,text=tclvalue(ShapeFile.val))
	ShapeFile.but <- tkbutton(ShapeFile.frm,text="Select Shapefile",command=function() tclvalue(ShapeFile.val)<-tclvalue(tkgetOpenFile(filetypes="{{ESRI Shapefiles} {.shp}} {{All files} *}")))
	ShapeFile.hlp <- tkbutton(ShapeFile.frm,text="?",command=function() tkmessageBox(title="Help",message="You can select an ESRI point shapefile (extension .shp) for analysis.",icon="info",type="ok"))
	tkconfigure(ShapeFile.txt,textvariable=ShapeFile.val)
	tkgrid(ShapeFile.txt,ShapeFile.but,ShapeFile.hlp)
	tkgrid(ShapeFile.frm,columnspan=2,sticky='nesw')
	
	R.txt <- tklabel(Input.frm,text='If no textfile or shapefile is selected, load data from your R environment:')
	tkgrid(R.txt,columnspan=2,sticky='w')
	R.frm <- tkframe(Input.frm,relief="groove",borderwidth=2)
	XY.txt <- tklabel(R.frm,text="Enter the source of the X-Y data to analyze:")
	XY.val <- tclVar("")
	XY.ent <-tkentry(R.frm,width="20",textvariable=XY.val)
	CY.hlp <- tkbutton(R.frm,text="?",command=function() tkmessageBox(title="Help",message="This is a two column matrix of data where the first column are the x-coordinate of animal sightings while the second column is the y-coordinate.",icon="info",type="ok"))
	tkgrid(XY.txt,XY.ent,CY.hlp)
	ID.txt <- tklabel(R.frm,text="Enter the source of the ID data for analysis (optional):")
	ID.val <- tclVar("")
	ID.ent <-tkentry(R.frm,width="20",textvariable=ID.val)
	ID.hlp <- tkbutton(R.frm,text="?",command=function() tkmessageBox(title="Help",message="The id's are a list of data tags--one for each data point. They can be used to tag data by the individual, or by the season of sighting, or by any other parameter. Each unique id will be analyzed separately.",icon="info",type="ok"))
	tkgrid(ID.txt,ID.ent,ID.hlp)
	tkgrid(R.frm,columnspan=2,sticky='nesw')
	
	#create Analysis frame
	Analysis.frm <- tkframe(Top.frm,relief="groove",borderwidth=2)
	
	Method.frm <- tkframe(Analysis.frm,relief="groove",borderwidth=2)
	Method.txt <- tklabel(Analysis.frm,text="Use which LoCoH algorithm:")
	k.rb <- tkradiobutton(Method.frm)
	r.rb <- tkradiobutton(Method.frm)
	a.rb <- tkradiobutton(Method.frm)
	Method.val <- tclVar("k")
	tkconfigure(k.rb,variable=Method.val,value="k")
	tkconfigure(r.rb,variable=Method.val,value="r")
	tkconfigure(a.rb,variable=Method.val,value="a")
	tkgrid(tklabel(Method.frm,text="Fixed k "),k.rb,tklabel(Method.frm,text="Fixed r "),r.rb,tklabel(Method.frm,text="Adaptive "),a.rb)
	Method.hlp <- tkbutton(Analysis.frm,text="?",command=function() tkmessageBox(title="Help",message="You can either use the Fixed k LoCoH method, Fixed r LoCoH method, or Adaptive LoCoH method. If the Fixed k method is selected, the hulls will be constructed out of the (k-1) nearest neighbors for each points. If the Fixed r method is selected, hulls will be created from all points within r distance from each point. If the Adaptive method is selected, hulls will be constructed from the maximum number of nearest neighbors such that the sum of their distances from the root point is less than or equal to a.",icon="info",type="ok"))
	tkgrid(Method.txt,Method.frm,Method.hlp)
	
	Var.txt <- tklabel(Analysis.frm,text="Enter the value of the variable (k, r, or d):")
	Var.val <- tclVar("10")
	Var.ent <-tkentry(Analysis.frm,width="20",textvariable=Var.val)
	tkgrid(Var.txt,Var.ent)
	
	Dups.frm <- tkframe(Analysis.frm,relief="groove",borderwidth=2)
	Dups.txt <- tklabel(Analysis.frm,text="When encountering duplicated points:")
	Displace.rb <- tkradiobutton(Dups.frm)
	Ignore.rb <- tkradiobutton(Dups.frm)
	Delete.rb <- tkradiobutton(Dups.frm)
	Dups.val <- tclVar("displace")
	tkconfigure(Displace.rb,variable=Dups.val,value="displace")
	tkconfigure(Ignore.rb,variable=Dups.val,value="ignore")
	tkconfigure(Delete.rb,variable=Dups.val,value="delete")
	Displace.val <- tclVar("1")
	Displace.ent <-tkentry(Dups.frm,width="4",textvariable=Displace.val)
	tkgrid(tklabel(Dups.frm,text="Displace the following amount "),Displace.rb,Displace.ent)
	tkgrid(tklabel(Dups.frm,text="Delete "),Delete.rb)
	tkgrid(tklabel(Dups.frm,text="Treat as Normal "),Ignore.rb)
	Dups.hlp <- tkbutton(Analysis.frm,text="?",command=function() tkmessageBox(title="Help",message="In animal location data, especially when data are generated with a radio or GPS collar at regular time intervals, there are often duplicate points. This can cause problems when constructing local hulls because you need at least three unique points to create a hull. There are a few options for how LoCoH handles duplicate points. If the \'displace\' option is selected, duplicated points will be displaced in a random direction by the amount specified. The displacement distance should represent an approximate radius of habitat used by the animal when stationary (e.g., a patch size or the distance they monitor for their \'safety zone\'). If the \'ignore\' option is selected, duplicate points will be included when searching for the k-1 nearest neighbors. And a hull might be formed with less that three unique points, resulting in a zero-area hull. The \'delete\' option will simply exclude duplicate points from hull creation and nearest neighbor searches. For example if your dataset has 429 points, but 5 of those points lie on top of each other, then 4 of the 5 duplicates will be basically excluded from the analysis.",icon="info",type="ok"))
	tkgrid(Dups.txt,Dups.frm,Dups.hlp)
	
	Isos.txt <- tklabel(Analysis.frm,text="Enter the isopleth levels:")
	Isos.val <- tclVar("c(100,90,80,70,60,50,40,30,20,10)")
	Isos.ent <-tkentry(Analysis.frm,width="20",textvariable=Isos.val)
	Isos.hlp <- tkbutton(Analysis.frm,text="?",command=function() tkmessageBox(title="Help",message="The isopleth levels tells LoCoH for what percentages it should generate isopleths. Isopleths contain a certain number of the data points and are sorted by density. Thus the 10% isopleth is the merger of the smallest hulls that contain 10% of all the original data points. The 100% isopleth will contain all the points. If you're just trying to see which value of k works best for your data, then you probably only need to view the 100% isopleth. The format for data entry is the same as that for the value of k.",icon="info",type="ok"))
	tkgrid(Isos.txt,Isos.ent,Isos.hlp)
	
	tkgrid.configure(Dups.frm,sticky="nsew")
	
	
	#create Output frame
	Output.frm <- tkframe(Top.frm,relief="groove",borderwidth=2)
	
	Save.txt <- tklabel(Output.frm,text="Save files to:")
	Folder.val <- tclVar(getwd())
	Folder.txt <- tklabel(Output.frm,text=tclvalue(Folder.val))
	Folder.hlp <- tkbutton(Output.frm,text="?",command=function() tkmessageBox(title="Help",message="This directory is where generated graphs and data files will be saved. Defaults to the current working directory.",icon="info",type="ok"))
	tkconfigure(Folder.txt,textvariable=Folder.val)
	tkgrid(Save.txt,Folder.txt,Folder.hlp)
	
	File.but <- tkbutton(Output.frm,text="Change Save Directory",command=function() tclvalue(Folder.val)<-tclvalue(tkchooseDirectory()))
	tkgrid(tklabel(Output.frm,text=" "),File.but)
	
	ExportShape.txt <- tklabel(Output.frm,text="Save a shapefile of the isopleths:")
	ExportShape.chk <- tkcheckbutton(Output.frm)
	ExportShape.val <- tclVar("0")
	tkconfigure(ExportShape.chk,variable=ExportShape.val)
	ExportShape.hlp <- tkbutton(Output.frm,text="?",command=function() tkmessageBox(title="Help",message="This option allows you to save a shapefile (compatable with ESRI's ArcGis) of number of polygon isopleths representing the homerange.",icon="info",type="ok"))
	tkgrid(ExportShape.txt,ExportShape.chk,ExportShape.hlp)
	
	GraphHR.txt <- tklabel(Output.frm,text="Save a pdf of the homerange plot:")
	GraphHR.chk <- tkcheckbutton(Output.frm)
	GraphHR.val <- tclVar("0")
	tkconfigure(GraphHR.chk,variable=GraphHR.val)
	GraphHR.hlp <- tkbutton(Output.frm,text="?",command=function() tkmessageBox(title="Help",message="This option allows you to save a pdf image of a graphical representation of the homeranges that LoCoh generates.",icon="info",type="ok"))
	tkgrid(GraphHR.txt,GraphHR.chk,GraphHR.hlp)
	
	GraphAreaIso.txt <- tklabel(Output.frm,text="Save a pdf of an area vs isopleth plot:")
	GraphAreaIso.chk <- tkcheckbutton(Output.frm)
	GraphAreaIso.val <- tclVar("0")
	tkconfigure(GraphAreaIso.chk,variable=GraphAreaIso.val)
	GraphAreaIso.hlp <- tkbutton(Output.frm,text="?",command=function() tkmessageBox(title="Help",message="This option allows you to save a pdf image of a graphical representation of how the area covered by the homerange as we look at successive isopleths.",icon="info",type="ok"))
	tkgrid(GraphAreaIso.txt,GraphAreaIso.chk,GraphAreaIso.hlp)
	
	GraphkArea.txt <- tklabel(Output.frm,text="Save a pdf of an area vs variable (k, r, or d) plot:")
	GraphkArea.chk <- tkcheckbutton(Output.frm)
	GraphkArea.val <- tclVar("0")
	tkconfigure(GraphkArea.chk,variable=GraphkArea.val)
	GraphkArea.hlp <- tkbutton(Output.frm,text="?",command=function() tkmessageBox(title="Help",message="This option allows you to save a pdf image of a graphical representation of how the area covered by the isopleths changes as we try out different values of k, r, or d.",icon="info",type="ok"))
	tkgrid(GraphkArea.txt,GraphkArea.chk,GraphkArea.hlp)
	
	# create the button frame
	Button.frm <- tkframe(Top.frm,relief="groove",borderwidth=0)
	Help.but <- tkbutton(Button.frm,text="Help",command=function() browseURL('http://locoh.palary.org/rtutorial'))
	Spacer1.txt <- tklabel(Button.frm,text="  ")
	Cancel.but <- tkbutton(Button.frm,text="Cancel",command=function() tclvalue(exitStatus)<-2)
	Spacer2.txt <- tklabel(Button.frm,text="  ")
	OK.but <- tkbutton(Button.frm,text="Proccess",    command=function() tclvalue(exitStatus)<-1)
	tkgrid(Help.but,Spacer1.txt,Cancel.but,Spacer2.txt,OK.but)
	
	
	# Draw the window
	tkgrid(tklabel(Top.frm,text="An implementation of the powerful LoCoH methods for generating homeranges."))
	
	tkgrid(tklabel(Top.frm,text=""))
	
	tkgrid(tklabel(Top.frm,text="Data Input"),sticky="w")
	tkgrid(Input.frm)
	tkgrid.configure(Input.frm,sticky="nsew")
	
	tkgrid(tklabel(Top.frm,text="Analysis Options"),sticky="w")
	tkgrid(Analysis.frm)
	tkgrid.configure(Analysis.frm,sticky="nsew")
	
	tkgrid(tklabel(Top.frm,text="Output"),sticky="w")
	tkgrid(Output.frm)
	tkgrid.configure(Output.frm,sticky="nsew")
	
	tkgrid(tklabel(Top.frm,text=""))
	tkgrid(Button.frm)
	tkgrid.configure(Button.frm,sticky="nsew")
	
	tkgrid(Top.frm)
	
	# Capture the event "Destroy" (e.g. Alt-F4 in Windows) and when this happens, assign 2 to done.
	tkbind(tt,"<Destroy>",function() tclvalue(exitStatus)<-2)
	
	tkfocus(tt)
	
	#Wait for completion
	tkwait.variable(exitStatus)
	doneVal <- as.integer(tclvalue(exitStatus))
	
	#The Window was canceled
	if (doneVal!=1) {
		tkdestroy(tt)
		return(NULL)
	}
	
	####################################
	#The Actual Proccessing of the Data
	####################################
	
	tkconfigure(tt,cursor="watch")
	dlg <- tktoplevel()
	tkwm.deiconify(dlg)
	tkwm.resizable(dlg, 0, 0) 
	tkgrab.set(dlg)
	tkfocus(dlg)
	tkwm.title(dlg,"Processing Data...")
	tkgrid(tklabel(dlg,text="Please be patient while we proccess your data and generate a homerange."))
	tkgrid(tklabel(dlg,text=""))
	tkgrid(tklabel(dlg,text="This procedure might take a few minutes."))
	 
	dups<-tclvalue(Dups.val)
	if(dups=='displace'){
		dups<-eval(parse(text=tclvalue(Displace.val)))
	}else{
		dups<-as.character(dups)
	}
	
	method<-tclvalue(Method.val)
	nVar<-eval(parse(text=tclvalue(Var.val)))
	if(method=='r'){
		aVal<-NULL
		kVal<-NULL
		rVal<-nVar
	}else if(method=='k'){
		aVal<-NULL
		rVal<-NULL
		kVal<-nVar
	}else if(method=='a'){
		rVal<-NULL
		kVal<-NULL
		aVal<-nVar
	}
	
	shpfile.name<-tclvalue(ShapeFile.val)
	txtfile.name<-tclvalue(TextFile.val)
	
	if(shpfile.name=='' && txtfile.name==''){
		xy<-eval(parse(text=tclvalue(XY.val)))
		
		IDs<-tclvalue(ID.val)
		if (IDs=='')
			IDs<-NULL
		else
			IDs<-eval(parse(text=IDs))
	}else{
		IDs<-NULL
		if(shpfile.name!=''){
			xy<-shapefile.points(shpfile.name)
		}else{
			xy<-read.delim(txtfile.name)[,1:2]
		}		
	}
	res<-NNCH(xy, id=IDs, k=kVal, r=rVal, a=aVal, status=TRUE, duplicates=dups)
	
	Iso.vect<-eval(parse(text=tclvalue(Isos.val)))
	
	if(as.character(tclvalue(GraphHR.val))=="1"){
		pdf(file=paste(tclvalue(Folder.val),'Homerange_Plot.pdf',sep='/'),width=10,height=10)
		plot(res,percent=Iso.vect)
		dev.off()
	}
	
	if(as.character(tclvalue(GraphAreaIso.val))=="1"){
		pdf(file=paste(tclvalue(Folder.val),'Area_vs_Iso_Plot.pdf',sep='/'))
		plot(NNCH.area(res,percent=Iso.vect))
		dev.off()
	}
	
	if(as.character(tclvalue(GraphkArea.val))=="1"){
		pdf(file=paste(tclvalue(Folder.val),'Area_vs_Variable(k|r|d).pdf',sep='/'))
		plot(NNCH.n.area(res,percent=Iso.vect,n=as.character(method)))
		dev.off()
	}
	
	if(as.character(tclvalue(ExportShape.val))=="1"){
		NNCH.export.shapefile(res,paste(tclvalue(Folder.val),'Isopleths',sep='/'),percent=Iso.vect)
	}
	
	tkconfigure(tt,cursor="arrow")
	tkdestroy(tt)
	tkdestroy(dlg)
	
	return(res)
}
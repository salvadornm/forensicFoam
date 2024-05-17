###############################################
# pseudo-python script to be used with Visit
#  SNM                          ------July 2020
# /Applications/VisIt.app/Contents/MacOS/VisIt -cli -nowin -s extract_images.py
############################################### 
# Step 1: Open the visit  database
OpenDatabase("/Volumes/LaCie/Chile_13July/DATA/solution.visit")
#OpenDatabase("/Users/snm/Data/a.visit")
#Step 2
# These are the variables' name the same way as they are in Visit
# For Instantaneous
variable=('U','V','W','F','DP','rho','P','Nfloor','Npart')
############################################## 
j='Nfloor'
min_plot  = 1
max_plot  = 100000
print "Variable : ",j
#--------------------------------------------------

## Annotation Options    ------------------------
b=AnnotationAttributes()
b.userInfoFlag=0
#b.axes2D.xAxis.title.title = "x"
#b.axes2D.yAxis.title.title = "y"
#b.axes2D.xAxis.label.font.scale = 2
#b.axes2D.yAxis.label.font.scale = 2
#b.axes2D.yAxis.title.font.font =  2
b.axes3D.visible = 0       # remove axes
b.axes3D.triadFlag = 0     # remove triad
b.axes3D.bboxFlag = 0      # remove bounding bpx
b.databaseInfoFlag = 0     # remove database info
b.legendInfoFlag = 0       # remove legend
SetAnnotationAttributes(b)

##  View point 3D ------------------------------
# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0, 0, +1)  # View Normal
View3DAtts.focus = (250, 0, 0)
View3DAtts.viewUp = (0, 1, 0)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 355.212
View3DAtts.nearPlane = -710.424
View3DAtts.farPlane = 710.424
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1.3 # Zoom 
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (250, 0, 0)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)


## Window -----------------
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "/Users/snm/Nfloor"
SaveWindowAtts.format = SaveWindowAtts.TIFF  
SaveWindowAtts.width = 1920
SaveWindowAtts.height = 1080
SetSaveWindowAttributes(SaveWindowAtts)


# DrawPlots()
# SaveWindow()
# exit()

#--------------------------------------------------
print "Printing All Times"
names = []
DeleteAllPlots()
for state in range(GetDatabaseNStates()):
      SetTimeSliderState(state)
      # plots---------------------------------------
      AddPlot("Pseudocolor",j)
      p = PseudocolorAttributes()
      p.colorTableName = "calewhite"
      p.min, p.minFlag=min_plot,1
      p.max, p.maxFlag=max_plot,1
      p.centering = p.Zonal #Natural,
      p.scaling = p.Log #Linear, Log, Skew
      SetPlotOptions(p)
      ## Threshold Nfloor > 1 ----------------------
      AddOperator("Threshold")
      ThresholdAtts = ThresholdAttributes()
      ThresholdAtts.listedVarNames = ("Nfloor")
      ThresholdAtts.zonePortions = (1)
      ThresholdAtts.lowerBounds = (0.1)
      ThresholdAtts.upperBounds = (100000)
      SetOperatorOptions(ThresholdAtts)
      ## Isosurface  FF=0.5  ------------------------
      AddOperator("Isosurface")
      IsosurfaceAtts = IsosurfaceAttributes()
      IsosurfaceAtts.contourValue = (0.5)
      IsosurfaceAtts.variable = "FF"
      SetOperatorOptions(IsosurfaceAtts)
      ### Output
      ## Window -----------------
#       SaveWindowAtts = SaveWindowAttributes()
#       SaveWindowAtts.outputToCurrentDirectory = 0
#       SaveWindowAtts.outputDirectory = "/Users/snm/Nfloor"
# SaveWindowAtts.format = SaveWindowAtts.TIFF  
# SaveWindowAtts.width = 1920
# SaveWindowAtts.height = 1080
# SetSaveWindowAttributes(SaveWindowAtts)


      ### -----------------------------
      DrawPlots()
      names = SaveWindow()
      DeleteAllPlots()
#--------------------------------------------------
print "End .. (press Ctrl + C)"
exit()

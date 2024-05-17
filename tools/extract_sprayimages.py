###############################################
# pseudo-python script to be used with Visit
#  SNM                          ------August 2020
# /Applications/VisIt.app/Contents/MacOS/VisIt -cli -nowin -s extract_sprayimages.py
############################################### 
### DATABASES
dbs = ("mesh_0.ply","VTK/a.visit")
# Step 0: Open solid mesh
#
OpenDatabase("mesh_0.ply")
# plot mesh quality/shape (not that important)
AddPlot("Pseudocolor","mesh_quality/shape")
meshplot = PseudocolorAttributes()
meshplot.colorTableName = "gray"
meshplot.lightingFlag = 0
# lighting off
SetPlotOptions(meshplot)

# Step 1: Open the visit  database
OpenDatabase(dbs[1])
# Step 2
# These are the variables' name the same way as they are in Visit
# For Instantaneous
variable=('velpart-x','velpart-y','velpart-z','n-parcels','mass','diameter','temppart','x-force','y-force','z-force','mass-term','heat-term')
############################################## 

j='diameter'
min_plot  = 1E-6
max_plot  = 5E-5
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
View3DAtts.viewNormal = (0.98, 0.1, -0.05)  # View Normal
View3DAtts.focus = (110, 4, 1)
View3DAtts.viewUp = (0, 0, 1)
View3DAtts.viewAngle = 100
View3DAtts.parallelScale = 500
View3DAtts.nearPlane = -989.57
View3DAtts.farPlane = 989.57
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 8 # Zoom 
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
#View3DAtts.centerOfRotation = (250, 0, 0)  # temp
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)

## Window -----------------
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "./sprayplots"
SaveWindowAtts.format = SaveWindowAtts.TIFF  
SaveWindowAtts.width = 3840
SaveWindowAtts.height = 2160
SetSaveWindowAttributes(SaveWindowAtts)


# black background TODO

#DrawPlots()
#SaveWindow()
#exit()

counter = 0
#--------------------------------------------------
print "Printing All Times"
names = []
#DeleteAllPlots()
for state in range(GetDatabaseNStates()):
      SetTimeSliderState(state)
      # plots---------------------------------------
      AddPlot("Pseudocolor",j)
      p = PseudocolorAttributes()
      p.min, p.minFlag=min_plot,1
      p.max, p.maxFlag=max_plot,1
      # geometry TODO
      # 
      SetPlotOptions(p)
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

      ##### free cache regularly 
      counter = counter + 1
      # every 10 states ?
      if (counter >10):
            ClearCache("localhost")
            counter = 0
#--------------------------------------------------
print "End .. (press Ctrl + C)"
exit()

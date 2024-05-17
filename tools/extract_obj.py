################################################################################
# pseudo-python script to be used with Visit
#  SNM                          ------May 2024
# execute with: visit -cli -nowin -s extract_obj.py
################################################################################ 
# Step 1: Open the database
OpenDatabase("./system/controlDict")
#Step 2
# These are the variables' names the same way as they appear in Visit
# (they may change often)
variable=('CO','CO2','H2O','F','DP','rho','P','Nfloor','Npart')
################################################################################
j='CO2'
min_plot  = 0
max_plot  = 1

varsmoke = "internalMesh/CO2"          # variable to indicate smoke
smokethreshold = 0.0002   # concentration that determines flame

print ("Variable Smoke: ",varsmoke)

#--------------------------------------------------
# The following liens are a text to show if reading/plotting/work
# they are not needed


## Annotation Options    ------------------------
b=AnnotationAttributes()
b.userInfoFlag=0
b.axes3D.visible = 0       # remove axes
b.axes3D.triadFlag = 0     # remove triad
b.axes3D.bboxFlag = 0      # remove bounding bpx
b.databaseInfoFlag = 0     # remove database info
b.legendInfoFlag = 0       # remove legend
SetAnnotationAttributes(b)

##  View point 3D ------------------------------
# Begin spontaneous state
# View3DAtts = View3DAttributes()
# View3DAtts.viewNormal = (0, 0, +1)  # View Normal
# View3DAtts.focus = (250, 0, 0)
# View3DAtts.viewUp = (0, 1, 0)
# View3DAtts.viewAngle = 30
# View3DAtts.parallelScale = 355.212
# View3DAtts.nearPlane = -710.424
# View3DAtts.farPlane = 710.424
# View3DAtts.imagePan = (0, 0)
# View3DAtts.imageZoom = 1.3 # Zoom 
# View3DAtts.perspective = 1
# View3DAtts.eyeAngle = 2
# View3DAtts.centerOfRotationSet = 0
# View3DAtts.centerOfRotation = (250, 0, 0)
# View3DAtts.axis3DScaleFlag = 0
# View3DAtts.axis3DScales = (1, 1, 1)
# View3DAtts.shear = (0, 0, 1)
# View3DAtts.windowValid = 1
# SetView3D(View3DAtts)


## Window -----------------
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = "./Test"
SaveWindowAtts.format = SaveWindowAtts.TIFF  
SaveWindowAtts.width = 1920
SaveWindowAtts.height = 1080
SetSaveWindowAttributes(SaveWindowAtts)


# uncoment the lines below to test reading/plotting
# DrawPlots()
# SaveWindow()
# print ("End Test .. (press Ctrl + C)")
# exit()

#--------------------------------------------------
print ("Printing All Times")
names = []
DeleteAllPlots()
for state in range(GetDatabaseNStates()):
      SetTimeSliderState(state)
      # plots---------------------------------------
      AddPlot("Pseudocolor",varsmoke)
      p = PseudocolorAttributes()
      #p.colorTableName = "calewhite"
      p.min, p.minFlag=min_plot,1
      p.max, p.maxFlag=max_plot,1
      p.centering = p.Zonal #Natural,
      p.scaling = p.Linear #Linear, Log, Skew
      p.renderPoints = 0
      p.smoothingLevel = 2
      SetPlotOptions(p)
      ## Isosurface   ------------------------
      AddOperator("Isosurface")
      IsosurfaceAtts = IsosurfaceAttributes()
      IsosurfaceAtts.contourValue = (smokethreshold)
      IsosurfaceAtts.contourMethod = IsosurfaceAtts.Value
      IsosurfaceAtts.variable = varsmoke
      SetOperatorOptions(IsosurfaceAtts)
      ### Output
      ## Window -----------------
      SaveWindowAtts = SaveWindowAttributes()
      SaveWindowAtts.outputToCurrentDirectory = 0
      SaveWindowAtts.outputDirectory = "./OUT"
      SaveWindowAtts.format = SaveWindowAtts.OBJ  
      SetSaveWindowAttributes(SaveWindowAtts)


      ### -----------------------------
      DrawPlots()
      names = SaveWindow()
      DeleteAllPlots()
#--------------------------------------------------
print ("End .. (press Ctrl + C)")
exit()

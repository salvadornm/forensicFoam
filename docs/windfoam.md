
## Grid

Terrain from STL file


The mesh is stored in **constant/polyMesh**

NOTE: The `./Allclean` script will delete the mesh, careful with it


## Patches and sets

Required for boundary conditions


## ABL


## Nudging Model


### Buffer Area

## Wind files (wind.dat)

It needs **wind.dat**,  four column file, with time, wind, wind direction and ground temperature.

```
#TIME   WIND(100 m)      WIND(degrees)  TGROUND
0       3.124100        140.194430      309.787196
1       3.929377        165.256440      308.007197
2       4.909175        176.496480      307.237196
3       4.925444        185.826340      305.457198
4       4.780168        195.780750      304.877198
```


## Parameter file

The **param.dat** contaisn parameter simulation

```
# parameter file for forensicWindFOAM
windsteps  20       
variable_wind 1     
currentwindstep  0  
nudging  1          
TIrmsbuf  0.05      
continuerun 0
```

### Problem definition



| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| **windsteps**                   | Int           |  0       | number of hours/step to read in the **wind.dat** file, by default hours                         |

## Run


Decompose the solution

```bash
$ decomposePar
```

Run in parallel

```bash
$ mpiexec -np 16 forensicWindFoam -parallel 
```

Reconstruct the solution

```bash
$ reconstructPar
```



### BC







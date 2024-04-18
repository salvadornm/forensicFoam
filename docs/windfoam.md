
## Grid

Terrain from STL file


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
windsteps  20       // number of wind hours read
variable_wind 1     // variable wind
currentwindstep  0  //  current wind step (for restarts)
nudging  1          // nudging model on/off
TIrmsbuf  0.05      // Turbulence intesnsity noise buffer
```










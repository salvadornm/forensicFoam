import math 
from pyproj import CRS
from pyproj import Transformer

#-------------------------------------
def xy_latlong(proj,lat,long,x0,y0):
    (x,y) = proj.transform(lat, long)
    xrel = x - x0
    yrel = y - y0
    return (xrel,yrel)
#-------------------------------------



# constants
sigma = 5.67e-8 # si units
Runi = 8.314 #J/mol K
kboltzmann =1.38066e-23 #J/K 
T0 = 273.15 #
NA = 6.023e23 # Avogadro
F = 96485.3 # C/mol
P0 = 103125  # 1atm
g =  9.81    # gravity
# unit conversion
unit_l2m = 0.001 
unit_cm2m = 0.01
unit_knts2ms = 0.514444

kappa = 0.41
Cmu  = 0.09

print('Lousiana numbers 2023')

print('  LOG WALL')

H=2000  # Height  1 km  !<------

#-----------------
#Uref = 3.12  # Knts
#Uref = 5*unit_knts2ms # reference velocity  (measured)
Uref = 3.12
zref = 100 # reference height
z0 = 0.5  # Surface roughness height 
zg =  0.0
zmax = 500
#-----------------

print(' [OpenFOAM] Uref for Uref=',Uref)

num = zref + z0
Ustar = kappa *Uref/math.log(num/z0)
print('U* (friction velocity)',Ustar, '[m/s]')


num = H - zg + z0
ln1 = math.log(num/z0)
print('U/U* ',ln1/kappa, " at z=",H)
U = Ustar*ln1/kappa
print('U  ',U, "[m/s] at z=",H)
print(' [OpenFOAM] number or top BC for 0/U=',U)

num = zmax - zg + z0
U = Ustar*math.log(num/z0)/kappa
print('U  ',U, "[m/s] at zmax=",zmax)


k = Ustar*Ustar/math.sqrt(Cmu)

print("k (turbulent kinetic energy) ",k, "[m2/s2] ")

print(' [OpenFOAM] number on inlet  for 0/k=',k)


# P=P0*exp(-z/z0) atmospheric drop
M = 0.028964 # kg/mol
T = 300.0
P = P0*math.exp(-g*M*H/(Runi*T))

zT = Runi*T/(g*M)

print("[OpenFOAM]   zT= ",zT," m")

print("P  (using chem pot)= ",P, "[Pa] at z= ",H," m")

rho0 = P0*M/(Runi*T)
rho  = P*M/(Runi*T)

print("rho0 = ",rho0, "[Pa] at z= ",0," m")
print("rho = ",rho, "[Pa] at z= ",H," m")


print(" press change=",100*(P-P0)/P0," %")
print(" rho change=",100*(rho-rho0)/rho0," %")

rho1 = rho0
P = P0 - rho1*g*H

print("P  (hydro)= ",P, "[Pa] at z= ",H," m")

print("---------------------------------------------------")
print(" X,Y coordinates of points of interest ")
print("---------------------------------------------------")

crs = CRS.from_epsg(3857)
crs.geodetic_crs
proj = Transformer.from_crs(crs.geodetic_crs, crs)

# ORIGIN 0,0 of simulation 30.038, -90.643
lat0= 30.038
long0=-90.643
(x0,y0) = proj.transform(lat0, long0)

# list of data from CTH
LOI = [(30.0616417,-90.582949) , (30.0633719, -90.583076) ]

print("------------------------------------")
print(" List of CTEH ")
print("------------------------------------")

for loc in LOI:
    (xrel,yrel) = xy_latlong(proj,loc[0],loc[1],x0,y0)
    print(f"{xrel:.5}",f"{yrel:.5}")

# TANKS
TANK = [(30.059,-90.595) , (30.058, -90.594) ]

print("------------------------------------")
print(" List of TANK locations ")
print("------------------------------------")

for loc in TANK:
    (xrel,yrel) = xy_latlong(proj,loc[0],loc[1],x0,y0)
    print(f"{xrel:.5}",f"{yrel:.5}")


# TANK 
Dtank1 = 70.5

Areatank1 = math.pi*70*70*0.25

nbarrels = 70
barrel = 159 #litre

volnafta = nbarrels*barrel*unit_l2m

print(" Vol=",volnafta, " [m3]")

time = 2 * 3600 # 2 h

vdot = volnafta/time

print(" vdot=",vdot, " [m3/2]")

vdotair = vdot*600*15

print(" vdotair=",vdotair, " [m3/2]")

print(" ---------------------------------------------------")












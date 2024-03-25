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
unit_lbs2kg = 0.453592
unit_h2sec = 3600

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

print("------------------------------------")

# TANK 
Dtank1 = 70.5  #m

Areatank1 = math.pi*70*70*0.25  #m2

nbarrels = 70
barrel = 159 #litre

volnafta = nbarrels*barrel*unit_l2m

#print(" Vol=",volnafta, " [m3]")

time = 2 * 3600 # 2 h

vdot = volnafta/time

#print(" vdot=",vdot, " [m3/2]")

vdotair = vdot*600*15

#print(" vdotair=",vdotair, " [m3/2]")


# EPA report
# flammable mixture taken as propane C3H8
mass_mixture = 6000000 #lbs
mass_mixture = mass_mixture*unit_lbs2kg

time_burn = 5 #h  from 8:30 to 13:30
time_burn = time_burn*unit_h2sec

mdot_mixture = mass_mixture/time_burn

print(" mass mixture burn", mass_mixture, " [kg] in t=",time_burn, " [secs]")
print(" mass flow mixture burn EPA report", mdot_mixture," [kg/s]")


# C3H8 + 5 (O2 + 3.72 N2) => 3 CO2 + 4H2O +  5*3.72*N2
# 1 mol propane reacts 5 mols air, 3 mols CO2, 4 mols H2O ,18.6 mols N2 
# 
n_air = 5
n_co2 = 3
n_h2o = 4
n_n2  = 5*3.72
ntot = n_co2 + n_h2o + n_n2

Mpropane = 44.094 # g/mol  or kg/kmol
M_co2    = 44 
M_h2o    = 18
M_n2     = 28
M_air    = 29


print(" assume mixture is propane mostly")

ndot_propane = mdot_mixture/Mpropane  # kmol

print(" ndot C3H8=",ndot_propane," [kmol/s]")


ndot_co2 = n_co2*ndot_propane
ndot_h2o = n_h2o*ndot_propane
ndot_n2  = n_n2*ndot_propane

print(" ndot CO2 H2O N2")
print(ndot_co2,ndot_h2o,ndot_n2," [kmol/s]")

X_co2 = n_co2/ntot
X_n2  = n_n2/ntot
X_h2o = n_h2o/ntot

MMW = X_co2*M_co2 + X_h2o*M_h2o + X_n2*M_n2
print(" X  CO2 H2O N2    MOLE FRACTION")
print(X_co2,X_h2o,X_n2," [-]")
Y_co2 = X_co2*M_co2/MMW
Y_n2  = X_n2*M_n2/MMW
Y_h2o = X_h2o*M_h2o/MMW

Y_o2 = 0  # all O2 consummed

print(" Y  CO2 H2O N2    MASS FRACTION")
print(Y_co2,Y_h2o,Y_n2," [-]")

T_unburn = 300
T_burn   = 1900 *0.8 # not adiabatic  (large heat losses)

rho_unburn = 1.808   # kg/m3 density
rho_burn   = rho_unburn*T_unburn/T_burn

print(" rho burn (propane)", rho_burn," [kg/m3]")

ndot_air = 5*ndot_propane
mdot_air = ndot_air*M_air

print(" mdot propane ", mdot_mixture," mdot_air",mdot_air, " [kg/s]")

mdot_mixture = mdot_mixture + mdot_air
print(" TOTAL MASS = ",mdot_mixture, " [kg/s]")

rhoU = mdot_mixture/Areatank1  # kg/m2 s

Vtank1 = rhoU/rho_burn

print(" Tank 1 velocity = ",Vtank1, " [m/s]")

print( " ---------------------- ")
print( " BC to set in OpenFOAM ")
print( " patch tank1 ")
print( " 0/U ",Vtank1)
print( " 0/T ",T_burn)
print( " 0/YCO2 ",Y_co2)
print( " 0/YN2  ",Y_n2)
print( " 0/YO2  ",Y_o2)
print( " 0/YH2O ",Y_h2o)
print( " ---------------------- ")



print(" ---------------------------------------------------")












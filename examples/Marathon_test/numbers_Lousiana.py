import math 
from pyproj import CRS
from pyproj import Transformer
from array import array

#-------------------------------------
def xy_latlong(proj,lat,long,x0,y0):
    (x,y) = proj.transform(lat, long)
    xrel = x - x0
    yrel = y - y0
    return (xrel,yrel)
#-------------------------------------

# Compute Mean molecular weight
#-------------------------------------
def comp_MMW(x,mw,nspecies):
    MMW = 0.0
    for i in range(nspecies):
        MMW = MMW + x[i]*mw[i]
    return MMW
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
unit_l2m = 0.001   # litre to m3
unit_cm2m = 0.01
unit_knts2ms = 0.514444
unit_lbs2kg = 0.453592
unit_h2sec = 3600
unit_barrels2lbs = 300

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

# TANK 1 (150-11)  Naphta
Dtank1 = 70.5  #m
# TANK 2 (300-7)   Low Sulfur Diesel
Dtank2 = 46.0  # m
# TANK 3 (200-7)   Gas-oil
Dtank3 = 52.0 #m

# AREAS
Areatank1 = math.pi*Dtank1*Dtank1*0.25  #m2
Areatank2 = math.pi*Dtank2*Dtank2*0.25  #m2
Areatank3 = math.pi*Dtank3*Dtank3*0.25  #m2

Areacont = 17000 # m2 (estimate from STL)

# quantities
barrels_tank1 = 22013 
barrel_tank2 = 40000
barrels_leaked = 70

leakrate = 8.75 # leak rate estimated per h



barrels_tank1_0 = barrels_tank1 - barrels_leaked

barrels_cont = barrels_leaked+17*leakrate

#volnafta = nbarrels*barrel*unit_l2m



# EPA report

time_release =  17 *3600  #  17 hfrom 6:30 to 23:30 

barrels_tank1_remain = barrels_tank1_0 - time_release*leakrate/3600


print(" Initial barrels", barrels_tank1 )
print(" Initial barrels 6:30 ",barrels_tank1_0)
print(" Initial barrels 8:30 ",barrels_tank1_0-2*leakrate)
print(" Barrels end  (BURN in tank)",barrels_tank1_remain)
print(" Barrels      (BURN in containment)",barrels_cont)



COreleased = 57496 #lbs during the event
NOreleased = 10567 
SO2released = 5536
Bzreleased = 1433
Xyreleased = 31530
COreleased = COreleased*unit_lbs2kg #kg
NOreleased = NOreleased*unit_lbs2kg #kg
SO2released = SO2released*unit_lbs2kg #kg
Bzreleased = Bzreleased*unit_lbs2kg #kg
Xyreleased = Xyreleased*unit_lbs2kg #kg

COrate_avg = COreleased/time_release
NOrate_avg = NOreleased/time_release
SO2rate_avg = SO2released/time_release
Bzrate_avg = Bzreleased/time_release
Xyrate_avg = Xyreleased/time_release


print(" In time =",time_release/3600, " [hours]")
print(" Average massrate CO (EPA) ",COrate_avg, "[kg/s]")
print(" Average massrate NO (EPA) ",NOrate_avg, "[kg/s]")
print(" Average massrate SO2 (EPA) ",SO2rate_avg, "[kg/s]")
print(" Average massrate Bz (EPA) ",Bzrate_avg, "[kg/s]")
print(" Average massrate Xy (EPA) ",Xyrate_avg, "[kg/s]")

t1 = 2*3600  # 2 h   6:30--> 8:30   Growth

t2 = 5*3600  # 5 h   8:30-->13:30   Qmax 

t3 = 17*3600 # 17 h   Decay

aux = (t1/3 + (t2 - t1) + (t3 - t2)/3)*t1**2

print(" aux ",aux, "[kg/s]")
alfa  = COreleased/aux
alfa2 = alfa*t1**2/(t3 - t2)**2 
print(" alfa CO ",alfa, " alfa2",alfa2)
COrate_max = alfa*t1**2
rmax = COrate_max/COrate_avg
print(" ratiop max/average = ",rmax)

print("COrate_max ",COrate_avg*rmax, " COrate_avg",COrate_avg)
print("NOrate_max ",NOrate_avg*rmax, " NOrate_avg",NOrate_avg)
print("SO2rate_max ",SO2rate_avg*rmax, " SO2rate_avg",SO2rate_avg)
print("Bzrate_max ",Bzrate_avg*rmax, " Bzrate_avg",Bzrate_avg)
print("Xyrate_max ",Xyrate_avg*rmax, " Xyrate_avg",Xyrate_avg)

print("------------------------------------")


T_unburn = 300
Tatm = 308.0  # atm

print("------------------------------------")




Mc = 12.011
Mh = 1.008
Mo = 15.999
Mn = 14.0067 

M_co2    = Mc + 2*Mo
M_h2o    = 2*Mh + Mo
M_n2     = 2*Mn
M_o2     = 2*Mo
M_air    = 0.21*M_o2 + 0.79*M_n2 
M_c7h16  = 7*Mc+ 16*Mh


print(" some molecular weights..") 
print(" M_co2 M_h2o M_n2  : ",M_co2,M_h2o,M_n2)

names_spec = ["O2","N2","CO2","H2O"]
iO2  = 0
iN2  = 1
iCO2 = 2
iH2O = 3
nspecies = 4
mol_weight = array('f',[0.0,0.0,0.0,0.0]) 

mol_weight[iO2] = M_o2
mol_weight[iN2] = M_n2
mol_weight[iCO2] = M_co2
mol_weight[iH2O] = M_h2o

# temporary arrays
n_mol = array('f',[0.0,0.0,0.0,0.0]) 

names_emiss = ["CO","NO","SO2","Benzene","Xylene"]

iCO  = 0
iNO  = 1
iSO2 = 2
iBz  = 3
iXy  = 4

nemissions = 5
mol_weight_e = array('f',[0.0,0.0,0.0,0.0,0.0]) 

mol_weight_e[iCO]  = 28.01
mol_weight_e[iNO]  = 30.01 
mol_weight_e[iSO2] = 64.066
mol_weight_e[iBz]  = 78.11
mol_weight_e[iXy]  = 106.16

release_e = array('f',[0.0,0.0,0.0,0.0,0.0]) 
y_e = array('f',[0.0,0.0,0.0,0.0,0.0]) 


release_e[iCO]  = 57496 #lbs during the event
release_e[iNO]  = 10567 
release_e[iSO2] = 5536
release_e[iBz]  = 1433
release_e[iXy]  = 31530

# change to kg
for i in range(nemissions):
   release_e[i] = release_e[i]*unit_lbs2kg

print("------------------------------------")
print("  TANK 1 (Naphta)")
print("------------------------------------")
x_tank1 = array('f',[0.0,0.0,0.0,0.0]) 
y_tank1 = array('f',[0.0,0.0,0.0,0.0]) 

# CnHm
# CnHn + (n + m/4)  (O2 + 3.76 N2) => n CO2 + m/2 H2O +  (n + m/4) *3.76*N2
# 1 mol alkene reacts n+m/4 mols air, n mols CO2, m/2 mols H2O , (n + m/4) *3.76  mols N2 


# Closest alkene to Naphta is propene (C3H6)
n=3
m=6

M_naphta = n*Mc+ m*Mh

print(" M_naphta    :",M_naphta)

n_air = n + m/4
n_co2 = n
n_h2o = m/2
n_n2  = n_air*79/21
n_o2  = 0
ntot = n_co2 + n_h2o + n_n2

n_mol[iO2]  = n_o2
n_mol[iN2]  = n_n2
n_mol[iH2O] = n_h2o
n_mol[iCO2] = n_co2

# mole fraction and mean molecular weight
for i in range(nspecies):
   x_tank1[i] = n_mol[i]/ntot

MMW = comp_MMW(x_tank1,mol_weight,nspecies)

# mass fraction
for i in range(nspecies):
   y_tank1[i] = x_tank1[i]*mol_weight[i]/MMW

# plot
for i in range(nspecies):
   print(names_spec[i]," = ",x_tank1[i],y_tank1[i])    


# Temperature and Density
#  2505 K is adiabatic T at Stoichiometric
#  images show orange flames appox 900 C 
Ttank1 = 900+ T0  


rho_unburn = P0*M_naphta/Runi/Tatm/1000

print(" rho Naphta unburn  ",rho_unburn," [kg/m3]")

rho_burn   = P0*MMW/Runi/Ttank1/1000

print(" rho burn (combustion gases)", rho_burn," [kg/m3]")


mass_tank1 = barrels_tank1_remain   # mass burned in barrels
mass_tank1 = mass_tank1*unit_barrels2lbs  # in[lbs]

print(" mass in lbs = ",mass_tank1, " in kgs ",mass_tank1*unit_lbs2kg)
mass_tank1 = mass_tank1*unit_lbs2kg  # in [kgs]

# mass_tank1 is the quantity of naphta burnt in the tank 1 over 15 hours
#

time_burn = 15 #h  from 8:30 to 23:30
time_burn = time_burn*unit_h2sec

mdot_tank1 = mass_tank1/time_burn
print(" mdot [TANK1 average] = ",mdot_tank1, " [kg/s] ")

ndot_tank1 = mdot_tank1/M_naphta
ndot_air = n_air*ndot_tank1
mdot_air = ndot_air*M_air


print(" mdot tank1 ", mdot_tank1," mdot_air",mdot_air, " [kg/s]")

mdot_mixture = mdot_tank1 + mdot_air
print(" TOTAL MASS = ",mdot_mixture, " [kg/s]")

massprod_tank1 = mdot_mixture*time_burn

rhoU = mdot_mixture/Areatank1  # kg/m2 s

Vtank1= rhoU/rho_burn

print(" Tank 1 (150-11) velocity  = ",Vtank1, " [m/s]")
print("                 velocity MAX = ",Vtank1*rmax, " [m/s]")



print("------------------------------------  ")
print("  Contaiment Area (Naphta in the soil)")
print("   **  burns lean  at LFL ")
print("------------------------------------  ")
x_cont = array('f',[0.0,0.0,0.0,0.0]) 
y_cont = array('f',[0.0,0.0,0.0,0.0]) 


# The contaiment area burns lean
# propene LFL is 2.75 %

LFL = 0.0275 # npropene/nair

#print(" original fuel/air ", 1/n_air," or ",1/n_air*100," %")

# CnHn + (n + m/4)  (O2 + 3.76 N2) + excess (O2 + 3.76 N2)
#   => n CO2 + m/2 H2O +  (n + m/4) *3.76*N2 + excess  (O2 + 3.76 N2)

n_air = 1/LFL
n_excess = n_air - n + m/4
n_co2 = n
n_h2o = m/2
n_n2  = n_air*79/21 
n_o2  = n_excess
ntot = n_co2 + n_h2o + n_n2 + n_o2

#print(" n_air ",n_air," (stoichio) ", n + m/4)


n_mol[iO2]  = n_o2
n_mol[iN2]  = n_n2
n_mol[iH2O] = n_h2o
n_mol[iCO2] = n_co2

# mole fraction and mean molecular weight
for i in range(nspecies):
   x_cont[i] = n_mol[i]/ntot

MMW = comp_MMW(x_cont,mol_weight,nspecies)

# mass fraction
for i in range(nspecies):
   y_cont[i] = x_cont[i]*mol_weight[i]/MMW

# plot
for i in range(nspecies):
   print(names_spec[i]," = ",x_cont[i],y_cont[i])    


# Temperature and Density

Tcont = 827 + T0  # 827 C experimental measuremnt pool fires

rho_unburn = P0*M_naphta/Runi/Tatm/1000

print(" rho Naphta unburn  ",rho_unburn," [kg/m3]")

rho_burn   = P0*MMW/Runi/Tcont/1000

print(" rho burn (combustion gases)", rho_burn," [kg/m3]")


mass_cont = barrels_cont   # mass burned in barrels
mass_cont = mass_cont*unit_barrels2lbs  # in[lbs]

print(" mass in lbs = ",mass_cont, " in kgs ",mass_cont*unit_lbs2kg)
mass_cont= mass_cont*unit_lbs2kg  # in [kgs]

# mass_cont is the quantity of naphta burnt in the cont area over 15 hours
#

time_burn = 17 #h  from 6:30 to 23:30
time_burn = time_burn*unit_h2sec

mdot_cont= mass_cont/time_burn
print(" mdot [CONT average] = ",mdot_cont, " [kg/s] ")

ndot_cont = mdot_cont/M_naphta
ndot_air = n_air*ndot_tank1
mdot_air = ndot_air*M_air

print(" mdot cont ", mdot_cont," mdot_air",mdot_air, " [kg/s]")

mdot_mixture = mdot_cont + mdot_air
print(" TOTAL MASS = ",mdot_mixture, " [kg/s]")

massprod_cont = mdot_mixture*time_burn

rhoU = mdot_mixture/Areacont  # kg/m2 s

Vcont = rhoU/rho_burn

print(" Containment velocity = ",Vcont, " [m/s]")
print("                 velocity MAX = ",Vcont*rmax, " [m/s]")



# pollutants



# emissions are released by Tank1 and containment

for i in range(nemissions):
   y_e[i] = release_e[i]/(massprod_tank1+massprod_cont)

print(" MASS FRACTIONS ******")

for i in range(nemissions):
   print(names_emiss[i]," = ",y_e[i]) 

print(" ******")



print("------------------------------------")
print("  TANK 2 (Diesel)")
print("------------------------------------")
x_tank2 = array('f',[0.0,0.0,0.0,0.0]) 
y_tank2 = array('f',[0.0,0.0,0.0,0.0]) 

# Diesel burn as n-heptane

M_diesel = M_c7h16
print(" M_diesel    :",M_diesel)


# C7H16 + 11 (O2 + 3.76 N2) => 7 CO2 + 8H2O +  11*3.76*N2
# 1 mol n-heptane reacts 11 mols air, producing 7 mols CO2, 8 mols H2O , 41.38mols N2 
n_air = 11
n_co2 = 7
n_h2o = 8
n_n2  = 11*79/21
n_o2 = 0
ntot = n_co2 + n_h2o + n_n2 + n_o2

n_mol[iO2]  = n_o2
n_mol[iN2]  = n_n2
n_mol[iH2O] = n_h2o
n_mol[iCO2] = n_co2

# mole fraction and mean molecular weight
for i in range(nspecies):
   x_tank2[i] = n_mol[i]/ntot

MMW = comp_MMW(x_tank2,mol_weight,nspecies)

# mass fraction
for i in range(nspecies):
   y_tank2[i] = x_tank2[i]*mol_weight[i]/MMW

# plot
for i in range(nspecies):
   print(names_spec[i]," = ",x_tank2[i],y_tank2[i])    

# Temperature and Density  TANK2

#  2469 K is adaibatic T at Stoichiometric
#  images show orange flames appox 900 C 
Ttank2 = 900 + T0 

rho_unburn = P0*M_diesel/Runi/Tatm/1000

print(" rho Diesel unburn  ",rho_unburn," [kg/m3]")

rho_burn   = P0*MMW/Runi/Ttank2/1000

print(" rho burn (combustion gases)", rho_burn," [kg/m3]")


mass_tank2 = barrel_tank2   # mass burned in barrels
mass_tank2 = mass_tank2*unit_barrels2lbs  # in[lbs]

print(" mass in lbs = ",mass_tank2, " in kgs ",mass_tank1*unit_lbs2kg)
mass_tank2 = mass_tank2*unit_lbs2kg  # in [kgs]

# mass_tank2 is the quantity of diesel burnt in the tank 1 over 15 hours
#

time_burn = 15 #h  from 8:30 to 23:30
time_burn = time_burn*unit_h2sec

mdot_tank2 = mass_tank2/time_burn
print(" mdot [TANK2 average] = ",mdot_tank2, " [kg/s] ")

ndot_tank2 = mdot_tank2/M_naphta
ndot_air = n_air*ndot_tank2
mdot_air = ndot_air*M_air

print(" mdot tank2 ", mdot_tank2," mdot_air",mdot_air, " [kg/s]")

mdot_mixture = mdot_tank2 + mdot_air
print(" TOTAL MASS = ",mdot_mixture, " [kg/s]")

rhoU = mdot_mixture/Areatank2  # kg/m2 s

Vtank2= rhoU/rho_burn

print(" Tank 2 (300-7) velocity  = ",Vtank2, " [m/s]")
print("                 velocity MAX = ",Vtank2*rmax, " [m/s]")

print("------------------------------------")
print("  TANK 3 (gasoil)")
print("------------------------------------")
x_tank3 = array('f',[0.0,0.0,0.0,0.0]) 
y_tank3 = array('f',[0.0,0.0,0.0,0.0]) 

# CnHm
# CnHn + (n + m/4)  (O2 + 3.76 N2) => n CO2 + m/2 H2O +  (n + m/4) *3.76*N2
# 1 mol alkene reacts n+m/4 mols air, n mols CO2, m/2 mols H2O , (n + m/4) *3.76  mols N2 
# Closest alkene to gasoil is C12H23
n=12
m=23

M_gasoil = n*Mc+ m*Mh

print(" M_gasoil    :",M_gasoil)

n_air = n + m/4
n_co2 = n
n_h2o = m/2
n_n2  = n_air*79/21
n_o2  = 0
ntot = n_co2 + n_h2o + n_n2

n_mol[iO2]  = n_o2
n_mol[iN2]  = n_n2
n_mol[iH2O] = n_h2o
n_mol[iCO2] = n_co2

# mole fraction and mean molecular weight
for i in range(nspecies):
   x_tank3[i] = n_mol[i]/ntot

MMW = comp_MMW(x_tank3,mol_weight,nspecies)

# mass fraction
for i in range(nspecies):
   y_tank3[i] = x_tank3[i]*mol_weight[i]/MMW

# plot
for i in range(nspecies):
   print(names_spec[i]," = ",x_tank3[i],y_tank3[i])    

## Temperature and Density  TANK3

#  2469 K is adaibatic T at Stoichiometric
#  images show orange flames appox 900 C 
Ttank3 = 900 + T0 

rho_unburn = P0*M_gasoil/Runi/Tatm/1000

print(" rho Gasoil unburn  ",rho_unburn," [kg/m3]")

rho_burn   = P0*MMW/Runi/Ttank3/1000

print(" rho burn (combustion gases)", rho_burn," [kg/m3]")

mass_tank3 = 2300000  # lbs  <-------------------
mass_tank3 = mass_tank3*unit_lbs2kg   # mass burned in kg

print(" mass tank3 = ",mass_tank3, " in kg ")

# mass_tank3 is the quantity of gasoil in tank 3 
# burn over 3 hours

time_burn = 3 #h  from 1:30 (approx) to 4:30
time_burn = time_burn*unit_h2sec

mdot_tank3 = mass_tank3/time_burn

print(" mdot [TANK3 average] = ",mdot_tank3, " [kg/s] ")

ndot_tank3 = mdot_tank3/M_gasoil
ndot_air = n_air*ndot_tank3
mdot_air = ndot_air*M_air

print(" mdot tank3 ", mdot_tank3," mdot_air",mdot_air, " [kg/s]")

mdot_mixture = mdot_tank3 + mdot_air
print(" TOTAL MASS = ",mdot_mixture, " [kg/s]")


rhoU = mdot_mixture/Areatank3  # kg/m2 s

Vtank3= rhoU/rho_burn

print(" Tank 3 (200-7) velocity  = ",Vtank3, " [m/s]")
print("                velocity MAX = ",Vtank3*rmax, " [m/s]")

# pollutants

print("------------------------------------------")


print( " ---------------------- ")
print( " BC to set in OpenFOAM ")
print( " patch tank1 (TANK 150-11)")
print( " 0/U ",Vtank1)
print( " 0/T ",Ttank1)
print( " 0/YCO2 ",y_tank1[iCO2])
print( " 0/YN2  ",y_tank1[iN2])
print( " 0/YO2  ",y_tank1[iO2])
print( " 0/YH2O ",y_tank1[iH2O])
print( " ---------------------- ")
print( " patch cont (contaiment area)")
print( " 0/U ",Vcont)
print( " 0/T ",Tcont)
print( " 0/YCO2 ",y_cont[iCO2])
print( " 0/YN2  ",y_cont[iN2])
print( " 0/YO2  ",y_cont[iO2])
print( " 0/YH2O ",y_cont[iH2O])
print( " ---------------------- ")
print( " patch tank2 (TANK 300-7)")
print( " 0/U ",Vtank2)
print( " 0/T ",Ttank2)
print( " 0/YCO2 ",y_tank2[iCO2])
print( " 0/YN2  ",y_tank2[iN2])
print( " 0/YO2  ",y_tank2[iO2])
print( " 0/YH2O ",y_tank2[iH2O])
print( " ---------------------- ")
print( " patch tank3 (TANK 200-7)")
print( " 0/U ",Vtank3)
print( " 0/T ",Ttank3)
print( " 0/YCO2 ",y_tank3[iCO2])
print( " 0/YN2  ",y_tank3[iN2])
print( " 0/YO2  ",y_tank3[iO2])
print( " 0/YH2O ",y_tank3[iH2O])
print( " ---------------------- ")

print( " 0/YCO ",y_e[iCO])
print( " 0/YNO  ",y_e[iNO])
print( " 0/YSO2  ",y_e[iSO2])
print( " 0/YBz ",y_e[iBz])
print( " 0/YXy ",y_e[iXy])

print( " ---------------------- ")






print(" ---------------------------------------------------")












import math 

# constants
sigma = 5.67e-8 # si units
Runi = 8.314 #J/mol K
kboltzmann =1.38066e-23 #J/K 
cm2m = 100.0#
T0 = 273.15 #
NA = 6.023e23 # Avogadro
F = 96485.3 # C/mol
P0 = 103125  # 1atm
g =  9.81    # garvity

kappa = 0.41
Cmu  = 0.09

print('Lousiana numbers 2023')

print('  LOG WALL')

H=2000  # Height  1 km

Knts2ms = 0.514444 # 1 knts = 0.51444 m/s

#-----------------
Uref = 5  # Knts
Uref = 5*Knts2ms # reference velocity  (measured)
zref = 9.1 # reference height
z0 = 0.5  # Surface roughness height 
zg =  0.0
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




print(" ---------------------------------------------------")












using PyCall
using Plots
@pyimport scipy.optimize as opt
@pyimport scipy.integrate as int
@pyimport CoolProp.CoolProp as CP

#Define known variables
T_inf = 300 #Air temperature in kelvin
p = 1 #Air pressure in atm
p = p*101325 #convert air pressure to Pa
rh_inf = 0 #Relative Humidity of air
v = 10 #Air speed in m/s
a = 1 #Surface area of object in crossflow, m^2
l = 1 #Surface characteristic lenth, m
D_ab = 0.26E-4 #Water-Air Diffusivity, m^2/s
M_a = 18.015 #Molar mass of water, kg/kmol
h_fg = 2256 #Enthalpy of vaporization of water, kj/kg
R = 8.3143 #Universal gas constant, kj/kmol*K
m_0 = 1 #initial water mass, kg
T_0 = 275 #initial sphere temperature


#A) Calculate the steady-state temperature for the water surface.


#Find relevant properties
rho_air = CP.PropsSI("D","T",T_inf,"P",p,"Air") #Compute the air density (kg/m^3)
cp_air = CP.HAPropsSI("cp","T",T_inf,"P",p,"R",rh_inf) #Compute the specific heat capacity of air (J/kg dry air/K)
cp_air = cp_air/1000 #convert to kJ/kg-K
k_air = CP.PropsSI("conductivity","T",T_inf,"P",p,"Air") #Compute the air thermal conductivity (W/m/K)
alpha = k_air/(rho_air*cp_air) #Compute the thermal diffusivity (m^2/s)
Le = alpha/D_ab #compute lewis number
#Steady state temperature is reqched when heat transfer in from convection equals heat transfer out from evaporative cooling
B = (M_a*h_fg*p)/(R*rho_air*cp_air*Le^2/3) #Compute the coefficient B (K^2) (see Incropera example 6.8)
T_ss = (T_inf + sqrt((T_inf^2)-(4*B)))/2 #Compute the steady state temperature of the sphere
print(T_ss)

viscosity = CP.HAPropsSI("Visc","T",T_inf,"P",p,"R",rh_inf) #Compute the dynamic viscosity of air (Pa-s)
Re = v*l/viscosity #Calculate the reynolds number of air
Pr = viscosity/alpha #Calculate the prandtl number
Nu = 0.43*(Re^0.58)*(Pr^0.4) #Calculate the nusselt number, this is geometry dependent

#Compute the heat transfer coefficient
rho_water = CP.PropsSI("D","T",T_inf,"P",p_inf,"Water") #Density of water
cp_water =  CP.PropsSI("CPMASS","T",T_inf,"P",p_inf,"Water") #CP of water
k = alpha*rho_water*cp_water   #Thermal conductivity of the object W/m-K
h = k*Nu/l #Heat transfer coefficient (W/m^2-K)

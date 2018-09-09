from CoolProp.CoolProp import HAPropsSI
import CoolProp

#Define known variables
t_inf = 300 #Air temperature in kelvin
p_inf = 1 #Air pressure in atm
p_inf = p_inf*101325 #convert air pressure to Pa
rh_inf = 0 #Relative Humidity of air
v_inf = 10 #Air speed in m/s
a = 1 #Surface area of object in crossflow, m^2
l = 1 #Surface characteristic lenth
alpha = 0.26E-4 #Water Air Diffusivity, m^2/s

#A) Calculate the steady-state temperature for the water surface.

#Find relevant properties
viscosity = HAPropsSI('Visc','T',t_inf,'P',p_inf,'R',rh_inf) #Compute the dynamic viscosity
Re = v_inf*l/viscosity #Calculate the reynolds number
Pr = viscosity/alpha #Calculate the prandtl number
Nu = 0.43*(Re**0.58)*(Pr**0.4) #Calculate the nusselt number, this is geometry dependent

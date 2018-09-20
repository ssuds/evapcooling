import CoolProp.CoolProp as CP

#Define known variables
t_inf = 300 #Air temperature in kelvin
p_inf = 1 #Air pressure in atm
p_inf = p_inf*101325 #convert air pressure to Pa
rh_inf = 0 #Relative Humidity of air
v_inf = 10 #Air speed in m/s
a = 1 #Surface area of object in crossflow, m^2
l = 1 #Surface characteristic lenth, m
alpha = 0.26E-4 #Water-Air Diffusivity, m^2/s

#A) Calculate the steady-state temperature for the water surface.

#Find relevant properties
viscosity = CP.HAPropsSI('Visc','T',t_inf,'P',p_inf,'R',rh_inf) #Compute the dynamic viscosity of air (Pa-s)
cp_air = CP.HAPropsSI('cp','T',t_inf,'P',p_inf,'R',rh_inf) #Compute the specific heat capacity of air (J/kg dry air/K)
rho_air = CP.PropsSI('D','T',t_inf,'P',p_inf,'Air') #Compute the air density (kg/m^3)
Re = v_inf*l/viscosity #Calculate the reynolds number of air
Pr = viscosity/alpha #Calculate the prandtl number
Nu = 0.43*(Re**0.58)*(Pr**0.4) #Calculate the nusselt number, this is geometry dependent

#Compute the heat transfer coefficient
rho_water = CP.PropsSI('D','T',t_inf,'P',p_inf,'Water') #Density of water
cp_water =  CP.PropsSI('CPMASS','T',t_inf,'P',p_inf,'Water') #CP of water
k = alpha*rho_water*cp_water   #Thermal conductivity of the object W/m-K 
k2 = CP.PropsSI('conductivity','T',t_inf,'P',p_inf,'Water') #Thermal conductivity of the object W/m-K 
h = k*Nu/l #Heat transfer coefficient (W/m^2-K)

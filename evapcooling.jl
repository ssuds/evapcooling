using PyCall
using Plots
using Statistics
@pyimport scipy.optimize as opt
@pyimport scipy.integrate as int
@pyimport CoolProp.CoolProp as CP

#Define known variables
T_inf = 300 #Air temperature in kelvin
p = 1 #Air pressure in atm
p = p*101325 #convert air pressure to Pa
rh_inf = 0 #Relative Humidity of air
v = 10 #Air speed in m/s
a_s = 1 #Surface area of object in crossflow, m^2
l = 1 #Surface characteristic lenth, m
D_ab = 0.26E-4 #Water-Air Diffusivity, m^2/s
M_a = 18.015 #Molar mass of water, kg/kmol
h_fg = 2256 #Enthalpy of vaporization of water, kj/kg
R = 8.3143 #Universal gas constant, kj/kmol*K
m_0 = 1 #initial water mass, kg #Assumption for part B
T_0 = 275 #initial sphere temperature, K  #Assumption for part B


#A) Calculate the steady-state temperature for the water surface.


#Find relevant properties
rho_air = CP.PropsSI("D","T",T_0,"Q",1,"Water") #Compute the air density (kg/m^3)
cp_air = CP.HAPropsSI("cp","T",T_0,"P",p,"R",rh_inf) #Compute the specific heat capacity of air (J/kg dry air/K)
cp_air = cp_air/1000 #convert to kJ/kg-K
k_air = CP.PropsSI("conductivity","T",T_0,"P",p,"Air") #Compute the air thermal conductivity (W/m/K)
alpha = k_air/(rho_air*cp_air) #Compute the thermal diffusivity (m^2/s)
Le = alpha/D_ab #compute lewis number
#Steady state temperature is reqched when heat transfer in from convection equals heat transfer out from evaporative cooling
B = (M_a*h_fg*p)/(R*rho_air*cp_air*Le^2/3) #Compute the coefficient B (K^2) (see Incropera example 6.8)
T_ss = (T_0 + sqrt((T_0^2)-(4*B)))/2 #Compute the steady state temperature of the sphere
print(T_ss)

#B) Plot the water surface temperature as a function of time
viscosity = CP.HAPropsSI("Visc","T",T_inf,"P",p,"R",rh_inf) #Compute the dynamic viscosity of air (Pa-s)
Re = v*l/viscosity #Calculate the reynolds number of air
Pr = viscosity/alpha #Calculate the prandtl number
Nu = 0.43*(Re^0.58)*(Pr^0.4) #Calculate the nusselt number, this is geometry dependent
Sh = Nu #The sherwood number is the nussselt number, heat and mass transfer analogy

h_m = Sh*(D_ab/l) #Compute the mass transfer coefficient, m/s
Sc = viscosity/D_ab #Compute the smith number



 function temp_vs_time(temp, time) #Function to tell us temperature at each time step
     dTdt = [1.0] #initialize with some value
     rho_air_s = CP.PropsSI("D","T",temp[1],"Q",1,"Water") #Compute the air density at the surface (kg/m^3)
     dTdt[1] = -6.0*(h_m*(temp[1]-T_inf))/(rho_air_s*cp_air*l) #change in temperature at each time step, taking the derivative of final_time function
     return dTdt
 end

 times = 0:2:100 #go from 0 by 2's to 100 seconds
 sol_t = int.odeint(temp_vs_time, T_0, times)

 plot(sol_t,
    ylabel = "Temperature (K)",
    xlabel = "Time (s)")
print(sol_t)

function mass_vs_time(mass, time) #Function to tell us water mass at each time step
    dMdt = [1.0] #initialize with some value
    rho_air_s = CP.PropsSI("D","T",sol_t[2],"Q",1,"Water") #Compute the air density at the surface (kg/m^3)
    dMdt[1] = -h_m*a_s*(rho_air_s-rho_air)  #change in mass at each time step
    return dMdt
end
times = 0:2:100 #go from 0 by 2's to 100 seconds
sol_m = int.odeint(mass_vs_time, m_0, times)

plot(sol_m,
    ylabel = "Mass (kg)",
    xlabel = "Time (s)")

mean(sol_t[2])

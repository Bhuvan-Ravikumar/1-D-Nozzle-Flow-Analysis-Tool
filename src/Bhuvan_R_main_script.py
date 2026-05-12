# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 17:09:32 2024

@author: bhuva
"""

#
#
#

import numpy as np
from math import *
import csv
import matplotlib.pyplot as plt
import func
import standard_atmosphere

 #Our bisection function
def my_bisection(f, a, b, tol, nmax): 
   
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("The scalars a and b do not bound a root")
    for n in range(1,nmax):
      m = (a + b)/2
      #print(n,a,b,m,f(m))
      if np.abs(f(m)) < tol:
        return m
      elif np.sign(f(a)) == np.sign(f(m)):
        a=m
      elif np.sign(f(b)) == np.sign(f(m)):
        b=m
  # to find supersonic mach no. from area ratio
def find_my_mach_sup(Ae_A_star, gamma, tol, nmax):
    f = lambda x: x - (1/Ae_A_star)*((2/(gamma+1))*(1+((gamma-1)/2)*x**2))**((gamma+1)/(2*(gamma-1)))
    a=1
    b=10
    Me = my_bisection(f,a,b,tol,nmax)
    return Me
 
    # to find subsonic mach no. from area ratio
def find_my_mach_sub(Ae_A_star, gamma, tol, nmax):
    f = lambda x: x - (1/Ae_A_star)*((2/(gamma+1))*(1+((gamma-1)/2)*x**2))**((gamma+1)/(2*(gamma-1)))
    a=0
    b=1
    Me = my_bisection(f,a,b,tol,nmax)
    return Me

   # to find area ratios from know mach no.
def find_my_area_ratio(Me,gamma):
  area_ratio= (1/Me)*((2/(gamma+1))*(1+((gamma-1)/2)*Me**2))**((gamma+1)/(2*(gamma-1)))
  return area_ratio

   # to find Po/P for known exit mach
def pressure_ratio(Me,gamma):
  pres_ratio= (1 + ((gamma-1)/2)*Me**2)**(gamma/(gamma-1))
  return pres_ratio

   # to find To/T for know mach no.
def temperature_ratio(Me,gamma):
    temp_ratio = (1 + ((gamma-1)/2)*Me**2)
    return temp_ratio

   # to find density ratio
def density_ratio(Me,gamma):
    rho_ratio = (1 + ((gamma-1)/2)*Me**2)**(1/(gamma-1))
    return rho_ratio

   # to find momentum thrust
def momentum_thrust(Pe_by_Po,gamma):
  return sqrt((2*gamma**2/(gamma-1))*((2/(gamma+1))**((gamma+1)/(gamma-1)))*(1-Pe_by_Po**((gamma-1)/gamma)))

   # to find pressure thrust
def pressure_thrust(Ae_A_star,Pe_by_Po,Pa_by_Po):
  pres = (Pe_by_Po - Pa_by_Po) * Ae_A_star
  return pres

   # to find the thrust co-efficient
def thrust_cof(Ae_A_star,Pe_by_Po,Pa_by_Po,gamma):
  momen_T = momentum_thrust(Pe_by_Po,gamma)
  pres_T = pressure_thrust(Ae_A_star,Pe_by_Po,Pa_by_Po)
  C_f= momen_T + pres_T
  return C_f



# Opening the given CSV file
with open('Ravi Kumar_Bhuvan_nozzle_geometry.csv', 'r') as file_one:
    csv_reader = list(csv.reader(file_one,delimiter = ","))
    
# Converting list to arrays
csv_data = np.array(csv_reader[1:],dtype=float)
axial_pos_x = csv_data[:,0]
norm_nozzle_radius = csv_data[:,1]


# Load data from a text file
with open('Ravi Kumar_Bhuvan_input.txt', 'r') as file_two:
    lines = file_two.readlines()

gamma = None
thrust_n = None
ambi_pressure = None
chamber_pressure = None
chamber_temperature = None
Mol_mass = None

# Processing and storing each data

for line in lines:
    div = line.split('!')
    value = div[0].strip()
    data = div[1].strip() if len(div) > 1 else None
    if data == 'Nominal Thrust (N)':
        thrust_n = float(value)
    elif data == 'Ambient Pressure (Pa)':
        ambi_pressure = float(value)
    elif data == 'Chamber Pressure (Pa)':
        chamber_pressure = float(value)
    elif data == 'Chamber Temperature (K)':
        chamber_temperature = float(value)
    elif data == 'Specific Heat Ratio gamma':
        gamma = float(value)
    elif data == 'Molar Mass (g/mol)':
        Mol_mass = float(value)

# Converting mol mass into kg/mol
Mol_mass = Mol_mass/1000

# Physical Constants
R_uni = 8.314
g_o = 9.80665

R = R_uni/Mol_mass

chamber_density = chamber_pressure/(chamber_temperature*R)
cp = gamma*R/(gamma-1)

# Task 1 computation

datapt = len(norm_nozzle_radius)

norm_area = np.zeros(norm_nozzle_radius.shape)
for i in range(len(norm_nozzle_radius)):
    norm_area[i]= pi*(norm_nozzle_radius[i])**2
A_star = min(norm_area)
A_by_A_star = np.zeros(norm_area.shape)
for i in range(datapt):
    A_by_A_star[i] = norm_area[i]/A_star

area_mach= np.zeros(A_by_A_star.shape)
for i in range(datapt):
    if axial_pos_x[i] <= 0:
        area_mach[i] = find_my_mach_sub(A_by_A_star[i], gamma,0.00001,1000)
    else:
        area_mach[i] = find_my_mach_sup(A_by_A_star[i], gamma,0.00001,1000)

Mach = area_mach

P_o = chamber_pressure*pressure_ratio(Mach[0], gamma)
T_o = chamber_temperature*temperature_ratio(Mach[0], gamma)
rho_o = chamber_density*density_ratio(Mach[0], gamma)
h_o = cp*T_o

# Task 2 computation

static_pres = np.zeros(Mach.shape)

iso_pres = np.zeros(Mach.shape)
for i in range(datapt):
    iso_pres[i] = P_o/pressure_ratio(Mach[i], gamma)

static_pres = iso_pres
P_2 = 0
P_2 = iso_pres[-1]*(1+(2*gamma/(gamma+1))*((Mach[-1]**2)-1))

if iso_pres[-1] < ambi_pressure:
    
    if P_2 > ambi_pressure:
        print('oblique shock outside nozzle')
    elif abs(P_2 - ambi_pressure) < 100:
        static_pres[-1] = P_2
        Mach[-1] = sqrt((1/(gamma-1)) + sqrt((1/(gamma-1)**2)+(2/(gamma-1))*((2/(gamma+1))**((gamma+1)/(gamma-1)))*((P_o/P_2)*(1/A_by_A_star[-1]))**2))
        print('normal shock at nozzle exit')
    else:
        Me_2 = sqrt((1/(gamma-1)) + sqrt((1/(gamma-1)**2)+(2/(gamma-1))*((2/(gamma+1))**((gamma+1)/(gamma-1)))*((P_o/P_2)*(1/A_by_A_star[-1]))**2))
        Po_2 = P_2*(1+((gamma-1)/1)*Me_2**2)**(gamma/(gamma-1))
        A_shock = func.normal_shock_Area(Po_2, P_o, gamma)
        A_star_2 = (P_o/Po_2)*A_star
        index = 0
        tol= 0.01
        for i in range(datapt):
            if abs(A_by_A_star[i]-A_shock) < tol:
                index = i+1
                break
        for index in range(datapt):
            A_by_A_star_2[index] = norm_area[index]/A_star_2
            Mach[index]= find_my_mach_sub(A_by_A_star_2[index], gamma,0.00001,1000)
            static_pres[index] = Po_2/pressure_ratio(Mach[index], gamma)
        print('normal shock inside nozzle')




static_temp = np.zeros(Mach.shape)
for i in range(datapt):
    static_temp[i] = T_o/temperature_ratio(Mach[i], gamma)

a = np.zeros(static_temp.shape)
for i in range(datapt):
    a[i] = sqrt(gamma*R*static_temp[i])

u = np.zeros(Mach.shape)
for i in range(datapt):
    u[i] = Mach[i]*a[i]

static_density = np.zeros(Mach.shape)
for i in range(datapt):
    static_density[i] = rho_o/density_ratio(Mach[i], gamma)

static_h = np.zeros(Mach.shape)
for i in range(datapt):
    static_h[i] = cp*static_temp[i]


# Sizing nozzle dimension and computing task 3 values


Pe_Po_ratio = 1/(pressure_ratio(area_mach[-1],gamma))
Pa_Po_ratio = ambi_pressure/P_o
Cf= thrust_cof(A_by_A_star[-1],Pe_Po_ratio,Pa_Po_ratio,gamma)
cee_star = sqrt((1/gamma)*(((gamma+1)/2)**((gamma+1)/(gamma-1)))*R*T_o)
choked_m_dot = thrust_n/(Cf*cee_star)
A_t_star = sqrt(T_o)*choked_m_dot/(P_o*sqrt((gamma/R)*(2/(gamma+1))**((gamma+1)/(gamma-1))))
throat_radius = sqrt(A_t_star/pi)


noz_radius = np.zeros(norm_nozzle_radius.shape)
for i in range(datapt):
    noz_radius[i] = norm_nozzle_radius[i]*throat_radius

Area = np.zeros(norm_nozzle_radius.shape)
for i in range(datapt):
    Area[i] = pi*(norm_nozzle_radius[i]*throat_radius)**2



# task 2 plotting graphs

plt.figure()
plt.plot(axial_pos_x,Area)
plt.title("Area along the nozzle")
plt.xlabel('normalized nozzle length X')
plt.ylabel('Area (m$^2$) ')

plt.figure()
plt.plot(axial_pos_x,Mach)
plt.title("Mach along the nozzle")
plt.xlabel('normalized nozzle length X')
plt.ylabel('Mach no.')

plt.figure()
plt.plot(axial_pos_x,u)
plt.title("flow velocity along the nozzle")
plt.xlabel('normalized nozzle length X')
plt.ylabel('flow velocity U (ms$^-1$)')

plt.figure()
plt.plot(axial_pos_x,static_pres)
plt.title("Static pressure along the nozzle")
plt.xlabel('normalized nozzle length X')
plt.ylabel('Static Pressure P (Pa)')
plt.axhline(P_o,ls='--',color="red", label='stagnation pressure')
plt.legend()

plt.figure()
plt.plot(axial_pos_x,static_temp)
plt.title("static temperature along the nozzle")
plt.xlabel('normalized nozzle length X')
plt.ylabel('Static Temperature T (K)')
plt.axhline(T_o,ls='--',color="yellow", label='stagnation temperature')
plt.legend()

plt.figure()
plt.plot(axial_pos_x,static_density)
plt.title("Static density along the nozzle")
plt.xlabel('normalized nozzle length X')
plt.ylabel('Static Density (kg/m$^3$)')
plt.axhline(rho_o,ls='--',color="black", label='stagnation density')
plt.legend()

plt.figure()
plt.plot(axial_pos_x,static_h)
plt.title("Static enthalpy along the nozzle")
plt.xlabel('normalized nozzle length X')
plt.ylabel('Static enthalpy h (J/kg)')
plt.axhline(h_o,ls='--',color="green", label='stagnation enthalpy')
plt.legend()

plt.show()

# task 3 tabulation
Pe = static_pres[-1]
Te = static_temp[-1]
h_e = static_h[-1]
Ae_by_A_star = A_by_A_star[-1]
L_div = 10*throat_radius
Me = Mach[-1]
I_sp = (Cf*cee_star)/g_o

print('Stagnation pressure P0 (Pa) = ' , '%.4f' % P_o)
print('Stagnation temperature T0 (K) = ' , '%.4f' % T_o)
print('Stagnation enthalpy h0 (J/kg) = ' , '%.4f' % h_o)
print('Nozzle exit pressure Pe (Pa) = ' , '%.4f' % Pe)
print('Nozzle exit temperature Te (K) = ' , '%.4f' % Te)
print('Nozzle exit enthalpy he (J/kg) = ' , '%.4f' % h_e)
print('Nozzle area ratio Ae/At = ', '%.4f' % Ae_by_A_star)
print('Nozzle length (diverging portion) from the throat to the exit Ldiv (m) = ' , '%.4f' % L_div)
print('Exit Mach number Me = ' , '%.4f' % Me)
print('Mass flow rate ˙m (kg/s)= ', '%.4f' % choked_m_dot)
print('Characteristic velocity c∗ (m/s) = ' , '%.4f' % cee_star)
print('Thrust Coefficient Cf = ' , '%.4f' % Cf)
print('Specific Impulse Isp (s) = ' , '%.4f' % I_sp)

# task 4

mach_sub = np.zeros(datapt)
for i in range(datapt):
    mach_sub[i] = find_my_mach_sub(A_by_A_star[i], gamma, 0.00001, 1000)

Pe_isen = iso_pres[-1]


under_exp_pressure = np.empty(len(area_mach))
isen_exp_pressure_sub = np.empty(len(mach_sub))
for i in range(datapt):
    isen_exp_pressure_sub[i] = P_o/pressure_ratio(mach_sub[i], gamma)
Pe_isen_sub = isen_exp_pressure_sub[-1]
isen_exp_pressure_sup = np.empty(len(area_mach))
over_exp_pressure_obq = np.empty(len(area_mach))
over_exp_pressure_normal_t = np.empty(len(area_mach))
over_exp_pressure_normal_1 = np.empty(len(area_mach))

A_by_A_star_2 = np.empty(len(norm_area))

Pb_range = np.array([0.0, Pe_isen, ambi_pressure, P_2, 311950.0132, Pe_isen_sub], dtype= np.float64)


for i in range(len(Pb_range)):
    if Pb_range[i] < Pe_isen:
        under_exp_pressure = iso_pres
    elif abs(Pb_range[i] - Pe_isen) < 10:
        isen_exp_pressure_sup = iso_pres
    elif Pb_range[i] > Pe_isen and Pb_range[i] < Pe_isen_sub:
        if Pb_range[i] < P_2:
            over_exp_pressure_obq = iso_pres
        elif abs(Pb_range[i] - P_2) < 10:
            over_exp_pressure_normal_t = iso_pres
            over_exp_pressure_normal_t[-1] = P_2
        elif Pb_range[i] > P_2:
            Me_2 = sqrt(-(1/(gamma-1)) + sqrt((1/(gamma-1)**2)+(2/(gamma-1))*((2/(gamma+1))**((gamma+1)/(gamma-1)))*((P_o/P_2)*(1/A_by_A_star[-1]))**2))
            Po_2 = Pb_range[i]*(1+((gamma-1)/2)*Me_2**2)**(gamma/(gamma-1))
            A_shock = func.normal_shock_Area(Po_2, P_o, gamma)
            A_star_2 = (P_o/Po_2)*A_star
            index = 0
            tol= 0.01
            for j in range(datapt):
                if abs(A_by_A_star[j]-A_shock) < tol:
                    index = j+1
                    break
            for k in range(datapt):
                if k == index:
                    A_by_A_star_2[index] = norm_area[index]/A_star_2
                    Mach[index]= find_my_mach_sub(A_by_A_star_2[index], gamma,0.00001,1000)
                    over_exp_pressure_normal_1[index] = Po_2/pressure_ratio(Mach[index], gamma)
    elif abs(Pb_range[i] - Pe_isen_sub) < 10:
        isen_exp_pressure_sub = isen_exp_pressure_sub



plt.figure(num=0,dpi=75)
plt.plot(axial_pos_x,under_exp_pressure, label='under expanded')
plt.plot(axial_pos_x,isen_exp_pressure_sub, label='isentropic subsonic')
plt.plot(axial_pos_x,isen_exp_pressure_sup, label='isntropic supersonic')
plt.plot(axial_pos_x,over_exp_pressure_normal_1, label='over exp w/ normal shock')
plt.plot(axial_pos_x,over_exp_pressure_normal_t, label='over exp w/ shock at exit')
plt.plot(axial_pos_x,over_exp_pressure_obq, label='over expanded oblique')
plt.title("Static nozzle pressure profile for different Pb")
plt.xlabel('normalized nozzle length X')
plt.ylabel('Pressure (Pa)')
plt.legend()
plt.show()

# task 5

alt_list = np.linspace(0, 86, 200)

Pa_vaccum = 0.0
Pa_list = np.zeros(len(alt_list))
for i in range(len(alt_list)):
    alt_data = standard_atmosphere.atmosphere(alt_list[i])
    Pa_list[i]= alt_data[0]

Pa_vaccum_Po = Pa_vaccum/P_o
Pa_by_Po_list = np.zeros(len(alt_list))
for i in range(len(alt_list)):
    Pa_by_Po_list[i]= Pa_list[i]/P_o

Cf_alt = np.zeros(len(alt_list))
for i in range(len(alt_list)):
    Cf_alt = thrust_cof(A_by_A_star[-1],Pe_Po_ratio,Pa_by_Po_list,gamma)

Cf_vaccum = thrust_cof(A_by_A_star[-1], Pe_Po_ratio, Pa_vaccum_Po, gamma)

Isp_alt = np.zeros(len(alt_list))
for i in range(len(alt_list)):
    Isp_alt[i]= Cf_alt[i]*cee_star/g_o

Isp_vaccum = Cf_vaccum*cee_star/g_o

plt.figure()
plt.plot(Cf_alt, alt_list)
plt.axvline(Cf_vaccum,ls='--',color="green", label='vaccum Cf')
plt.title("Cf vs Altitude")
plt.xlabel('Cf')
plt.ylabel('Altitude (km)')
plt.legend()
plt.show

plt.figure()
plt.plot(Isp_alt, alt_list)
plt.axvline(Isp_vaccum,ls='--',color="red", label='vaccum Isp')
plt.title("Isp vs Altitude")
plt.xlabel('Isp (s)')
plt.ylabel('Altitude (km)')
plt.legend()
plt.show














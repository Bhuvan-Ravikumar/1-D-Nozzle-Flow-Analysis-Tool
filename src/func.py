# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 22:40:20 2024

@author: bhuva
"""

import numpy as np
from math import *

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


def normal_shock_Area(Po_2,Po_1,gamma):
    f = lambda x: Po_2/Po_1 - ((((gamma+1)*x**2)/(((gamma-1)*x**2)+2))**(gamma/(gamma-1))) * (((gamma+1)/((2*gamma*x**2)-(gamma-1)))**(1/(gamma-1)))
    a = 1
    b = 10
    M1 = my_bisection(f, a, b, 0.00001, 1000)
    A_shock_A_star = (1/M1)*(((2/(gamma+1))*(1+((gamma-1)/2)*M1**2))**((gamma+1)/(2*(gamma-1))))
    return A_shock_A_star
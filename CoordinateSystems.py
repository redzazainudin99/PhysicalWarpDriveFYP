import numpy as np
import sympy as sp
from sympy import Matrix,shape,symbols, diff, sin, cos, tan, sinh, tanh, acos, atan, sqrt, limit, oo

#this file specifies the coordinate bases used in the programs.

#symbols for coordinate bases

t, x, y, z = sp.symbols('t x y z')
t, r, theta, phi = sp.symbols('t r theta phi')


#Minkowski coordinates

xminkow = [t, x, y, z]
XMinkowFuncSpher = [t,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta)]

#Spherical coordinates

XSpherFuncMinkow = [t, sqrt((x**2)+(y**2)+(z**2)), acos(z/(sqrt((x**2)+(y**2)+(z**2)))),atan(y/x)]
xspher = [t, r, theta, phi]

#creating pairs for utility (0 is a placeholder term) (the first term is an index denoting the basis)
Minkowski = [1,xminkow, XMinkowFuncSpher]
Spherical = [2,XSpherFuncMinkow, xspher]

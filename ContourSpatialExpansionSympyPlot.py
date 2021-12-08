import sympy as sp
import numpy as np
from sympy import Matrix, Array, shape,symbols, diff, sin, cos, tan, sinh, tanh, acos, atan, sqrt, limit, oo, ln, evalf

#importing the previous module
import CalculationsForEinsteinTensor as ein
from TensorModule import *
from CoordinateSystems import *

import matplotlib.pyplot as plt
from sympy.plotting import plot3d


#warp drive parameters
from WarpDriveParameters import *

print('Overall tensor calculations start here')

#substituting values of coordinates
ScalarCWarp =  ein.ScalarCOverall.subs(r,sqrt((x**2)+(y**2)+(z**2))).subs(theta, acos(z/(sqrt((x**2)+(y**2)+(z**2))))).subs(phi, atan(y/x))
#print(ScalarCWarp)

def Hotelling(Ricci,n,E,Z,o_o,R_o,m_o,G_o,c_o,pi_o):
    ans = (1 + Ricci*(E**2)/(6*n))
    answer = ans.subs(ein.z,Z).subs(ein.o,o_o).subs(ein.R,R_o).subs(ein.m,m_o).subs(ein.G,G_o).subs(ein.c,c_o).subs(ein.pi,pi_o)
    return answer
    print('Success')
    pass

print('Expansion calculations start here')

output = Hotelling(ScalarCWarp,n,E,Z,o_o,R_o,m_o,G_o,c_o,pi_o).evalf()

print('Starting to plot')

plott = plot3d(output, (x, -outputrange, outputrange), (y, -outputrange, outputrange),  title = 'Plot pengembangan ruang pemacu', xlabel = 'x',ylabel = 'y') #(example 300 in all directions)
plott[0].surface_color = 'r'
plott.show()

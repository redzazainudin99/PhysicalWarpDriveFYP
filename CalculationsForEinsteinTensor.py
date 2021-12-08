import numpy as np
import sympy as sp

from sympy.printing.mathml import print_mathml, mathml
from sympy import Matrix, Array, shape,symbols, diff, sin, cos, tan, sinh, cosh, tanh, acos, atan, sqrt, limit, oo, latex
sp.init_printing(use_unicode = True)

#establishing special symbols.
v = sp.symbols('v')                                     #this is a symbol denoting the velocity of the warp drive
R, m =sp.symbols('R m')                                 #these symbols represent the radius and mass of the shell
c, G = sp.symbols('c G')                                #these symbols represent the speed of light and the gravitational constant
pi = sp.symbols('pi')

#custom modules
from CoordinateSystems import *
from TensorModule import *

#defining the inner and outer tensors of the metric
#establishing inner metric tensor
G_in = Matrix([[-(c**2)*(1 - (2*G*m/((c**2)*R))), 0, 0, 0],[0, 1, 0, 0],[0, 0, ( r**2),0],[0,0,0,((r**2)*(sin(theta))**2)]])
InnerMetric = TwoRankTensor(G_in, Spherical)

#establishing the outer metric tensor
G_out = Matrix([[-(c**2)*(1 - (2*G*m/((c**2)*R))), 0, 0, 0],[0, (1 - (2*m/r))**-1, 0, 0],[0, 0, (r**2),0],[0,0,0,((r**2)*(sin(theta))**2)]])
OuterMetric = TwoRankTensor(G_out, Spherical)

ChrisSymIn = ChristoffelSym(InnerMetric)
RCTIn = RiemannCurvature(InnerMetric, ChrisSymIn)   #for Schwarzschild coordinates, there is a single unsimplified term that, when simplified, is equal to zero
RicciTIn = RicciTensor(InnerMetric,RCTIn)
ScalarCIn = RicciScalar(InnerMetric, RicciTIn)

ChrisSymOut = ChristoffelSym(OuterMetric)
RCTOut = RiemannCurvature(OuterMetric, ChrisSymOut)
RicciTOut = RicciTensor(OuterMetric,RCTOut)
ScalarCOut = RicciScalar(OuterMetric, RicciTOut)

#with this, proof that the inner region is flat
RCTIn

#the outer region is flat at the limit of r = oo
RCTOutLimit = ZeroArrayCreate (4,4,4,4)
for a in range(0,4):
    for b in range(0,4):
        for c in range(0,4):
            for d in range(0,4):
                RCTOutLimit[a,b,c,d] = limit(RCTOut[a,b,c,d], r, oo)



#creating the approximation matrix
o, R = sp.symbols('o R')                                                    #defining the terms that are parameters of the shape function
ShapeFunc = (tanh(o*(r+R)) - tanh(o*(r-R))) / 2*tanh(o*R)                 #this function denotes the internal, basic shape function

C = sp.exp(o*R)
A = (C**2) - (1/C)**2
B = (C**2) + (1/C)**2
ShapeFunc_0 = ((2/((C+(1/C))**2))*((B/2) + cosh(2*o*r)))**-1                      #this is a simplification of ShapeFunc
ShapeFunc_1 = 2*A / (B + sp.exp(2*o*r) + sp.exp(-2*o*r))                            #another simplification
SSF = OuterMetric.arr*(1 - ShapeFunc_1) + InnerMetric.arr * ShapeFunc_1     #defining the selective shape function

#this describes the overall metric tensor of the entire warp drive
OverallMetric = TwoRankTensor(SSF, Spherical)

#determining characteristic tensors
ChrisSymOverall = ChristoffelSym(OverallMetric)
RCTOverall = RiemannCurvature(OverallMetric, ChrisSymOverall, simplifyState = 'off')
RicciTOverall = RicciTensor(OverallMetric,RCTOverall, simplifyState = 'off')
ScalarCOverall = RicciScalar(OverallMetric, RicciTOverall, simplifyState = 'off')

#calculating the Einstein tensor
EinsteinT = ZeroMatrixCreate(4,4)
for u in range(0,4):
    for v in range(0,4):
        EinsteinT[u,v] = (RicciTOverall[u,v] - 0.5*ScalarCOverall*OverallMetric.arr[u,v])

MassEnergyTensor = EinsteinT/((8*pi*G)/(c**4))

MassDistribution = MassEnergyTensor[0,0]
MassDistribution

out_file = open("EinsteinTensor00.txt","w")
out_file.write(latex(sp.factor(MassDistribution)))
out_file.close()

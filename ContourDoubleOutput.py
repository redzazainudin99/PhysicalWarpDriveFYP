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

#importing outputs of spatial distortion and mass density
import ContourSpatialExpansionSympyPlot as space
import ContourDensitySympyPlot as mass


plott = plot3d(space.output, mass.output, (x, -outputrange, outputrange), (y, -outputrange, outputrange),title = 'Plot gabungan pengembangan ruang dan ketumpatan jisim', xlabel = 'x',ylabel = 'y') #(example 300 in all directions)
plott[0].surface_color = 'g'
plott[1].surface_color = 'r'
plott.show()

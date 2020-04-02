from casadi import *
from numpy import *
import time

from fourierControl import *
from systemModel import *
from lagrangianjacobian import *

#Pseudoinverse version

model = systemModel()
inverse = LagrangianJacobian(model)

time = 20
yd = np.array([0,0.8,0.5])
e = model.getOutput() - yd

i = 0

while numpy.linalg.norm(e) > 0.001:
    print(numpy.linalg.norm(e))
    inverseJ = inverse.inverseLagrangianJacobian(time)
    print("--------------------------------------------------------------")
    print("Step: ",i+1)
    print("Control: ",model.getControl())
    print("Configuration:", model.getConfiguration())
    e = model.getOutput() - yd
    print("Error: ",e)
    inverseJe = mtimes(inverseJ,e)
    gamma = inverse.getGamma(inverseJe)
    change = gamma*inverseJe
    print(gamma)
    model.fourierControl.actualizeControlVector(model.getControl()-change)
    olde = e
    i = i+1

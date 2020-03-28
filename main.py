from casadi import *
from numpy import *
import time

from fourierControl import *
from systemModel import *
from lagrangianjacobian import *

model = systemModel()
inverse= LagrangianJacobian(model)

time = 20
yd = np.array([0.5,-0.8,2])
e = model.getOutput() - yd

i = 0
while abs(e[0]) > 0.001 or abs(e[1]) > 0.001 or abs(e[2]) > 0.001:
    print("--------------------------------------------------------------")
    print("Step: ",i+1)
    print("Control: ",model.getControl())
    print("Configuration:", model.getConfiguration())
    e = model.getOutput() - yd
    print("Error: ",e)
    inverseJe = mtimes(inverse.inverseLagrangianJacobian(time),e)
    gamma = inverse.getGamma(inverseJe)
    change = gamma[0]*inverseJe
    change = np.array(change).astype(np.float64).flatten()
    model.fourierControl.actualizeControlVector(model.getControl()-change)
    i = i+1

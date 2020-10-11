from lagrangianjacobian import *

# Lagrangian inverse version
from lagrangianjacobian import LagrangianJacobian

model = systemModel()
inverse: LagrangianJacobian = LagrangianJacobian(model)

yd: ndarray = desPos
e = model.getOutput - yd
i: int = 0

while numpy.linalg.norm(e) > maxError:
    print(numpy.linalg.norm(e))
    print('--------------------------------------------------------------')
    print('Step: ', i + 1)
    print('Control: ', model.getControl)
    print('Control value at the beginning: ', model.getControlValue(0))
    print('Control value at the end: ', model.getControlValue(Th))
    print("Configuration:", model.getConfiguration)
    e = model.getOutput - yd
    print('Error: ', e)
    inverseJ = inverse.inverseLagrangianJacobian(Th)
    inverseJe = mtimes(inverseJ, vertcat(e, DM([[0], [0]])))
    gamma = inverse.getGamma(inverseJe)
    print('Gamma: ', gamma)
    change = gamma * inverseJe
    model.fourierControl.actualizeControlVector(model.getControl - change)
    i = i + 1

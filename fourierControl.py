import sys
import math
from casadi import *

t = SX.sym("t")
omega = 2*math.pi/20

class fourierControl:
    def __init__(self,init_coefficients,controlVectorSize):
        self.size = len(init_coefficients)
        self.controlVectorSize = controlVectorSize
        self.coefficients = SX(init_coefficients)
        self.P = SX.zeros(self.controlVectorSize,self.size)
        self.FourierFunctions = list([1,sin(omega*t),cos(omega*t)])
        for i in range(controlVectorSize):
            for j in range(self.size):
                if (int)(j/3) == i:
                    self.P[i,j] = self.FourierFunctions[j%3]
        self.Us = mtimes(self.P,self.coefficients)
        self.u = Function('u',[t],[self.Us])

    def Pmatrix(self):
        return self.P

    def actualizeControlVector(self,changed_coefficients):
        self.coefficients = SX(changed_coefficients)
        self.Us = mtimes(self.P,self.coefficients)
        self.u = Function('u',[t],[self.Us])

    def getControlCoefficients(self):
        return self.coefficients

    def evaluate(self,time):
        return self.u(time)

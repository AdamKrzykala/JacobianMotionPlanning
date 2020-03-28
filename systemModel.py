import sys
import math
from sympy import *
from fourierControl import *
from casadi import *
from numpy import *

ConfDim = 3
ControlDim = 2
OutputDim = 3
controlVectorSize = 6

q = SX.sym("q",ConfDim)

class systemModel:
    def __init__(self):
        global ConfDim, ConfDim, OutputDim
        #-----------------------INPUT DATA---------------------------
        self.initialConfiguration = np.array([0,0,0])
        self.qTab = self.initialConfiguration
        self.initialControl = np.array([2,2,2,2,2,2])
        self.fourierControl = fourierControl(self.initialControl,ControlDim)
        self.Gs = vertcat(
            horzcat(cos(q[2]),    0),
            horzcat(sin(q[2]),    0),
            horzcat(        0,    1))
        self.Ys = vertcat(
            q[0],
            q[1],
            q[2])
        self.p = Function('p',[t],[self.fourierControl.Pmatrix()])
        self.g = Function('g',[q],[self.Gs])
        self.y = Function('y',[q],[self.Ys])

    def getConfiguration(self):
        return self.qTab

    def getOutput(self):
        return self.y(self.qTab)

    def getControl(self):
        return self.fourierControl.getControlCoefficients()

    def conf(self):
        return self.qTab

    def actualizeConfiguration(self,new_configuration):
        self.qTab = new_configuration

    def backToInitialConfiguration(self):
        self.qTab = self.initialConfiguration

    def Gsolver(self,tempq):
        return self.g(tempq)

    def Asolver(self,tempq,time):
        temp = mtimes(self.Gs,self.fourierControl.evaluate(time))
        func = jacobian(temp,q)
        f = Function('f',[q],[func])
        return f(tempq)

    def Bsolver(self,tempq):
        return self.g(tempq)

    def Csolver(self,tempq):
        func = jacobian(self.Ys,q)
        f = Function('f',[q],[func])
        return f(tempq)

    def Ysolver(self,tempq):
        return self.y(tempq)

    def Psolver(self,time):
        return self.p(time)

    def EndogenousConfigurationODE(self,q,t):
        return mtimes(self.Gsolver(q),self.fourierControl.evaluate(t))

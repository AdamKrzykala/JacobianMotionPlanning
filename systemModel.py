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
        M = 60.017
        m1 = 4.5
        m2 = 1.5
        l1 = 0.619
        l2 = 0.6
        d1 = 0.313
        d2 = 0.287
        I = 2.384
        B = (m1*m2*(l1-d1)**2+M*(m1*d1**2+m2*l1**2))/(M+m1+m2)
        C = ((M+m1)*m2*d2**2)/(M+m1+m2)
        D = (m1*m2*(l1-d1)*d2+M*m2*l1*d2)/(M+m1+m2)

        self.initialConfiguration = np.array([0,0,pi/8])
        self.qTab = self.initialConfiguration
        self.initialControl = np.array([0.05,0,0,0,0,0])
        self.fourierControl = fourierControl(self.initialControl,ControlDim)
        self.Gs = vertcat(
            horzcat(-(B+C+2*D*cos(q[2]))/(I+B+C+2*D*cos(q[2])),-(C+D*cos(q[2]))/(I+B+C+2*D*cos(q[2]))),
            horzcat(        1,    0),
            horzcat(        0,    1))
        self.Ys = vertcat(
            q[0],
            l1*cos(q[1])+l2*cos(q[1]+q[2]),
            l1*sin(q[1])+l2*sin(q[1]+q[2]))
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

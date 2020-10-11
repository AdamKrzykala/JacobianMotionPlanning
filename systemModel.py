from sympy import *

from fourierControl import *
from taskDef import *


class systemModel:
    qTab: ndarray

    def __init__(self):
        """
        @rtype: systemModel
        """
        global ConfDim, ConfDim, OutputDim
        # -----------------------INPUT DATA---------------------------

        self.initialConfiguration = InitConfiguration
        self.qTab = self.initialConfiguration
        self.initialControl = InitControl
        self.fourierControl = fourierControl(self.initialControl, ControlDim)
        self.Gs = generatorsMatrix
        self.Fs = initialDrift
        self.Ys = outputSys
        self.p = Function('p', [t], [self.fourierControl.Pmatrix()])
        self.g = Function('g', [q], [self.Gs])
        self.f = Function('f', [q], [self.Fs])
        self.y = Function('y', [q], [self.Ys])

    @property
    def getConfiguration(self):
        """
        @rtype: object
        """
        return self.qTab

    def getControlValue(self, time: float):
        """
        @type time: float
        @rtype: object
        """
        return self.fourierControl.evaluate(time)

    @property
    def getOutput(self):
        return self.y(self.qTab)

    @property
    def getControl(self):
        return self.fourierControl.getControlCoefficients()

    def conf(self):
        return self.qTab

    def actualizeConfiguration(self, new_configuration):
        self.qTab = new_configuration

    def backToInitialConfiguration(self):
        self.qTab = self.initialConfiguration

    def Gsolver(self, tempq):
        return self.g(tempq)

    def Fsolver(self, tempq):
        return self.f(tempq)

    def Asolver(self, tempq, time):
        temp = mtimes(self.Gs, self.fourierControl.evaluate(time))
        func = jacobian(temp, q)
        f = Function('f', [q], [func])
        return f(tempq)

    def Bsolver(self, tempq):
        return self.g(tempq)

    def Csolver(self, tempq):
        func = jacobian(self.Ys, q)
        f = Function('f', [q], [func])
        return f(tempq)

    def Ysolver(self, tempq):
        return self.y(tempq)

    def Psolver(self, time):
        return self.p(time)

    def EndogenousConfigurationODE(self, q, t):
        return mtimes(self.Gsolver(q), self.fourierControl.evaluate(t)) + self.Fsolver(q)

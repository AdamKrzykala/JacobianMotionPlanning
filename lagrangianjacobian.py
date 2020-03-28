import numpy as np
import sys
import math
from casadi import *
from systemModel import *

U = SX.sym("U",ConfDim,controlVectorSize)
V = SX.sym("V",controlVectorSize,controlVectorSize)

class LagrangianJacobian:
    def __init__(self,model):
        self.baseModel = model
        self.baseModel.backToInitialConfiguration()
        self.Gamma0 = 0.1

    def UlambdaODE(self,q,t,U):
        return mtimes(self.baseModel.Bsolver(q),self.baseModel.Psolver(t))+mtimes(self.baseModel.Asolver(q,t),U)

    def VlamdaODE(self,Ulamda,q,t):
        Vt =  mtimes(mtimes(Ulamda.T,self.Q(q,t)),Ulamda)
        Vt += mtimes(mtimes(Ulamda.T,self.S(q,t)),self.baseModel.Psolver(t))
        Vt += mtimes(mtimes(self.baseModel.Psolver(t).T,self.S(q,t).T),Ulamda)
        Vt += mtimes(mtimes(self.baseModel.Psolver(t).T,self.R(q)),self.baseModel.Psolver(t))
        return Vt

    def Vlambda(self,time):
        self.baseModel.backToInitialConfiguration()
        Vtinit = mtimes(mtimes(self.baseModel.Psolver(0).T,self.R(self.baseModel.getConfiguration())),self.baseModel.Psolver(0))
        a = controlVectorSize*controlVectorSize
        X = vertcat(casadi.reshape(V,a,1),casadi.reshape(U,ConfDim*controlVectorSize,1),q)
        ode = vertcat(casadi.reshape(self.VlamdaODE(U,q,t),a,1),
                    casadi.reshape(self.UlambdaODE(q,t,U),ConfDim*controlVectorSize,1),
                    self.baseModel.EndogenousConfigurationODE(q,t))
        dae = {'t':t,'x':X,'ode':ode}
        opts = {}
        opts["t0"] = 0
        opts["tf"] = 1
        F = integrator("F", "cvodes", dae, opts)
        sol = F(x0=vertcat(casadi.reshape(Vtinit,a,1),
                            DM.zeros(ConfDim*controlVectorSize,1),
                            [0,0,0]))
        VT = casadi.reshape(sol['xf'][0:a],controlVectorSize,controlVectorSize)
        UT = casadi.reshape(sol['xf'][a:a+ConfDim*controlVectorSize],ConfDim,controlVectorSize)
        qT = casadi.reshape(sol['xf'][a+ConfDim*controlVectorSize:a+ConfDim*controlVectorSize+3],3,1)
        return VT,UT,qT

    def getGamma(self,dLam):
        gmm = 2*self.Gamma0
        nrm = np.linalg.norm(np.array(dLam).astype(np.float64), ord=2)
        return [min(sqrt(gmm/nrm),1)]

    def Q(self,q,t):
        a = self.baseModel.Asolver(q,t)
        return mtimes(a.T,a)

    def S(self,q,t):
        return mtimes(self.baseModel.Asolver(q,t).T,self.baseModel.Bsolver(q))

    def R(self,q):
        return mtimes(self.baseModel.Bsolver(q).T,self.baseModel.Bsolver(q))

    def inverseLagrangianJacobian(self,timeHorizon):
        Vt,Ut,q = self.Vlambda(timeHorizon)
        Ct = self.baseModel.Csolver(q)
        self.baseModel.actualizeConfiguration(q)
        V = inv(Vt)
        sol = mtimes(V,Ut.T)
        sol = mtimes(sol,Ct.T)
        sol = mtimes(sol,Ct)
        sol = mtimes(sol,Ut)
        sol = mtimes(sol,V)
        sol = mtimes(sol,Ut.T)
        sol = mtimes(sol,Ct.T)
        return sol
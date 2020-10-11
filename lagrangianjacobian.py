from systemModel import *

U: casadi.SX = SX.sym('U', ConfDim, controlVectorSize)
V: casadi.SX = SX.sym('V', controlVectorSize, controlVectorSize)


class LagrangianJacobian:
    baseModel: systemModel
    Gamma0: float

    def __init__(self, model: systemModel):
        """
        @type model: systemModel
        """
        assert isinstance(model, systemModel)
        self.baseModel = model
        self.baseModel.backToInitialConfiguration()
        assert isinstance(InitGamma, float)
        self.Gamma0 = InitGamma

    def UlambdaODE(self, q, t, U):
        return mtimes(self.baseModel.Bsolver(q), self.baseModel.Psolver(t)) + mtimes(self.baseModel.Asolver(q, t), U)

    def VlamdaODE(self, Ulamda, q, t):
        Vt = mtimes(mtimes(Ulamda.T, self.Q(q, t)), Ulamda)
        Vt += mtimes(mtimes(Ulamda.T, self.S(q, t)), self.baseModel.Psolver(t))
        Vt += mtimes(mtimes(self.baseModel.Psolver(t).T, self.S(q, t).T), Ulamda)
        Vt += mtimes(mtimes(self.baseModel.Psolver(t).T, self.R(q)), self.baseModel.Psolver(t))
        return Vt

    def Vlambda(self, time):
        self.baseModel.backToInitialConfiguration()
        Vtinit = mtimes(mtimes(self.baseModel.Psolver(0).T, self.R(self.baseModel.getConfiguration)),
                        self.baseModel.Psolver(0))
        a = controlVectorSize * controlVectorSize
        X = vertcat(casadi.reshape(V, a, 1), casadi.reshape(U, ConfDim * controlVectorSize, 1), q)
        ode = vertcat(casadi.reshape(self.VlamdaODE(U, q, t), a, 1),
                      casadi.reshape(self.UlambdaODE(q, t, U), ConfDim * controlVectorSize, 1),
                      self.baseModel.EndogenousConfigurationODE(q, t))
        dae = dict(t=t, x=X, ode=ode)
        F = integrator("F", "cvodes", dae, dict(t0=0, tf=time))
        sol = F(x0=vertcat(casadi.reshape(Vtinit, a, 1),
                           DM.zeros(ConfDim * controlVectorSize, 1),
                           InitConfiguration))
        VT = casadi.reshape(sol['xf'][0:a], controlVectorSize, controlVectorSize)
        UT = casadi.reshape(sol['xf'][a:a + ConfDim * controlVectorSize], ConfDim, controlVectorSize)
        qT = casadi.reshape(sol['xf'][a + ConfDim * controlVectorSize:a + ConfDim * controlVectorSize + ConfDim],
                            ConfDim, 1)
        return VT, UT, qT

    def getGamma(self, dLam):
        gmm = 2 * self.Gamma0
        nrm = np.linalg.norm(np.array(dLam), ord=2)
        return [min(sqrt(gmm / nrm) / 100, 1)]

    def Q(self, q, t):
        a = self.baseModel.Asolver(q, t)
        return mtimes(a.T, a)

    def S(self, q, t):
        return mtimes(self.baseModel.Asolver(q, t).T, self.baseModel.Bsolver(q))

    def R(self, q):
        return mtimes(self.baseModel.Bsolver(q).T, self.baseModel.Bsolver(q))

    def inverseLagrangianJacobian(self, timeHorizon: float) -> 'inversedJacobian':
        """
        @rtype: 'inversedJacobian'
        @type timeHorizon: float
        """
        Vt, Ut, q = self.Vlambda(timeHorizon)
        # Utt = mtimes(Ut.T, inv(mtimes(Ut, Ut.T)))
        # print(Utt.size())
        # print(Vt.size())
        K = DM([[1,0,1,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,1,0,1]])
        Ct = self.baseModel.Csolver(q)
        self.baseModel.actualizeConfiguration(q)
        Vinv = inv(Vt)
        J = mtimes(Ct,Ut)
        JK = vertcat(J, K)
        sol = mtimes(Vinv, JK.T)
        w = mtimes(JK, Vinv)
        w = mtimes(w, JK.T)
        sol = mtimes(sol, inv(w))
        return sol

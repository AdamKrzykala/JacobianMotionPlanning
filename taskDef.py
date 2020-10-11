from casadi import *
from numpy import *

from numpy.core._multiarray_umath import ndarray

# ------------------------------------- USER TASK AND MODEL CONFIGURATION ----------------------------------------------

# Time horizon
Th: float = 10.0

# Dimension of configuration space
ConfDim: int = 3

q: casadi.SX = SX.sym("q", ConfDim)

# Amount of controls
ControlDim: int = 2

# Amount of outputs
OutputDim: int = 3

# Dimension of control vector
controlVectorSize: int = 10

# Initial configuration

InitConfiguration: ndarray = np.array([0, 0, pi / 8])

# Initial control vector
InitControl: ndarray = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

# USER VARIABLES
# -----------------------------------------------------------------------------------------------------------------------

M = 60.017
m1 = 4.5
m2 = 1.5
l1 = 0.619
l2 = 0.6
d1 = 0.313
d2 = 0.287
Ib = 2.384
A = m1 + m2 + M
a = 0.3
b = 0.2
Ip = M*(m1+m2)*(a**2 + b**2)/A
B = (m1*m2*(l1-d1)**2+M*(m1*d1**2+m2*l1**2))/A
C = (m1+M)*m2*d2**2/A
D = (m1*m2*(l1-d1)*d2+M*m2*l1*d2)/A
E = M*(m1*d1+m2*l1)/A
F = M*m2*d2/A
p = 0.1

# -----------------------------------------------------------------------------------------------------------------------

#   Generators - G Matrix

generatorsMatrix = \
    vertcat(
        horzcat(-(B+C+2*D*cos(q[2])+E*(a*cos(q[1])+b*sin(q[1]))+F*(a*cos(q[1]+q[2])+b*sin(q[1]+q[2])))/(Ib+Ip+B+C+2*D*cos(q[2])+2*E*(a*cos(q[1])+b*sin(q[1]))+2*F*(a*cos(q[1]+q[2])+b*sin(q[1]+q[2]))),
                -(C+D*cos(q[2])+F*(a*cos(q[1]+q[2])+b*sin(q[1]+q[2]))) / (Ib+Ip+B+C+2*D*cos(q[2])+2*E*(a*cos(q[1])+b*sin(q[1]))+2*F*(a*cos(q[1]+q[2])+b*sin(q[1]+q[2])))),
        horzcat(1, 0),
        horzcat(0, 1)
    )

initialDrift = \
    vertcat(
        p/(Ib+Ip+B+C+2*D*cos(q[2])+2*E*(a*cos(q[1])+b*sin(q[1]))+2*F*(a*cos(q[1]+q[2])+b*sin(q[1]+q[2]))),
        0,
        0
    )

#   Output
outputSys = \
    vertcat(
        q[0],
        a + l1 * cos(q[1]) + l2 * cos(q[1] + q[2]),
        b + l1 * sin(q[1]) + l2 * sin(q[1] + q[2])
    )

# Initial gamma coefficient
InitGamma: float = 0.1

# Desired output when time horizon elapsed
desPos: ndarray = np.array([0, 0.9, -0.5])

# Maximal allowed error
maxError: float = 0.0001

# ------------------------------------------- END OF USER CONFIGURATION ------------------------------------------------

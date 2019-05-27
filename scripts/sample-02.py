import numpy as np
from genosolver import minimize

"""
Sample for constrained problem
"""


def fg(x):
    f = (x[0] - 1)**2 + 2 * (x[1] - 2)**2
    g = np.array([2 * (x[0] - 1), 4 * (x[1] - 2)])
    return f, g


def c(x):
    return np.array([x[0] + 4 * x[1]])


def cjprod(x, v):
    return np.array([1, 4]) * v

x0 = np.array([0, 0])
# y0 is optional
y0 = 0
bounds = None

cl = 3
cu = 3
constraints = {'type' : 'bnds',
               'fun' : c,
               'jacprod' : cjprod,
               'cl' : cl, 
               'cu' : cu}

tol = 1E-8
# any option is optionable
options = {'maxiter' : 50,
           'constraintsTol' : 1E-5,
           'disp' : False,
           'debug_fg' : False,  
           'debug_c' : False, 
           'debug_cjprod' : False}

#constraints = []
result = minimize(fg, x0, y0, bounds = bounds, tol = tol, 
                  constraints = constraints, options = options)
print(result)

"""
runfile('C:/Users/Soeren/Desktop/soeren/c++/genosolver/src/pygeno/scripts/sample-02.py', wdir='C:/Users/Soeren/Desktop/soeren/c++/genosolver/src/pygeno/scripts')
Reloaded modules: genosolver, pygeno
 elapsed: 0.00839401303525733
     fun: 3.9999904459346016
     jac: array([-1.33333174, -5.33332696])
  nInner: 30
    nfev: 45
     nit: 7
 success: True
       x: array([ 0.33333413,  0.66666826])
       y: array([ 1.33332458])
"""
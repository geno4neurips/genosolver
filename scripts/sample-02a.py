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
    M = np.array([[1, 4], 
                  [2, 0],
                  [0, 3]])
    return M.dot(x)


def cjprod(x, v):
    M = np.array([[1, 4], 
                  [2, 0],
                  [0, 3]])
    return v.dot(M)

x0 = np.array([0, 0])
y0 = np.array([0, 0, 0])
bounds = None

cl = np.array([3, 2, -10])
cu = np.array([3, 10, 10])
constraints = {'type' : 'bnds',
               'fun' : c,
               'jacprod' : cjprod,
               'cl' : cl, 
               'cu' : cu}

tol = 1E-8
options = {'maxiter' : 20,
           'constraintsTol' : 1E-5,
           'disp' : False,
           'debug_fg' : False,  
           'debug_c' : False, 
           'debug_cjprod' : False}

#constraints = []
result = minimize(fg, x0, y0, bounds = bounds, tol = tol, 
                  constraints = constraints, options = options)
print(result)

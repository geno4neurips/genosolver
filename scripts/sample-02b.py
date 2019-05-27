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
    return np.array([x[0] + 4 * x[1]]) - 3


def cjprod(x, v):
    return np.array([1, 4]) * v

x0 = np.array([0, 0])
# y0 is optional
y0 = 0
bounds = None

# constraint == 0
constraints = {'type' : 'eq',
               'fun' : c,
               'jacprod' : cjprod}
# constraint <= 0
constraints = {'type' : 'ineq',
               'fun' : c,
               'jacprod' : cjprod}

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

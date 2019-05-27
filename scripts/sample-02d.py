import numpy as np
from genosolver import minimize

"""
Sample for constrained problem

min c'*x
st. A*x <= b
    x >= 0

"""


def fg(x):
    c = np.array([-1., -2.])
    f = c.dot(x)
    g = c
    return f, g

A = np.array([[1., 3.],
              [2., 1.],
              [4., 3.]])
b = np.array([3., 3., 6.])

def c(x):
    return A.dot(x)


def cjprod(x, v):
    return v.dot(A)

x0 = np.array([0, 0])
bounds = [(0, np.inf)] * 2
y0 = 0

constraints = {'type' : 'bnds',
               'fun' : c,
               'jacprod' : cjprod,
               'cl' : -np.inf,
               'cu' : b}

tol = 1E-8
# any option is optionable
options = {'maxiter' : 500,
           'constraintsTol' : 1E-5,
           'disp' : False,
           'debug_fg' : False,
           'debug_c' : False,
           'debug_cjprod' : False}

#constraints = []
result = minimize(fg, x0, y0, bounds = bounds, tol = tol,
                  constraints = constraints, options = options)
print(result)
x = result.x
print(np.dot(A, x) - b)
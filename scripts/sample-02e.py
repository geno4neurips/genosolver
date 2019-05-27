import numpy as np
from genosolver import minimize

"""
Sample for constrained problem

min c'*x
st. A*x <= b
    sum(x) == 1
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

def c1(x):
    return A.dot(x)


def c1jprod(x, v):
    return v.dot(A)

def c2(x):
    return np.sum(x)

def c2jprod(x, v):
    u = np.array([1., 1.])
    return v * u

x0 = np.array([0, 0])
bounds = [(0, np.inf)] * 2

constraints = ({'type' : 'bnds',
                'fun' : c1,
                'jacprod' : c1jprod,
                'cl' : -np.inf,
                'cu' : b},
               {'type' : 'bnds',
                'fun' : c2,
                'jacprod' : c2jprod,
                'cl' : 1.4,
                'cu' : 1.4})


tol = 1E-12
# any option is optionable
options = {'maxiter' : 100,
           'constraintsTol' : 1E-8,
           'disp' : False,
           'debug_fg' : False,
           'debug_c' : False,
           'debug_cjprod' : False}

#constraints = []
result = minimize(fg, x0, bounds = bounds, tol = tol,
                  constraints = constraints, options = options)
print(result)
x = result.x
print(np.dot(A, x) - b)
print(np.sum(x))

import numpy as np
from genosolver import minimize

"""
Sample for unconstraint problem and verbose options
"""


def fg(x):
    return sum(x**2), 2 * x


x0 = np.array([100, 100], dtype=np.double)

bounds = [(1, np.inf)] * 2
tol = 1E-6
options = {'disp' : True,
           'debug_fg' : False}

result = minimize(fg, x0, bounds = bounds, tol = tol, options = options)
print(result)

# -*- coding: utf-8 -*-
import warnings
import genosolver.genointerface as genointerface
from scipy.optimize import OptimizeResult
import numpy as np
from timeit import default_timer

def minimize(fg, x0, y0 = None, bounds = None, tol = 1E-8, constraints = (),
             options = {}):
    start = default_timer()
    options = dict(options)
    options['tol'] = tol
    options.setdefault('num_cor', 10)
    options.setdefault('maxiter', 10000)
    options.setdefault('constraintsTol', 1E-5)
    options.setdefault('disp', False)
    options.setdefault('debug_fg', False)
    options.setdefault('debug_c', False)
    options.setdefault('debug_cjprod', False)

    message = ['Optimal solution found.',
               'Suboptimal solution found.',
               'Problem is unbounded.',
               'Problem is infeasible.',
               'Solver encountered internal numerical error.',
               'Starting point causes numerical error.',
               'Stopped by user.']

    x0 = np.ascontiguousarray(x0, dtype = np.float64)
    n = x0.size

    if bounds is not None:
        bnds = np.ascontiguousarray(bounds, dtype = np.float64)
        lb, ub = np.ascontiguousarray(bnds[:, 0]), np.array(bnds[:, 1])
    else:
        lb = np.full(n, -np.inf)
        ub = np.full(n, np.inf)

    if constraints is None:
        constraints = ()

    if isinstance(constraints, dict):
        constraints = (constraints, )

    solver = genointerface.Geno(options['debug_fg'],
                         options['debug_c'],
                         options['debug_cjprod'])
    numSymbolicConstraints = len(constraints)
    def secureFG(x):
        f, g = fg(x)
        f = np.float64(f)
        g = np.ascontiguousarray(g, dtype = np.float64)
        return f, g

    if numSymbolicConstraints == 0:
        res = solver.solve_unconstrained(secureFG, x0, lb, ub,
                                               options['tol'],
                                               options['num_cor'],
                                               options['maxiter'],
                                               options['disp'])
        (x, f, g, status, nIter, funEval) = res
        elapsed = default_timer() - start
        return OptimizeResult(x = x, fun = f, jac = g, success = (status <= 2),
                              nit = nIter, nfev = funEval, elapsed = elapsed,
                              status = status, message=message[status])
    else:
        numConstraints = []
        shapeConstraints = []
        offset = [0]
        clAll = []
        cuAll = []
        for c in constraints:
            dummyFConstraints = c['fun'](x0)
            dummyFConstraints = np.ascontiguousarray(dummyFConstraints,
                                                     dtype = np.float64)
            shapeConstraints.append(dummyFConstraints.shape)
            m = len(dummyFConstraints.reshape(-1))
            numConstraints.append(m)
            offset.append(offset[-1] + m)

            # check the type of constraints
            if c['type'] == 'eq':
                cl = 0
                cu = 0
            elif c['type'] == 'ineq':
                cl = -np.inf
                cu = 0
            elif c['type'] == 'bnds':
                cl = c['cl']
                cu = c['cu']
            else:
                assert(False)

            cl = np.ascontiguousarray(cl, dtype = np.float64)
            cu = np.ascontiguousarray(cu, dtype = np.float64)
            if len(cl) == 1 and not m == 1:
                cl = np.full(m, cl[0], dtype = np.float64)
            if len(cu) == 1 and not m == 1:
                cu = np.full(m, cu[0], dtype = np.float64)

            clAll.append(cl)
            cuAll.append(cu)

        mTotal = offset[-1]
        clAll = np.concatenate(clAll)
        cuAll = np.concatenate(cuAll)

        def allC(x):
            l = [np.ascontiguousarray(c['fun'](x).reshape(-1), dtype = np.float64) for c in constraints]
            f = np.concatenate(l)
            assert(f.flags['C_CONTIGUOUS'])
            return f
        def allCjprod(x, v):
            g = np.zeros_like(x)
            for i, c in enumerate(constraints):
                g = g + np.ascontiguousarray(c['jacprod'](x, v[offset[i]:offset[i+1]].reshape(shapeConstraints[i])).reshape(-1), dtype = np.float64)
            assert(g.flags['C_CONTIGUOUS'])
            return g

        if y0 is None:
            y0 = np.full(mTotal, 0, dtype = np.float64)
        y0 = np.ascontiguousarray(y0, dtype = np.float64)
        if len(y0) == 1 and not mTotal == 1:
            y0 = np.full(mTotal, y0[0], dtype = np.float64)

        # add jacprod field if it is missing but jac is present
        # TODO: maybe issue a warning
        # TODO: works only for scalars and vectors, not matrices
        jacprodMissing = False
        for i, c in enumerate(constraints):
            if not 'jacprod' in c and 'jac' in c:
                jacprodMissing = True
                if numConstraints[i] == 1:
                    c['jacprod'] = lambda x, v, jac=c['jac']: v * np.ascontiguousarray(jac(x))
                else:
                    c['jacprod'] = lambda x, v, jac=c['jac']: np.dot(v, jac(x))

        if jacprodMissing:
            warnings.warn('jacprod missing')

        assert(x0.flags['C_CONTIGUOUS'])
        assert(lb.flags['C_CONTIGUOUS'])
        assert(ub.flags['C_CONTIGUOUS'])
        assert(y0.flags['C_CONTIGUOUS'])
        assert(cl.flags['C_CONTIGUOUS'])
        assert(cu.flags['C_CONTIGUOUS'])
        res = solver.solve_constrained(secureFG,
                                       x0,
                                       lb, ub,
                                       y0,
                                       allC,
                                       allCjprod,
                                       clAll,
                                       cuAll,
                                       options['tol'],
                                       options['constraintsTol'],
                                       options['num_cor'],
                                       options['maxiter'],
                                       options['disp'])
        (x, y, f, g, slack, status, nIter, funEval, nInner) = res

        elapsed = default_timer() - start
        return OptimizeResult(x = x, y = y, fun = f,
                              jac = g, success = (status <= 2),
                              nit = nIter, nfev = funEval, nInner = nInner,
                              elapsed = elapsed, status=status,
                              message=message[status],
                              slack=slack)

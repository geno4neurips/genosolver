#defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
"""
 Note: Not thread safe ;)
"""
cimport cython
import numpy as np
cimport numpy as np
from libcpp.limits cimport numeric_limits
from libcpp cimport bool
from libc.string cimport memcpy
np.import_array()

#from timeit import default_timer as timer

cdef object global_fg
cdef object global_c
cdef object global_cjprod
cdef int global_num_var
cdef int global_num_const
ctypedef void (*callback_fg_type) (const double *x, double &f, double *g)
ctypedef void (*callback_c_type) (const double *x, double *c)
ctypedef void (*callback_cjprod_type) (const double *x, const double *duals,
                                       double *g)
ctypedef bool (*callback_iter_type) (const double *x)

cdef bool verbose_fg = False
cdef bool verbose_c = False
cdef bool verbose_cjprod = False

cdef bool stopped = False

cdef double timeFAndG
cdef double timeFAndG2

cdef bool callback_iter(const double *x):
    return stopped

cdef void callback_fg(const double *x, double &f, double *g):
#    global timeFAndG
#    global timeFAndG2
#    start = timer()
#    start2 = timer()
    cdef double[:] x_mv
    cdef np.ndarray pyx
    cdef np.ndarray[double, ndim=1, mode='c'] pyg2
    try:
        x_mv = <double[:global_num_var]> x
        pyx = np.asarray(x_mv)
        if verbose_fg:
            print("callback_fg(*x,&f,*g)")
            print("x=",pyx)
#        start2 = timer()
        pyf, pyg = (<object>global_fg)(pyx)
#        elapsed2 = timer() - start2
#        timeFAndG2 += elapsed2
        if verbose_fg:
            print("f=", pyf)
            print("g=", pyg)
        (&f)[0] = <double> pyf
        pyg2 = pyg
        memcpy(g, &pyg2[0], global_num_var * sizeof(double))
    except KeyboardInterrupt:
        print('User pressed CTRL-C.')
        global stopped
        stopped = True
#    elapsed = timer() - start
#    timeFAndG += elapsed

cdef void callback_c(const double *x, double *c):
    cdef double[:] x_mv = <double[:global_num_var]> x
    cdef np.ndarray pyx = np.asarray(x_mv)
    cdef np.ndarray[double, ndim=1, mode='c'] pyc2
    try:
        x_mv = <double[:global_num_var]> x
        pyx = np.asarray(x_mv)
        if verbose_c:
            print("callback_c(*x, *c)")
            print("x=",pyx)
        pyc = (<object>global_c)(pyx)
        if verbose_c:
            print("c=",pyc)
        pyc2 = pyc
        memcpy(c, &pyc2[0], global_num_const * sizeof(double))
    except KeyboardInterrupt:
        print('User pressed CTRL-C.')
        global stopped
        stopped = True

cdef void callback_cjprod(const double *x, const double *v, double *g):
    cdef double[:] x_mv = <double[:global_num_var]> x
    cdef np.ndarray pyx = np.asarray(x_mv)

    cdef double[:] v_mv = <double[:global_num_const]> v
    cdef np.ndarray pyv = np.asarray(v_mv)
    cdef np.ndarray[double, ndim=1, mode='c'] pyg2
    try:
        x_mv = <double[:global_num_var]> x
        pyx = np.asarray(x_mv)

        v_mv = <double[:global_num_const]> v
        pyv = np.asarray(v_mv)

        if verbose_cjprod:
            print("callback_cjprod(*x, *v, *g)")
            print("x=", pyx)
            print("v=", pyv)
        pyg = (<object>global_cjprod)(pyx, pyv)
        if verbose_cjprod:
            print("g=", len(pyg), global_num_var, pyg)
        pyg2 = pyg
        memcpy(g, &pyg2[0], global_num_var * sizeof(double))
    except KeyboardInterrupt:
        print('User pressed CTRL-C.')
        global stopped
        stopped = True



cdef class Geno:
    def __init__(self, debug_fg=False, debug_c=False, debug_cjprod=False):
        global verbose_fg
        verbose_fg = debug_fg
        global verbose_c
        verbose_c = debug_c
        global verbose_cjprod
        verbose_cjprod = debug_cjprod
    def solve_constrained(self,
                         fg,
                         double[:] x0,
                         double[:] lb,
                         double[:] ub,
                         double[:] y0,
                         c,
                         cjprod,
                         double[:] cl,
                         double[:] cu,
                         double tol,
                         double constraintsTol,
                         int num_corr,
                         int maxiter,
                         bool verbose):
        global global_num_var
        global global_num_const
        global global_fg
        global global_cjprod
        global global_c
        global_num_var = x0.size
        global_num_const = y0.size
        global_fg = fg
        global_c = c
        global_cjprod = cjprod
        global stopped
        stopped = False
        cdef np.ndarray x = np.ones(global_num_var, dtype=np.double)
        cdef np.ndarray y = np.ones(global_num_const, dtype=np.double)
        cdef double f = 0
        cdef np.ndarray g = np.ones(global_num_var, dtype=np.double)
        cdef double slack = 0
        cdef int status = 1
        cdef int nIter = 0
        cdef int funEval = 0
        cdef int nInner = 0
        solveConstrained(<callback_fg_type> callback_fg,
                         <callback_iter_type> callback_iter,
                         global_num_var,
                         &x0[0],
                         &lb[0],
                         &ub[0],
                         global_num_const,
                         &y0[0],
                         &cl[0],
                         &cu[0],
                         <callback_c_type> callback_c,
                         <callback_cjprod_type> callback_cjprod,
                         num_corr,
                         <double *> x.data,
                         <double *> y.data,
                         f,
                         <double *> g.data,
                         slack,
                         status,
                         nIter,
                         funEval,
                         nInner,
                         tol,
                         constraintsTol,
                         maxiter,
                         verbose)
        return x, y, f, g, slack, status, nIter, funEval, nInner


    def solve_unconstrained(self, fg, double[:] x0,
                            double[:] lb, double[:] ub,
                            double tol,
                            int num_corr, int maxiter,
                            bool verbose):
        global global_num_var
        global_num_var = x0.size
        cdef np.ndarray x = np.ones(global_num_var, dtype=np.double)
        cdef double f = 0
        cdef np.ndarray g = np.ones(global_num_var, dtype=np.double)
        cdef int status = 1
        cdef int nIter = 0
        cdef int funEval = 0
        global global_fg
        global_fg = fg
        global stopped
        stopped = False
#        global timeFAndG
#        timeFAndG = 0.
#        global timeFAndG2
#        timeFAndG2 = 0.
        solve(<callback_fg_type> callback_fg,
              <callback_iter_type> callback_iter,
              global_num_var,
              &x0[0],
              &lb[0],
              &ub[0],
              num_corr,
              <double *> x.data,
              f,
              <double *> g.data,
              status,
              nIter,
              funEval,
              tol,
              maxiter,
              verbose)
#        print('cython fAndG  %.3f' % timeFAndG)
#        print('cython fAndG2 %.3f' % timeFAndG2)
        return x, f, g, status, nIter, funEval
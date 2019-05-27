from libcpp cimport bool

cdef extern from "pygenointerface.hpp":
    void solve(void (*fg)(const double *x, double &f, double *g),
               bool (*callback)(const double *x),
		       int num_var,
		       double *x0,
		       double *lb,
		       double *ub,
		       int num_corr, # num correction pairs
		       double *x,	 # solution
            double &f,
            double *g,
		       int &status,
		       int &iter,
		       int &funEval,
		       double tol,
		       int maxiter,
		       bool display);

    void solveConstrained(void (*fg)(const double* x, double &f, double *g),
        bool (*callback)(const double *x),
		   int num_var,
		   double *x0,
		   double *lb,
		   double *ub,
		   int num_const,
		   double *y0,
		   double *lbC,
		   double *ubC,
		   void (*cf) (const double*, double*),
		   void (*cJprod)(const double*, const double*, double*),
		   int num_corr,
        double *x,	 # solution
        double *y,	 # dual solution
        double &f,
        double *g,
        double &slack,
		   int &status,
		   int &iter,
		   int &funEval,
        int &inner,
		   double tol,
        double constraintsTol,
		   int maxiter,
		   bool display);
#include <iostream>
#include <functional>
#include "pygenointerface.hpp"
#include "pygenonlp.hpp"
#include "lbfgsb.hpp"
#include "augmentedLagrangian.hpp"


void solve(std::function<void(const Scalar*, Scalar&, Scalar*)> fg,
           std::function<bool(const Scalar*)> callback,
		   int num_var,
		   double *x0,
		   double *lb,
		   double *ub,
		   int num_corr,
		   double *x,
		   double &f,
		   double *g,
		   int &status,
		   int &iter,
		   int &funEval,
		   double tol,
		   int maxiter,
		   bool verbose)
{
	PyGenoNLP gpnlp(num_var, 0, x0, lb, ub, fg, callback);
	LBFGSB solver(gpnlp, num_corr, verbose);
	solver.setParameter("pgtol", tol);
	solver.setParameter("maxiter", maxiter);
	status = solver.solve();
	f = solver.f();
	const double *xConst = solver.x();
	const double *gConst = solver.g();
	for (int i = 0; i < num_var; ++i) {
		x[i] = xConst[i];
		g[i] = gConst[i];
	}
	iter = solver.iter();
	funEval = solver.funEval();
}


void solveConstrained(std::function<void(const Scalar*, Scalar&, Scalar*)> fg,
        std::function<bool(const Scalar*)> callback,
		   int num_var,
		   double *x0,
		   double *lb,
		   double *ub,
		   int num_const,
		   double *y0,
		   double *lbC,
		   double *ubC,
		   std::function<void(const Scalar*, Scalar*)> cf,
		   std::function<void(const Scalar*, const Scalar*, Scalar*)> cJprod,
		   int num_corr, // num correction pairs
		   double *x,
		   double *y,
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
		   bool verbose) {

	PyGenoNLP gpnlp(num_var, num_const, x0, lb, ub, fg, callback,
		            y0, lbC, ubC, cf, cJprod);
	AugmentedLagrangian solver(gpnlp, num_corr, verbose);
	solver.setParameter("pgtol", tol);
	solver.setParameter("maxiter", maxiter);
	solver.setParameter("constraintsTol", constraintsTol);
	status = solver.solve();
	f = solver.f();
	const double *xConst = solver.x();
	const double *gConst = solver.g();
	for (int i = 0; i < num_var; ++i) {
		x[i] = xConst[i];
		g[i] = gConst[i];
	}
	const double *yConst = solver.y();
	for (int i = 0; i < num_const; ++i) {
		y[i] = yConst[i];
	}
	slack = solver.slack();
	iter = solver.iter();
	funEval = solver.funEval();
	inner = solver.inner();
}

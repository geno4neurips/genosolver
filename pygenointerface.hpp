#ifndef __PYGENOINTERFACE_HPP__
#define __PYGENOINTERFACE_HPP__

#include <iostream>
#include <functional>
#include "pygenonlp.hpp"
#include "lbfgsb.hpp"
#include "augmentedLagrangian.hpp"

void solve(std::function<void(const Scalar*, Scalar&, Scalar*)> fg,
           std::function<bool(const Scalar*)> callback,
		   int num_var,
		   double *x0,
		   double *lb,
		   double *ub,
		   int num_corr, // num correction pairs
		   double *x,
		   double &f,
		   double *g,
		   int &status,
		   int &iter,
		   int &funEval,
		   double tol,
		   int maxiter,
		   bool verbose);

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
		   bool verbose);
#endif
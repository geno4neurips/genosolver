#pragma once
#ifndef AUGMENTED_LAGRANGIAN_H
#define AUGMENTED_LAGRANGIAN_H

#include "lbfgsb.hpp"

class AugmentedLagrangian
{
public:
  AugmentedLagrangian(GenoNLP& genoNLP, size_t correctionPairs,
                      bool verbose = false);

  SolverStatus solve();

  // Set x for warm start.
  //  const Vector& x(const Vector& x);
  const Scalar* x() const;
  const Scalar* y() const;
  const Scalar* g() const;
  Scalar f() const;
  Scalar slack() const;
  int iter() const;
  int funEval() const;
  int inner() const;

  // Set parameters
  bool setParameter(std::string parameter, Scalar value);

private:
  // disable default constructor
  AugmentedLagrangian();
  // disable copy constructor
  AugmentedLagrangian(const AugmentedLagrangian& other);

//  double computeConstraintErrorNorm(const AugmentedNLP& augmentedNLP);
  GenoNLP& _genoNLP;
  double _constraintsTol;
  double _c;
  size_t _correctionPairs;
  Scalar _f;
  Vector _x;
  Vector _augLagG;
  Vector _y;
  Scalar _slack;
  size_t _maxIter;
  size_t _iter;
  size_t _funEval;
  size_t _inner;

  Scalar _tol;
  Scalar _tolFun;
  bool _verbose;
};

#endif // AUGMENTED_LAGRANGIAN_H

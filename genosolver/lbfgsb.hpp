#pragma once
#ifndef LBFGSB_HPP_
#define LBFGSB_HPP_

#include <utility>

#include <Eigen/Dense>

#include "genoNLP.hpp"

typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

enum SolverStatus {SOLVED = 0, SUBOPTIMAL, UNBOUNDED, INFEASIBLE, NUM_ERROR,
                   INIT_POINT_ERROR, CALLBACK_STOP};

class LBFGSB
{
public:
  LBFGSB(GenoNLP& genoNLP, Index m, bool verbose = true);
//  ~LBFGSB();

  void restart();
  // TODO
  // write function that does one iteration
  // and call it within solve()
  // iterate()

  // TODO
  // pass continueSolve call back
  SolverStatus solve();

  // Set x for warm start.
  const Scalar* x(const Scalar* x);
  const Scalar* x() const;
  const Scalar* g() const;
  Scalar f() const;
  int iter() const;
  int funEval() const;

  // Set solver parameters.
  bool setParameter(std::string parameter, Scalar value);

private:
  void boundInit();
  void updateB();
  void cauchyPoint();
  void choleskyFactorK();
  void formK();
  void subsm();
  int lineSearchMT(const Vector& deltaX, size_t iter, Scalar& stp);
  void projectedGradient();
  Scalar helper();
  void refresh();
  bool initFeasible();
  void projectFeasible(Vector& v);
  bool hasNanInf(Scalar f);
  bool hasNanInf(const Vector& g);

  GenoNLP& _genoNLP;
  Scalar _f;
  Vector _x;
  Vector _lb;
  Vector _ub;
  Vector _g;

  size_t _iter;
  size_t _funEval;
  bool _verbose;
  size_t _maxIter;
  Scalar _tol;
  Scalar _tolFun;
  Index _m;
  Index _index;
  Vector _s;
  Vector _y;
  Matrix _S;
  Matrix _Y;
  Scalar _theta;
  Matrix _SY;
  Matrix _SS;
  Matrix _M;
  Matrix _Q;
  Matrix _K;
  Matrix _L;
  Matrix _E;

  Matrix _SYws;
  Matrix _SSws;
  Matrix _YYws;

  Vector _xCP;
  Vector _workingSet;
  Vector _workingSetOld;
  Vector _c;
  Vector _xSubspace;
  Vector _projectedGradient;
  Scalar _fNew;
  Vector _xNew;
  Vector _gNew;
  bool _constrained;
  bool _boxed;
};

#endif // LBFGSB_HPP_

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include "lbfgsb.hpp"
#include "augmentedLagrangian.hpp"

static const Scalar INF = std::numeric_limits<Scalar>::infinity();

class AugmentedNLP : public GenoNLP
{
public:
  AugmentedNLP(GenoNLP& genoNLP, Scalar rho, const Vector& y)
    :
    _genoNLP(genoNLP),
    _n(genoNLP.getN()),
    _m(genoNLP.getM()),
    _rho(rho),
    _y(y),
    _cl(Vector(genoNLP.getM())),
    _cu(Vector(genoNLP.getM())),
    _constraintError(Vector(genoNLP.getM())),
    _augmentedConstraintError(Vector(genoNLP.getM()))
  {
    _genoNLP.getBoundsConstraints(_cl.data(), _cu.data());
  }

  virtual ~AugmentedNLP()
  {
  }

  virtual bool getInfo(Index& n, Index& m)
  {
    n = _n;
    m = 0;
    return true;
  }

  virtual bool getBounds(Scalar* lb, Scalar* ub)
  {
    _genoNLP.getBounds(lb, ub);
    return true;
  }

  virtual bool getStartingPoint(Scalar* x)
  {
    _genoNLP.getStartingPoint(x);
    return true;
  }

  void computeConstraintError(const Scalar* x)
  {
    Vector constraintValues(_m);
    _genoNLP.functionValueConstraints(x, constraintValues.data());

    _augmentedConstraintError =
      (constraintValues - _cl + _y / _rho).cwiseMin(0.0) +
      (constraintValues - _cu + _y / _rho).cwiseMax(0.0);
    _constraintError =
      (constraintValues - _cl).cwiseMin(0.0) +
      (constraintValues - _cu).cwiseMax(0.0);
  }

  virtual bool functionValueAndGradient(const Scalar* x,
                                        Scalar& functionValue,
                                        Scalar* gradient)
  {
    _genoNLP.functionValueAndGradient(x, functionValue, gradient);
    computeConstraintError(x);

    functionValue += _rho/2 * _augmentedConstraintError.squaredNorm();

    Vector::MapType gradientMapX(gradient, _n);
    Vector gradientX(_n);

    Vector v = _rho * _augmentedConstraintError;
    _genoNLP.gradientConstraintsTimesVector(x, v.data(), gradientX.data());

    gradientMapX += gradientX;
    return true;
  }

  virtual bool callback(const Scalar *x)
  {
    return _genoNLP.callback(x);
  }

  const Vector& constraintError()
  {
    return _constraintError;
  }

  const Vector& augmentedConstraintError()
  {
    return _augmentedConstraintError;
  }

  void setParams(Scalar rho, const Vector& y)
  {
    _rho = rho;
    _y = y;
  }

  size_t m()
  {
    return _m;
  }

private:
  AugmentedNLP(const AugmentedNLP&);
  void operator=(const AugmentedNLP&);

  GenoNLP& _genoNLP;
  size_t _n;
  size_t _m;
  Scalar _rho;
  Vector _y;
  Vector _cl;
  Vector _cu;
  Vector _constraintError;
  Vector _augmentedConstraintError;
};


AugmentedLagrangian::AugmentedLagrangian(GenoNLP& genoNLP,
    size_t correctionPairs, bool verbose)
  : _genoNLP(genoNLP),
    _constraintsTol(1e-4),
    _c(1),
    _correctionPairs(correctionPairs),
    _f(0),
    _x(Vector(genoNLP.getN())),
    _augLagG(Vector(genoNLP.getN())),
    _y(Vector::Zero(_genoNLP.getM())),
    _slack(0),
    _maxIter(50),
    _iter(0),
    _funEval(0),
    _inner(0),
    _tol(-1),
    _tolFun(-1),
    _verbose(verbose)
{
}

SolverStatus AugmentedLagrangian::solve()
{
  static const double gamma = 2;
  static const double big = 1E2;
  Scalar rho = 1.0;

  _genoNLP.getStartingPointDual(_y.data());
  AugmentedNLP augmentedNLP(_genoNLP, rho, _y);
  size_t m = augmentedNLP.m();


  /*
  _genoNLP.getStartingPoint(_x.data());
  Vector dummyG(_genoNLP.getN());
  _genoNLP.functionValueAndGradient(_x.data(), _f, dummyG.data());
  double dummyF;
  augmentedNLP.functionValueAndGradient(_x.data(), dummyF, dummyG.data());
  Vector constraintError = augmentedNLP.constraintError();

  for (size_t i = 0; i < m; ++i)
    if (constraintError(mEqualities + i) < 0)
      constraintError(mEqualities + i) = 0;
  std::cout << "constraint Error = " << constraintError.transpose() << std::endl;
  rho = 2 * (std::abs(_f) + 1) / constraintError.squaredNorm();
  if (rho < rhoMin)
    rho = rhoMin;
  if(rho > rhoMax)
    rho = rhoMax;
  */

  LBFGSB solver(augmentedNLP, _correctionPairs, _verbose);
  if (_tol != -1) solver.setParameter("pgtol", _tol);
  if (_tolFun != -1) solver.setParameter("factr", _tolFun);
  double oldFactor = INF;
  double tau = 0.5;
  double muMax = 1E20;
  double muMin = -1E20;
  SolverStatus status = SOLVED;

  _funEval = 0;
  _inner = 0;
  int consecutiveRhoIncreases = 0;
  for (_iter = 0; _iter < _maxIter; ++_iter)
  {
    if (_verbose)
    {
      std::cout << "\niteration " << _iter << std::endl;
      std::cout << "rho = " << rho << std::endl;
    }

    augmentedNLP.setParams(rho, _y);
    solver.restart();
    status = solver.solve();
    _funEval += solver.funEval();
    _inner += solver.iter();

    _x = Vector::Map(solver.x(), _x.size());
    _augLagG = Vector::Map(solver.g(), _x.size());
    if (_verbose)
    {
      std::cout << "problem solved with " << status << std::endl;
      std::cout << "y was " << _y.transpose() << std::endl;
      std::cout << "x = " << _x.transpose() << std::endl;
    }

    augmentedNLP.computeConstraintError(_x.data());
    Vector constraintError = augmentedNLP.constraintError();
    double constraintErrorNorm =
      constraintError.lpNorm<Eigen::Infinity>();

//    std::cout << "constraintError " << constraintErrorNorm << std::endl;
// TODO: check this
    if ((status == SUBOPTIMAL) || (status == UNBOUNDED))
    {
      // check if rho was too small, i.e., constraints are heavily violated
        if (constraintErrorNorm > big)
        {
          std::cout << "rho was most likely too small" << std::endl;
          rho *= 2;
          // reset x
          augmentedNLP.getStartingPoint(_x.data());
          solver.x(_x.data());
          _genoNLP.getStartingPointDual(_y.data());
          continue;
         }
    }
    // TODO: make this more robust
    if ((status != SOLVED) && (status != SUBOPTIMAL))
      break;

    Vector augmentedConstraintError = augmentedNLP.augmentedConstraintError();
    // compute new Lagrange multiplier and safeguard it
    _y = rho * augmentedConstraintError;
    _y = _y.cwiseMax(muMin);
    _y = _y.cwiseMin(muMax);

    //TODO: check this
    double factor = 0; //constraintErrorNorm;

    for (size_t i = 0; i < m; ++i)
    {
      factor = std::max(std::abs(std::max(std::abs(constraintError(i)), -_y(i) / rho)), factor);
    }
//    std::cout << "factor " << factor << std::endl;
    if (factor > tau * oldFactor)
    {
      rho *= gamma;
      consecutiveRhoIncreases++;
    }
    else
        consecutiveRhoIncreases = 0;
    oldFactor = factor;

    // TODO: make this correct

//    std::cout << "y = " << _y.transpose() << std::endl;
//    std::cout << "c error = " << constraintErrorNorm << std::endl;
    if (constraintErrorNorm < _constraintsTol)
    {
      _slack = (_y.array() * constraintError.array()).abs().maxCoeff();
//      std::cout << "c = " << _slack << std::endl;
      if (_slack < _constraintsTol)
        break;
    }

    if (consecutiveRhoIncreases > 20)
    {
      if (_verbose)
        std::cout << "Problem seems infeasible." << std::endl;
      status = INFEASIBLE;
      break;
    }
  }

  Vector dummyG(_genoNLP.getN());
  _genoNLP.functionValueAndGradient(_x.data(), _f, dummyG.data());
  return status;
}

const Scalar* AugmentedLagrangian::x() const
{
  return _x.data();
}

const Scalar* AugmentedLagrangian::y() const
{
  return _y.data();
}

Scalar AugmentedLagrangian::f() const
{
  return _f;
}

const Scalar* AugmentedLagrangian::g() const
{
  return _augLagG.data();
}

Scalar AugmentedLagrangian::slack() const
{
  return _slack;
}

int AugmentedLagrangian::iter() const
{
  return (int) _iter;
}

int AugmentedLagrangian::funEval() const
{
  return (int) _funEval;
}

int AugmentedLagrangian::inner() const
{
  return (int) _inner;
}



bool AugmentedLagrangian::setParameter(std::string parameter, Scalar value)
{
  // lower case parameter
  std::transform(parameter.begin(), parameter.end(),
                 parameter.begin(), ::tolower);
  if (parameter == "pgtol")
  {
    _tol = value;
    return true;
  }

  if (parameter == "factr")
  {
    _tolFun = value;
    return true;
  }

  if (parameter == "constraintstol")
  {
    _constraintsTol = value;
    return true;
  }

  if (parameter == "maxiter")
  {
    _maxIter = (size_t)std::round(value);
    return true;
  }

  if (parameter == "verbose")
  {
    _verbose = value;
    return true;
  }

  return false;
}

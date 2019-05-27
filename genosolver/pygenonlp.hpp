#ifndef __PYGENONLP_HPP__
#define __PYGENONLP_HPP__


#include <functional>
#include "genoNLP.hpp"

class PyGenoNLP : public GenoNLP
{
public:
    PyGenoNLP(Index n, Index m, Scalar* x0,
  		      Scalar* lb, Scalar* ub,
              std::function<void(const Scalar*, Scalar&, Scalar*)> fg,
              std::function<bool(const Scalar*)> callback);

    PyGenoNLP(Index n, Index m, Scalar* x0,
              Scalar* lb, Scalar* ub,
              std::function<void(const Scalar*, Scalar&, Scalar*)> fg,
              std::function<bool(const Scalar*)> callback,
              Scalar* x0_dual,
              Scalar* c_lb, Scalar* c_ub,
              std::function<void(const Scalar*, Scalar*)> cf,
              std::function<void(const Scalar*, const Scalar*,
                                 Scalar*)> cJprod);

    ~PyGenoNLP();

    bool getInfo(Index& n, Index& m) override;
    bool getBounds(Scalar* lb, Scalar* ub) override;
    bool getBoundsConstraints(Scalar* cl, Scalar* cu) override;
    bool getStartingPoint(Scalar* x) override;
    bool getStartingPointDual(Scalar* x_dual) override;
    // TODO
    // add bool
    // if bool is true gradient needs to be computed
    bool functionValueAndGradient(const Scalar* x,
                                  Scalar& functionValue,
                                  Scalar* gradient) override;
    bool functionValueConstraints(const Scalar* x,
                                  Scalar* constValues) override;

    bool gradientConstraintsTimesVector(const Scalar* x,
                                        const Scalar* dualVars,
                                        Scalar* gradValues) override;
    bool callback(const Scalar* x) override;

private:
    Index _n;               // num var
    Index _m;               // num constraints
    Scalar *_lb, *_ub;      // bounds on variables
    Scalar *_x0;            // starting point
    // (x, f_value, g_values)
    std::function<void(const Scalar*, Scalar&, Scalar*)> _fg;
    std::function<bool(const Scalar*)> _callback;
    // (x, c_values)
    Scalar *_x0_dual;       // starting point dual
    Scalar *_c_lb, *_c_ub;  // bounds on constraints
    std::function<void(const Scalar*, Scalar*)> _cf;
    // (x, dual_vars, c_grad_values)
    std::function<void(const Scalar*, const Scalar*, Scalar*)> _cJprod;
    PyGenoNLP(const PyGenoNLP&);
    void operator=(const PyGenoNLP&);
};

#endif

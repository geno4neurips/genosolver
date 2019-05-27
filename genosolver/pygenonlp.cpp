#include <functional>
#include "genoNLP.hpp"
#include "pygenonlp.hpp"

PyGenoNLP::PyGenoNLP(Index n, Index m, Scalar* x0,
            Scalar* lb, Scalar* ub,
            std::function<void(const Scalar*, Scalar&, Scalar*)> fg,
            std::function<bool(const Scalar*)> callback)
  	:_n(n),
  	 _m(0),
  	 _lb(lb),
  	 _ub(ub),
  	 _x0(x0),
  	 _fg(fg),
  	 _callback(callback){}
PyGenoNLP::~PyGenoNLP(){}
PyGenoNLP::PyGenoNLP(Index n, Index m, Scalar* x0,
            Scalar* lb, Scalar* ub,
            std::function<void(const Scalar*, Scalar&, Scalar*)> fg,
            std::function<bool(const Scalar*)> callback,
            Scalar* x0_dual,
            Scalar* c_lb, Scalar* c_ub,
            std::function<void(const Scalar*, Scalar*)> cf,
            std::function<void(const Scalar*, const Scalar*, Scalar*)> cJprod)
    :_n(n),
     _m(m),
     _lb(lb),
     _ub(ub),
     _x0(x0),
     _fg(fg),
     _callback(callback),
     _x0_dual(x0_dual),
     _c_lb(c_lb),
     _c_ub(c_ub),
     _cf(cf),
     _cJprod(cJprod){}


bool PyGenoNLP::getInfo(Index& n, Index& m)
{
    n = _n;
	m = _m;
    return true;
}

bool PyGenoNLP::getBounds(Scalar* lb, Scalar* ub)
{
    for (Index i = 0; i < _n; ++i) {
        lb[i] = _lb[i];
        ub[i] = _ub[i];
    }
    return true;
}


bool PyGenoNLP::getBoundsConstraints(Scalar* cl, Scalar* cu)
{
    for (Index i = 0; i < _m; ++i) {
        cl[i] = _c_lb[i];
        cu[i] = _c_ub[i];
    }
    return true;
};

bool PyGenoNLP::getStartingPoint(Scalar* x)
{
    for (Index i = 0; i < _n; ++i){
        x[i] = _x0[i];
    }
    return true;
}


bool PyGenoNLP::getStartingPointDual(Scalar* x_dual)
{
    for(Index i = 0; i < _m; ++i){
        x_dual[i] = _x0_dual[i];
    }
    return true;
}

bool PyGenoNLP::functionValueAndGradient(const Scalar* x,
			Scalar& functionValue,
			Scalar* gradient)
{
    _fg(x, functionValue, gradient);
    return true;
}

bool PyGenoNLP::functionValueConstraints(const Scalar* x, Scalar* constValues)
{
    _cf(x, constValues);
    return true;
}

bool PyGenoNLP::gradientConstraintsTimesVector(const Scalar* x,
                                               const Scalar* dualVars,
                                               Scalar* gradValues)
{
    _cJprod(x, dualVars, gradValues);
    return true;
}

bool PyGenoNLP::callback(const Scalar* x)
{
    return _callback(x);
}
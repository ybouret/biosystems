#ifndef LI_COMMON_INCLUDED
#define LI_COMMON_INCLUDED 1


#include "y/ios/ocstream.hpp"
#include "y/math/timings.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"
#include "y/math/fit/ls.hpp"

#include "y/os/progress.hpp"

using namespace upsylon;
using namespace math;

typedef Fit::Variables                Variables;
typedef array<double>                 Array;
typedef ODE::DriverCK<double>::Type   ODE_Driver;
typedef ODE::Field<double>::Equation  ODEquation;


struct __lithium
{
    static inline
    double h_ini(const Array &aorg, const Variables &vars)
    {

        return pow(10.0,-vars(aorg,"pH_ini"));
    }


    static inline
    double h_end(const Array &aorg, const Variables &vars)
    {
        return pow(10.0,-vars(aorg,"pH_end"));
    }

    static inline
    double h(const double t, const Array &aorg, const Variables &vars)
    {
        const double  _ini = h_ini(aorg,vars);
        if(t<=0)
        {
            return _ini;
        }
        else
        {
            const double   _end = h_end(aorg,vars);
            return _ini + (_end-_ini) * t/(t+vars(aorg,"t_h"));
        }
    }

    static inline
    double pH(const double t, const Array &aorg, const Variables &vars)
    {
        return -log10( h(t,aorg,vars) );
    }

    static inline
    double eta( double h, const Array &aorg, const Variables &vars )
    {
        const double   pH   = -log10(h);
        const double   u    = vars(aorg,"pw_eta") * ( vars(aorg,"pH_eta") - pH );
        const double   v    = pow(10.0,u);
        return v/(1.0+v);
    }

    static inline
    void save(ios::ostream &fp, const double lt, const Array &Y, const double *extra=NULL)
    {
        fp("%.15g", lt);
        for(size_t i=1;i<=Y.size();++i)
        {
            fp(" %.15g", Y[i]);
        }
        fp(" %.15g", exp(lt));
        if( extra )
        {
            fp(" %.15g", *extra);
        }
        fp << "\n";
    }

};

#endif


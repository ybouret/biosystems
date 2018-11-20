#ifndef FCNLI_INCLUDED
#define FCNLI_INCLUDED 1

#include "y/math/types.hpp"

using namespace upsylon;
using namespace math;

struct Lithium
{
    static inline double BumpSeries(const double lambda, const double u )
    {
        const double fac = u * (lambda-1);
        int    sgn = 1;
        double den = 1;
        double num = fac;

        double ans = u;
        double old_ans = ans;
        for(size_t i=2;;++i)
        {
            sgn = -sgn;
            num *= fac;
            den *= i;
            ans += sgn * num/den;
            if( fabs(ans-old_ans) <= numeric<double>::ftol * fabs(ans) )
            {
                break;
            }
            old_ans = ans;
        }

        return ans * exp(-u);
    }

    static const double ltol;

    static inline double Bump( const double lambda, const double u )
    {

        if(fabs(lambda-1)<ltol)
        {
            return BumpSeries(lambda,u);
        }
        else
        {
            return (exp(-lambda*u)-exp(-u))/(1.0-lambda);
        }
    }

    static inline double UmaxSeries(const double lambda)
    {
        const double fac = lambda-1;
        double num = 1;
        double ans = 1;
        double old_ans = ans;
        int    sgn = 1;
        for(size_t i=2;;++i)
        {
            sgn = -sgn;
            num *= fac;
            ans += (sgn * num)/i;

            if( fabs(ans-old_ans) <= numeric<double>::ftol * fabs(ans) )
            {
                break;
            }
            old_ans = ans;
        }
        return ans;
    }

    static  inline double Umax(const double lambda)
    {
        if(fabs(lambda-1)<ltol)
        {
            return UmaxSeries(lambda);
        }
        else
        {
            return log(lambda)/(lambda-1.0);
        }
    }

    static double Bmax(const double lambda)
    {
        return Bump(lambda,Umax(lambda));
    }

    static inline double Grow(const double tau) { return 1.0-exp(-tau);} 
};

const double Lithium::ltol = 1e-3;

#endif


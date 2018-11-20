#include "y/math/fit/ls.hpp"
#include "y/sequence/vector.hpp"
#include "y/program.hpp"
#include "y/ios/ocstream.hpp"

using namespace upsylon;
using namespace math;

namespace
{
    static double BumpSeries(const double lambda, const double u )
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

    static const double ltol = 1e-3;

    static double Bump( const double lambda, const double u )
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

    static double UmaxSeries(const double lambda)
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

    static  double Umax(const double lambda)
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


    struct BumpApprox
    {
        double Compute( double u, const array<double> &aorg, const Fit::Variables &var)
        {
            const double u0  = var(aorg,"u0");
            const double a  =  var(aorg,"a");
            const double b   = var(aorg,"b");
            const double U   = u-u0;
            return 1.0-1.0/(1.0+exp(-(a*U+b*U*U)));
        }
    };


}
Y_PROGRAM_START()
{
    {
        static const double lam[] = { 0.1, 0.5, 1-ltol*2, 1-ltol/2, 1, 1+ltol/2, 1+ltol*2, 2, 10 };
        static const size_t num = sizeof(lam)/sizeof(lam[0]);
        ios::ocstream fp("bumps.dat");
        {
            for(size_t i=0;i<num;++i)
            {
                const double l = lam[i];
                const double umax = Umax(l);
                const double bmax = Bump(l,umax);
                for(double u=0;u<=10;u += 0.01)
                {
                    fp("%.15g %.15g %.15g\n", u, Bump(lam[i],u), bmax);
                }
                fp << '\n';
            }
        }
    }

    vector<double> L(10000,as_capacity);
    vector<double> B(10000,as_capacity);

    {
        ios::ocstream fp("bmax.dat");
        const double LamMax    = 100;
        const double LogLamMax = log(LamMax);
        const size_t N         = 1000;
        for(size_t i=0;i<=N;++i)
        {
            const double logLam = -LogLamMax + (i<<1)*LogLamMax/N;
            const double lam    = exp(logLam);
            fp("%.15g %.15g %.15g\n", logLam, Bmax(lam), lam );

            L.push_back(logLam);
            B.push_back(Bmax(lam));

        }
    }
    Fit::LeastSquares<double>   ls;
    const size_t                M = L.size();
    vector<double>              F(M);
    Fit::Sample<double>         sample(L,B,F);
    BumpApprox                  approx;
    Fit::Type<double>::Function Bfit( &approx, & BumpApprox::Compute );
    Fit::Variables             &vars = sample.variables;
    vars << "u0" << "a" << "b";
    vector<double> aorg( vars.size() );
    vector<double> aerr( vars.size() );
    vector<bool>   used( vars.size(), true );
    vars(aorg,"u0") = -0.7;
    vars(aorg,"a")  = 0.8;
    vars(aorg,"b")  = 0.0;

    if( !ls.fit(sample, Bfit, aorg, aerr, used ) )
    {
        throw exception("unable to fit Bmax");
    }
    else
    {
        vars.diplay(std::cerr, aorg, aerr);
        ios::ocstream fp("bfit.dat");
        for(size_t i=1;i<=M;++i)
        {
            fp("%g %g %g\n", L[i], B[i], F[i]);
        }
    }

}
Y_PROGRAM_END()


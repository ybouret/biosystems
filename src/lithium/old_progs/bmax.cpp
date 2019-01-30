#include "y/math/fit/ls.hpp"
#include "y/sequence/vector.hpp"
#include "y/program.hpp"
#include "y/ios/ocstream.hpp"
#include "fcnli.hpp"

using namespace upsylon;
using namespace math;

namespace
{

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
        static const double lam[] = { 0.1, 0.5, 1-Lithium::ltol*2, 1-Lithium::ltol/2, 1, 1+Lithium::ltol/2, 1+Lithium::ltol*2, 2, 10 };
        static const size_t num = sizeof(lam)/sizeof(lam[0]);
        ios::ocstream fp("bumps.dat");
        {
            for(size_t i=0;i<num;++i)
            {
                const double l = lam[i];
                const double umax = Lithium::Umax(l);
                const double bmax = Lithium::Bump(l,umax);
                for(double u=0;u<=10;u += 0.01)
                {
                    fp("%.15g %.15g %.15g\n", u, Lithium::Bump(lam[i],u), bmax);
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
            fp("%.15g %.15g %.15g\n", logLam, Lithium::Bmax(lam), lam );

            L.push_back(logLam);
            B.push_back(Lithium::Bmax(lam));

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


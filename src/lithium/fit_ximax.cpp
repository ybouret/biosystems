
#include "yocto/program.hpp"
#include "yocto/math/fit/fit.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/sort/ysort.hpp"
#include "yocto/math/core/tao.hpp"

using namespace yocto;
using namespace math;

static inline double Xi(double u, double p)
{
    return (exp(-p*u)-exp(-u))/(1.0-p);
}

static inline double u_max(double p)
{
    return log(p)/(p-1.0);
}

static inline double XiMax(double p)
{
    return Xi(u_max(p),p);
}

struct Approx
{

    size_t ncall;
    double Compute( double t, const array<double> &a, const Fit::Variables &var)
    {
        const double slope  = a[ var["slope"]  ];
        const double offset = a[ var["offset"] ];
        const double coeff  = a[ var["coeff"]  ];
        const double tt     = t+offset;
        const double xtra   = (tt>0) ? coeff*square_of(tt) : 0;
        return 0.5*(1-tanh(slope*tt+xtra));
    }

};

YOCTO_PROGRAM_START()
{
    int    n     = 8192;
    double width = 5.0;

    const int nn = n+n;
    vector<double> X(nn,as_capacity);
    vector<double> Y(nn,as_capacity);
    vector<double> Z(nn,as_capacity);

    Fit::Sample<double> sample(X,Y,Z);


    size_t N=0;
    for(int i=1;i<=n;++i)
    {
        const double lnp = (i*width)/n;
        X.__push_back(lnp);
        Y.__push_back(XiMax(exp(lnp)));
        Z.__push_back(0);
        ++N;

        if(true)
        {
            X.__push_back(-lnp);
            Y.__push_back(XiMax(exp(-lnp)));
            Z.__push_back(0);
            ++N;
        }
    }
    yCoSort(X,Y,__compare<double>);

    {
        ios::wcstream fp("ximax.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g\n", X[i], Y[i], exp(X[i]));
        }
    }

    Fit::Variables &var = sample.variables;
    var << "slope" << "offset" << "coeff";
    const size_t   nvar = var.size();
    vector<double> aorg(nvar);
    vector<bool>   used(nvar,true);
    vector<double> aerr(nvar);

    aorg[ var["slope"]  ] = 0.3955;
    aorg[ var["offset"] ] = 0.7169;
    aorg[ var["coeff"]  ] = 0.01;

    used[ var["coeff"] ] = true;

    Approx approx;
    Fit::Sample<double>::Function F( & approx, & Approx::Compute );
    Fit::LS<double> fit;

    if( fit.run(sample,F,aorg,used,aerr) )
    {
        sample.display(std::cerr, aorg, aerr, "\t" );
    }
    else
    {
        std::cerr << "couldn't fit" << std::endl;
    }

    {
        ios::wcstream fp("ximax.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g %g\n", X[i], Y[i], exp(X[i]), Z[i]);
        }
    }

}
YOCTO_PROGRAM_END()



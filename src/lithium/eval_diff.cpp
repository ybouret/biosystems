#include "yocto/code/rand32.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/program.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/point3d.hpp"

using namespace yocto;

namespace
{
    static inline
    double compute_variance(const array<double> &u, double &average)
    {
        const size_t n = u.size();
        average = 0;
        for(size_t i=n;i>0;--i)
        {
            average += u[i];
        }
        average/=n;
        double var = 0;
        for(size_t i=n;i>0;--i)
        {
            var += math::Square( u[i] - average );
        }
        return var/n;
    }
}

YOCTO_PROGRAM_START()
{
    size_t       N = 10000;
    rand32_kiss ran;
    ran.wseed();

    const double D = 1.0;
    const double sigma = 2*D;

    // 1D
    {
        vector<double> X(N);
        for(size_t i=1;i<=N;++i)
        {
            X[i] = ran.sym1<double>();
        }
        double ave=0;
        double var=compute_variance(X,ave);
        std::cerr << "1D uniform: average=" << ave << ", variance=" << var << std::endl;

        const double fac = sqrt(sigma);
        for(size_t i=1;;)
        {
            double a,b;
            ran.gaussian(a,b);
            a *= fac;
            b *= fac;
            if(i<=N)
            {
                X[i++] = a;
            }
            else
            {
                break;
            }

            if(i<=N)
            {
                X[i++] = b;
            }
            else
            {
                break;
            }
        }
        double aveG = 0;
        double varG = compute_variance(X,aveG);
        std::cerr << "1D gaussian: average=" << aveG << ", variance=" << varG << std::endl;


    }
}
YOCTO_PROGRAM_END()


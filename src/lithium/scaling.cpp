#include "yocto/program.hpp"
#include "yocto/math/types.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;

YOCTO_PROGRAM_START()
{
    {
        ios::wcstream fp("grow.dat");
        const double dtau = 0.005;
        for(double tau=0.005;tau<=5;tau+=dtau)
        {
            fp("%g %g %g\n", tau, 1.0-exp(-tau), log(tau) );
        }
    }
}
YOCTO_PROGRAM_END()


#include "y/program.hpp"
#include "y/ios/ocstream.hpp"
#include "fcnli.hpp"
#include "y/lua++/state.hpp"
#include "y/type/physics.hpp"

using namespace upsylon;
using namespace math;


#define _INI(NAME) NAME( vm->get<double>(#NAME) )
class Linear
{
public:
    Lua::VM      vm;
    const double Theta;
    const double sigma;
    const double mu7;
    const double mu6;
    Linear( Lua::VM &_vm ) :
    vm(_vm),
    _INI(Theta),
    _INI(sigma),
    _INI(mu7),
    mu6( sigma*mu7 )
    {
        std::cerr << "Theta=" << Theta << std::endl;
        std::cerr << "sigma=" << sigma << std::endl;
    }


    double Beta(double tau, double mu, double theta, double eta, double rho)
    {
        const double c2 = square_of( cos(theta) );
        const double s2 = 1.0 - c2;
        const double t2 = s2/c2;

        return Theta * Lithium::Grow(tau*mu) + s2 * ( Lithium::Grow(tau*mu) ) * t2 / mu;
    }



private:
    Y_DISABLE_COPY_AND_ASSIGN(Linear);
};

Y_PROGRAM_START()
{
    Lua::VM vm = new Lua::State();
    for(int i=1;i<argc;++i)
    {
        vm->doFile(argv[i]);
    }
    Linear lin( vm );

    const string outname = "ratio.dat";
    {
        ios::ocstream fp(outname);

        for(double theta=0;theta<=1.5;theta+=0.1)
        {
            for(double ltau=-5;ltau<=5;ltau += 0.01)
            {
                const double tau = exp(ltau);
                //const double r   = lin.Beta(tau,1.0,theta) / lin.Beta(tau,lin.sigma,theta);
                const double r = 1;
                fp("%g %g %g\n", ltau, r , tau);
            }
            fp("\n");

        }
    }

}
Y_PROGRAM_END()


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
    const double kappa;
    const double fac7;
    const double fac6;

    Linear( Lua::VM &_vm ) :
    vm(_vm),
    _INI(Theta),
    _INI(sigma),
    _INI(mu7),
    mu6( sigma*mu7 ),
    _INI(kappa),
    _INI(fac7),
    fac6( kappa*fac7 )
    {
        std::cerr << "Theta=" << Theta << std::endl;
        std::cerr << "sigma=" << sigma << std::endl;
    }


    double Beta(double tau, double mu, double theta, double fac)
    {
        const double c2  = square_of( cos(theta) );
        const double s2  = 1.0 - c2;
        const double t2  = s2/c2;
        const double lam = 1.0/mu/c2;

        return
        Theta * Lithium::Grow(tau*mu)
        + fac * s2 *  ( Lithium::Grow(tau*mu)  + Lithium::Bump(lam, mu*tau) * t2) / mu;
    }

    double beta7(double tau,double theta)
    {
        return Beta(tau,mu7,theta,fac7);
    }

    double beta6(double tau, double theta)
    {
        return Beta(tau,mu6,theta,fac6);
    }

    double ratio(double tau, double theta)
    {
        return beta7(tau,theta)/beta6(tau,theta);
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

    const string outname = vformat("mu%gkappa%g.dat",lin.mu7,lin.kappa);

    {
        std::cerr << "Saving into " << outname << std::endl;
        ios::ocstream fp(outname);

        for(double theta=0;theta<=1.5;theta+=0.1)
        {
            for(double ltau=-5;ltau<=5;ltau += 0.01)
            {
                const double tau = exp(ltau);
                const double r   = lin.ratio(tau, theta);
                fp("%g %g %g %g\n", ltau, r , lin.beta6(tau,theta), lin.beta7(tau,theta) );
            }
            fp("\n");

        }
    }

}
Y_PROGRAM_END()


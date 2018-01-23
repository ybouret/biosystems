#include "yocto/program.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;

#define NVAR 3
#define EH   1
#define Li6  2
#define Li7  3

typedef numeric<double>::function    Function;
typedef ode::Field<double>::Equation Equation;
typedef ode::Field<double>::Callback Callback;
typedef ode::driverCK<double>::type  ODE_Intg;

class LiSystem
{
public:
    Equation       diffeq;
    ODE_Intg       odeint;
    vector<double> lambda;
    matrix<double> mu;
    double         rho_h;
    double         rho_7;
    explicit LiSystem() :
    diffeq(this, &LiSystem::equations),
    odeint(1e-7),
    lambda(NVAR),
    mu(NVAR,NVAR),
    rho_h(2),
    rho_7(0)
    {
        odeint.start(NVAR);
        lambda[EH]  = 1;
        lambda[Li6] = 1;
        lambda[Li7] = 3;

        mu[1][1] = 1+rho_h;     mu[1][2] = 1.0-rho_7;   mu[1][3] = rho_7;
        mu[2][1] = lambda[Li6]; mu[2][2] = lambda[Li6]; mu[2][3] = 0;
        mu[3][1] = lambda[Li7]; mu[3][2] = 0;           mu[3][3] = lambda[Li7];

    }

    virtual ~LiSystem()
    {
    }

    void equations( array<double> &dZdtau, double tau, const array<double> &Z)
    {
        tao::set(dZdtau,lambda);
        tao::mul_sub(dZdtau,mu,Z);
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(LiSystem);
};

YOCTO_PROGRAM_START()
{
    LiSystem       LS;
    vector<double> Z(NVAR);
    std::cerr << "lam=" << LS.lambda << std::endl;
    std::cerr << "mu="  << LS.mu  << std::endl;

    ios::wcstream fp("output.dat");
    fp("%g %g %g %g\n", 0.0, 0.0, 0.0, 0.0);
    double dtau = 0.001;
    for(double tau=0;tau<=100;tau+=0.01)
    {
        double ctrl = dtau/10;
        LS.odeint(LS.diffeq,Z,tau,tau+dtau,ctrl,NULL);
        fp("%g %g %g %g\n", tau+dtau, Z[1], Z[2], Z[3]);
    }

}
YOCTO_PROGRAM_END()


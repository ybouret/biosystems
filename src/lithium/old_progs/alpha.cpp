
#include "y/program.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"
#include "y/sequence/vector.hpp"
#include "y/math/round.hpp"
#include "y/ios/ocstream.hpp"
#include "y/os/progress.hpp"
#include "y/string/convert.hpp"

using namespace upsylon;
using namespace math;

typedef array<double>                 Array;
typedef ODE::DriverCK<double>::Type   ODE_Driver;
typedef ODE::Field<double>::Equation  ODEquation;

class HatAlpha
{
public:
    static const size_t NVAR = 1;

    double Upsylon;
    double h_ini;
    double h_end;
    double scale;
    double Q;
    double Final;


    HatAlpha() : Upsylon(0.1),
    h_ini(1),
    h_end(1),
    scale(0.5),
    Q(0.5),
    Final((1.0+Q)/(1.0+Q+Upsylon*h_end))
    {
    }

    ~HatAlpha()
    {
    }

    void Compute( Array &dYdt, double tau, const Array &Y)
    {
        const double ac = Y[1];
        const double h  = getH(tau);
        dYdt[1] = 1 - ac*(1+h*Upsylon) + (1-ac)*Q*(1-exp(-scale*tau));
    }

    void Setup( Array &Y)
    {
        Y[1] = 1;
    }

    double getH(double tau)
    {
        return h_ini+(h_end-h_ini)*(1-exp(-scale*tau));
    }
    
private:
    Y_DISABLE_COPY_AND_ASSIGN(HatAlpha);
};

static inline void save( ios::ostream &fp, const double t, const Array &fields, const double r=0)
{
    fp("%.15g",t);
    for(size_t i=1;i<=fields.size();++i)
    {
        fp(" %.15g", fields[i]);
    }
    fp(" %.15g",r);
    fp << '\n';
}


Y_PROGRAM_START()
{
    HatAlpha HA;

    ODEquation      diffeq( &HA, & HatAlpha::Compute );
    ODE_Driver      odeint;
    odeint.start( HA.NVAR );
    vector<double> Y( HA.NVAR );
    odeint.eps = 1e-5;

    Y[1] = 1;

    double t_max   = 10;
    double dt      = 0.001;
    double dt_save = 0.01;
    const size_t every    = simulation_save_every(dt,dt_save);
    const size_t iters    = simulation_iter(t_max,dt,every);
    double h = dt/10;

    const string filename = "alpha.dat";

    std::cerr << "running " << iters << " steps up to " << iters * dt << " time units, saving every " << every << " steps" << std::endl;
    {
        ios::ocstream fp(filename);
        save(fp,0,Y,HA.Final);
    }

    progress bar;

    bar.start();
    for(size_t i=1;i<=iters;++i)
    {
        const double t0 = (i-1)*dt;
        const double t1 =     i*dt;
        odeint( diffeq, Y, t0, t1, h, NULL);
        if( 0 == (i%every) )
        {
            ios::ocstream fp(filename,true);
            save(fp,t1,Y,HA.Final);
            bar.update(i,iters);
            bar.display(std::cerr) << '\r';
        }
    }
    std::cerr << std::endl;
}
Y_PROGRAM_END()


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

static const double beta_s = 12.0192;

class Intake
{
public:
    static const size_t NVAR = 3;
    static const size_t I_AC   = 1;
    static const size_t I_B6   = 2;
    static const size_t I_B7   = 3;

    const double Theta;
    const double d7in;
    const double d7out;
    const double sigma;
    const double eps;
    const double omeps; //!< 1-eps
    const double k7;
    const double k6;
    const double k_eps;
    const double catalytic;
    const double kh;
    const double Ua;
    const double Ub;

    static double getSigma( const double d7In, const double d7Out )
    {
        return (1.0+0.001*d7Out)/(1.0+0.001*d7In);
    }

    inline Intake(const double _Theta,
                  const double _d7in,
                  const double _d7out,
                  const double _catalytic,
                  const double _Ua,
                  const double _Ub) :
    Theta( _Theta ),
    d7in(_d7in),
    d7out(_d7out),
    sigma( getSigma(d7in,d7out) ),
    eps(1.0/(1.0+beta_s*(1.0+0.001*d7out))),
    omeps(1.0-eps),
    k7(1.0),
    k6(sigma*k7),
    k_eps(k6*eps+k7*omeps),
    catalytic(_catalytic),
    kh(catalytic*k7),
    Ua(_Ua),
    Ub(_Ub)
    {
        std::cerr << "Theta   = " << Theta << std::endl;
        std::cerr << "sigma   = " << sigma << std::endl;
        std::cerr << "epsilon = " << eps   << std::endl;
        std::cerr << "k7      = " << k7    << std::endl;
        std::cerr << "k6      = " << k6    << std::endl;
        std::cerr << "k_eps   = " << k_eps << std::endl;
    }

    void setup( Array &Y )
    {
        Y[I_AC] = 1;
        Y[I_B6] = 0;
        Y[I_B7] = 0;
    }

    inline
    double ratio( const Array &Y, const double t ) const
    {
        if(t<=0)
        {
            return 1.0/sigma;
        }
        else
        {
            return Y[I_B7]/Y[I_B6];
        }
    }

    inline
    double delta7( const Array &Y, const double t ) const
    {
        return 1000.0 * ( (1+0.001*d7out) * ratio(Y,t) - 1.0 );
    }

    void compute( Array &dYdt, double, const Array &Y)
    {
        const double ac   = Y[I_AC];
        const double b6   = Y[I_B6];
        const double b7   = Y[I_B7];
        //const double h    = get_h(t);

        dYdt[I_AC]   = kh - ac * (kh+k_eps*Ua);
        dYdt[I_B6]   = k6*( (Theta-b6) + ac * Ub );
        dYdt[I_B7]   = k7*( (Theta-b7) + ac * Ub );
    }

    double get_h(double)
    {
        return pow(1.0,-5.9);
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Intake);
};

static inline void save( ios::ostream &fp, const double t, const Array &fields, const double r=0)
{
    fp("%.15g %.15g",t,r);
    for(size_t i=1;i<=fields.size();++i)
    {
        fp(" %.15g", fields[i]);
    }
    fp << '\n';
}

Y_PROGRAM_START()
{
    double Theta     = 4;
    double catalytic = 2;
    double Ua        = 0;
    double Ub        = 0;
    if(argc>1)
    {
        Theta = string_convert::to<double>(argv[1],"Theta");
    }

    if( argc > 2 )
    {
        catalytic = string_convert::to<double>(argv[2],"catalytic");
    }

    if( argc > 3 )
    {
        Ua = string_convert::to<double>(argv[3],"Ua");
    }

    if( argc > 4 )
    {
        Ub = string_convert::to<double>(argv[4],"Ub");
    }



    const double d7Out = 15.2;
    const double d7In  = 1.05;


    Intake          intake(Theta,d7In,d7Out,catalytic,Ua,Ub);

    ODEquation      diffeq( &intake, & Intake::compute );
    ODE_Driver      odeint;
    odeint.start( Intake::NVAR );
    vector<double> Y( Intake::NVAR );
    odeint.eps = 1e-5;

    intake.setup(Y);
    double t_max   = 10;
    double dt      = 0.001;
    double dt_save = 0.01;
    const size_t every    = simulation_save_every(dt,dt_save);
    const size_t iters    = simulation_iter(t_max,dt,every);
    double h = dt/10;
    const string filename = vformat("Theta%g_catalytic%g_Ua%g_Ub%g.dat", intake.Theta, intake.catalytic, intake.Ua, intake.Ub);
    std::cerr << "running " << iters << " steps up to " << iters * dt << " seconds, saving every " << every << " steps" << std::endl;
    {
        ios::ocstream fp(filename);
        save(fp,0,Y, intake.delta7(Y,0) );
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
            save(fp,t1,Y,intake.delta7(Y,t1));
            bar.update(i,iters);
            bar.display(std::cerr) << '\r';
        }
    }
    std::cerr << std::endl;

}
Y_PROGRAM_END()



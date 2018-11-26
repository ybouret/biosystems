#include "y/program.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"
#include "y/sequence/vector.hpp"
#include "y/math/round.hpp"
#include "y/ios/ocstream.hpp"
#include "y/os/progress.hpp"
#include "y/string/convert.hpp"
#include "y/lua++/state.hpp"

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


#if 0
    static double getSigma( const double d7In, const double d7Out )
    {
        return (1.0+0.001*d7Out)/(1.0+0.001*d7In);
    }
#endif

    const double Theta;
    const double d7out;
    const double eps6;
    const double eps7;
    const double U6;
    const double U7;
    const double Ua;
    const double mu6;
    const double mu7;
    const double eta;
    const double Q6;
    const double Q7;
    const double h_ini;
    const double h_end;
    const double scale;


#define __INI(VAR) VAR( vm.get<double>( #VAR ) )
    inline Intake(Lua::State &vm) :
    __INI(Theta),
    __INI(d7out),
    eps6(1.0/(1.0+beta_s*(1.0+0.001*d7out))),
    eps7(1.0-eps6),
    __INI(U6),
    __INI(U7),
    Ua(eps6*U6+eps7*U7),
    __INI(mu6),
    __INI(mu7),
    __INI(eta),
    __INI(Q6),
    __INI(Q7),
    h_ini(1),
    h_end(1),
    scale(0.5)
    {
        std::cerr << "Theta   = " << Theta << std::endl;
        std::cerr << "eps6    = " << eps6  << std::endl;
        std::cerr << "eps7    = " << eps7  << std::endl;
        std::cerr << "U6      = " << U6    << std::endl;
        std::cerr << "U7      = " << U7    << std::endl;
        std::cerr << "Ua      = " << Ua    << std::endl;
        std::cerr << "mu6     = " << mu6   << std::endl;
        std::cerr << "mu7     = " << mu7   << std::endl;
        std::cerr << "eta     = " << eta   << std::endl;
        std::cerr << "Q6      = " << Q6    << std::endl;
        std::cerr << "Q7      = " << Q7    << std::endl;

    }

    void setup( Array &Y )
    {
        Y[I_AC] = 1;
        Y[I_B6] = 0;
        Y[I_B7] = 0;
    }


    inline
    double ratio( const Array &Y, const double tau )
    {
        if(tau<=0)
        {
            vector<double> dY( NVAR );
            compute(dY,0,Y);
            return dY[I_B7]/dY[I_B6];
        }
        else
        {
            return Y[I_B7]/Y[I_B6];
        }
    }

    inline
    double delta7( const Array &Y, const double t )
    {
        return 1000.0 * ( (1+0.001*d7out) * ratio(Y,t) - 1.0 );
    }

    void compute( Array &dYdt, double tau, const Array &Y)
    {
        const double ac   = Y[I_AC];
        const double b6   = Y[I_B6];
        const double b7   = Y[I_B7];
        const double h    = getH(tau);

        dYdt[I_AC]   = 1.0 -  ac*(1.0+h*Ua) + (1.0-ac)*(eps6*b6*Q6+eps7*b7*Q7);
        dYdt[I_B6]   = mu6*(Theta - b6) + eta*( ac*h*U6 - (1-ac) * Q6 * b6 );
        dYdt[I_B7]   = mu7*(Theta - b7) + eta*( ac*h*U7 - (1-ac) * Q7 * b7 );
    }

    double getH(double tau)
    {
        return h_ini + (h_end-h_ini)*(1.0-exp(-scale*tau));
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
    Lua::State vm;
    for(int i=1;i<argc;++i)
    {
        vm.doFile(argv[i]);
    }


    //const double d7Out = 15.2;
    //const double d7In  = 1.05;


    Intake          intake(vm);

    ODEquation      diffeq( &intake, & Intake::compute );
    ODE_Driver      odeint;
    odeint.start( Intake::NVAR );
    vector<double> Y( Intake::NVAR );
    odeint.eps = 1e-5;

    intake.setup(Y);
    double t_max   = 20;
    double dt      = 0.001;
    double dt_save = 0.01;
    const size_t every    = simulation_save_every(dt,dt_save);
    const size_t iters    = simulation_iter(t_max,dt,every);
    double h = dt/10;
    const string filename = "intake.dat"; //vformat("Theta%g_catalytic%g_Ua%g_Ub%g.dat", intake.Theta, intake.catalytic, intake.Ua, intake.Ub);
    std::cerr << "running " << iters << " steps up to " << iters * dt << " seconds, saving every " << every << " steps" << std::endl;
    {
        ios::ocstream fp(filename);
        //save(fp,0,Y, intake.delta7(Y,0) );
        save(fp,0,Y, intake.ratio(Y, 0) );
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
            //save(fp,t1,Y,intake.delta7(Y,t1));
            bar.update(i,iters);
            bar.display(std::cerr) << '\r';
            save(fp,t1,Y, intake.ratio(Y,t1) );
        }
    }
    std::cerr << std::endl;

}
Y_PROGRAM_END()



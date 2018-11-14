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

class Intake
{
public:
    static const size_t NVAR = 3;
    static const size_t I_AC   = 1;
    static const size_t I_B6   = 2;
    static const size_t I_B7   = 3;

    double k6;
    double k7;
    double Theta;

    inline Intake(const double _Theta) :
    k6(1.1e-2),
    k7(1e-2),
    Theta(_Theta)
    {
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
            return k7/k6;
        }
        else
        {
            return Y[I_B7]/Y[I_B6];
        }
    }

    void compute( Array &dYdt, double t, const Array &Y)
    {
        const double ac   = Y[I_AC];
        const double b6   = Y[I_B6];
        const double b7   = Y[I_B7];
        //const double h    = get_h(t);

        dYdt[I_AC]   = 0;
        dYdt[I_B6]   = k6*(Theta-b6);
        dYdt[I_B7]   = k7*(Theta-b7);
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
    double Theta = 4;
    if(argc>1)
    {
        Theta = string_convert::to<double>(argv[1],"Theta");
    }
    Intake          intake(Theta);
    ODEquation      diffeq( &intake, & Intake::compute );
    ODE_Driver      odeint;
    odeint.start( Intake::NVAR );
    vector<double> Y( Intake::NVAR );
    odeint.eps = 1e-5;

    intake.setup(Y);
    double t_max   = 300;
    double dt      = 0.01;
    double dt_save = 0.1;
    const size_t every    = simulation_save_every(dt,dt_save);
    const size_t iters    = simulation_iter(t_max,dt,every);
    double h = dt/10;
    const string filename = vformat("output_Theta%g.dat", intake.Theta);
    std::cerr << "running " << iters << " steps up to " << iters * dt << " seconds, saving every " << every << " steps" << std::endl;
    {
        ios::ocstream fp(filename);
        save(fp,0,Y, intake.ratio(Y,0) );
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
            save(fp,t1,Y,intake.ratio(Y,t1));
            bar.update(i,iters);
            bar.display(std::cerr) << '\r';
        }
    }
    std::cerr << std::endl;

}
Y_PROGRAM_END()



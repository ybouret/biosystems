
#include "y/program.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"
#include "y/sequence/vector.hpp"
#include "y/lua++/function.hpp"
#include "y/math/kernel/tao.hpp"
#include "y/math/timings.hpp"
#include "y/ios/ocstream.hpp"
#include "y/os/progress.hpp"
#include "fcnli.hpp"

using namespace upsylon;
using namespace math;


typedef array<double>                 Array;
typedef ODE::DriverCK<double>::Type   ODE_Driver;
typedef ODE::Field<double>::Equation  ODEquation;

static const double beta_s = 12.0192;

#define _INI(VAR) VAR( vm->get<double>(#VAR) )
class LiSim
{
public:
    static const size_t I_AC = 1;
    static const size_t I_B6 = 2;
    static const size_t I_B7 = 3;
    static const size_t NVAR = 3;

    Lua::VM      vm;
    const double Theta;
    const double sigma;
    const double d7out;
    const double eps6;
    const double eps7;
    const double d7in;
    const double r0;

    const double mu7;
    const double mu6;


    inline LiSim( const Lua::VM & _vm ) :
    vm( _vm ),
    _INI(Theta),
    _INI(sigma),
    _INI(d7out),
    eps6(1.0/(1.0+beta_s*(1.0+0.001*d7out))),
    eps7(1.0-eps6),
    _INI(d7in),
    r0( (1.0+d7in/1000.0)/(1.0+d7out/1000.0) ),
    _INI(mu7),
    mu6( sigma*mu7 )
    {
        std::cerr << "Theta = " << Theta << std::endl;
        std::cerr << "sigma = " << sigma << std::endl;
        std::cerr << "mu7   = " << mu7   << std::endl;
        std::cerr << "my6   = " << mu6   << std::endl;
        std::cerr << "d7out = " << d7out << std::endl;
        std::cerr << "eps6  = " << eps6  << std::endl;
        std::cerr << "eps7  = " << eps7  << std::endl;
        std::cerr << "d7in  = " << d7in  << std::endl;
#if 0
        std::cerr << "r0    = " << r0    << std::endl;
        std::cerr << "kappa = " << kappa << std::endl;
        std::cerr << "f0    = " << f0    << std::endl;
        std::cerr << "eta   = " << eta   << std::endl;
#endif
    }



    void Compute( array<double> &dY, double, const array<double> &Y )
    {
#if 0
        const double ac = Y[I_AC];
        const double b6 = Y[I_B6];
        const double b7 = Y[I_B7];

        const double h   = 1.0;
        const double ach = ac*h;
        const double aa  = 1.0-ac;

        const double QB6 = Q6*b6;
        const double QB7 = Q7*b7;

        dY[I_AC] = 1.0 - ac * (1.0 + h * Ua ) + aa * (eps6*QB6+eps7*QB7);
        dY[I_B6] = mu6*(Theta-b6) + eta * ( ach * U6 - aa * QB6);
        dY[I_B7] = mu7*(Theta-b7) + eta * ( ach * U7 - aa * QB7);
#endif
    }

    double computeDelta( const array<double> &Y )
    {
        return 1000.0 * ( (1+d7out/1000.0) * Y[I_B7]/Y[I_B6] - 1.0);
    }

    double computeLiAllRatio( const array<double> &Y )
    {
        return eps6 * Y[I_B6] + eps7 * Y[I_B7];
    }


    void setup( array<double> &Y )
    {
        Y[I_AC] = 1;
        Y[I_B6] = 0;
        Y[I_B7] = 0;
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(LiSim);

};

namespace
{
    static inline
    void save( ios::ostream &fp, const double mark, const Array &Y)
    {
        fp("%.15g",mark);
        for(size_t i=1;i<=Y.size();++i)
        {
            fp(" %.15g", Y[i]);
        }
        fp << '\n';
    }
}

Y_PROGRAM_START()
{
    ODE_Driver driver;
    driver.eps = 1e-5;
    Lua::VM vm = new Lua::State;
    for( int i=1; i<argc; ++i )
    {
        vm->doFile(argv[i]);
    }

    LiSim       Li(vm);
    ODEquation  diffeq( &Li, & LiSim::Compute );

    return 0;
    

    const double lt_min = -6;
    double       lt_max =  8;
    double       lt_amp = lt_max-lt_min;
    double       lt_stp = 0.002;
    double       lt_sav = 0.05;
    size_t       every  = 0;
    const size_t iters  = timings::setup(lt_amp, lt_stp, lt_sav, every);
    lt_max = lt_amp + lt_min;
    std::cerr << "#iters=" << iters  << ", saving every " << every << " lt_stp=" << lt_stp << " from " << lt_min << " to " << lt_max << std::endl;

    vector<double> Y( LiSim::NVAR );

    driver.start( Y.size() );
    progress bar;
    bar.start();

    Li.setup(Y);
    {
        ios::ocstream::overwrite("output.dat");
        ios::ocstream::overwrite("li.dat");
    }

    {
        vector<double> dY(Y.size());
        Li.Compute(dY,0,Y);
        std::cerr << "Y0="   << Y << std::endl;
        std::cerr << "dY0=" << dY << std::endl;
        const double r0 = dY[ Li.I_B7 ]/ dY[Li.I_B6];
        std::cerr << "r0=" << r0 << " / " << Li.r0 << std::endl;
    }


    double t0   = 0;
    double ctrl = exp(lt_min)/10;
    for(size_t i=1;i<=iters;++i)
    {
        const double lt1 = lt_min + ( (i-1)*lt_amp )/(iters-1);
        const double t1  = exp(lt1);
        driver( diffeq, Y, t0, t1, ctrl, NULL);
        if(1==i||0==(i%every))
        {
            bar.update(i,iters);
            bar.display(std::cerr) << '\r';
            {
                ios::ocstream fp("output.dat",true);
                save(fp,lt1,Y);
            }
            {
                ios::ocstream fp("li.dat",true);
                fp("%g %g %g\n",lt1,Li.computeDelta(Y),Li.computeLiAllRatio(Y));
            }
        }
        t0 = t1;
    }
    std::cerr << std::endl;
    std::cerr << "Final Y=" << Y << std::endl;
}
Y_PROGRAM_END()

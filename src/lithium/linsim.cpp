#include "y/program.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"
#include "y/sequence/vector.hpp"
#include "y/lua++/function.hpp"
#include "y/math/kernel/tao.hpp"
#include "y/math/timings.hpp"
#include "y/ios/ocstream.hpp"
#include "y/os/progress.hpp"

using namespace upsylon;
using namespace math;


typedef array<double>                 Array;
typedef ODE::DriverCK<double>::Type   ODE_Driver;
typedef ODE::Field<double>::Equation  ODEquation;

static const double beta_s = 12.0192;

#define _INI(VAR) VAR( vm->get<double>(#VAR) )
class LinSim
{
public:
    static const size_t I_AC = 1;
    static const size_t I_B6 = 2;
    static const size_t I_B7 = 3;
    static const size_t NVAR = 3;

    Lua::VM      vm;
    const double Theta;
    const double sigma;
    const double theta;
    const double c2;
    const double s2;
    const double t2;
    const double d7out;
    const double eps6;
    const double eps7;
    Lua::Function<double> get_h;
    double       Ua;

    inline LinSim(const Lua::VM &_vm) :
    vm( _vm ),
    _INI(Theta),
    _INI(sigma),
    _INI(theta),
    c2( square_of( cos(theta) ) ),
    s2( 1.0-c2 ),
    t2( s2/c2 ),
    _INI(d7out),
    eps6(1.0/(1.0+beta_s*(1.0+0.001*d7out))),
    eps7(1.0-eps6),
    get_h("h",vm),
    Ua(0)
    {
        std::cerr << "Theta=" << Theta << std::endl;
        std::cerr << "sigma=" << sigma << std::endl;
        std::cerr << "theta=" << theta << " => c2=" << c2 << std::endl;
        std::cerr << "d7out=" << d7out << " => eps6=" << eps6 << ", eps7=" << eps7 << std::endl;
        const double h0 = get_h(0);
        Ua = t2/h0;
    }

    inline ~LinSim() throw()
    {
    }

    inline void Compute(Array       &dY,
                        const double tau,
                        const Array &Y)
    {
        const double ac = Y[I_AC];
        tao::ld(dY,0);
        const double h = get_h(tau);
        dY[I_AC] = 1.0 - ac * (1.0+Ua*h);
    }

    inline void setup(Array &Y)
    {
        tao::ld(Y,0);
        Y[I_AC] = 1;
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(LinSim);
};

namespace
{
    static inline
    void save( ios::ostream &fp, const double tau, const Array &Y )
    {
        fp("%.15g",tau);
        for(size_t i=1;i<=Y.size();++i)
        {
            fp(" %.15g\n", Y[i]);
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

    LinSim      lin(vm);
    ODEquation  diffeq( &lin, & LinSim::Compute );
    double t_run = 5;
    double dt    = 0.001;
    double emit  = 0.01;
    size_t every = 0;
    const size_t iters = timings::setup(t_run, dt, emit, every);
    std::cerr << "#iters=" << iters << " for t_run=" << t_run << " at dt=" << dt << ", saving every " << emit << std::endl;
    vector<double> Y( LinSim::NVAR );

    driver.start( Y.size() );
    progress bar;
    bar.start();

    lin.setup(Y);
    {
        ios::ocstream fp("output.dat");
        save(fp,0,Y);
    }
    double ctrl = dt/10;
    for(size_t i=1;i<=iters;++i)
    {
        const double t0  = (i-1) * dt;
        const double t1  = i     * dt;
        driver( diffeq, Y, t0, t1, ctrl, NULL);
        if(0==(i%every))
        {
            bar.update(i,iters);
            bar.display(std::cerr) << '\r';
            ios::ocstream fp("output.dat",true);
            save(fp,t1,Y);
        }
    }
    std::cerr << std::endl;


}
Y_PROGRAM_END()



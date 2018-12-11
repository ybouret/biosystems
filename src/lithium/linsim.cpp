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
    const double          h0;
    double       Ua;
    const double mu7;
    const double mu6;
    const double kappa;
    const double phi7;
    const double phi6;
    const double tau_shift;

    inline LinSim(const Lua::VM &_vm) :
    vm( _vm ),
    _INI(Theta),
    _INI(sigma),
    _INI(theta),
    c2( square_of( cos(theta) ) ),
    s2( square_of( sin(theta) ) ),
    t2( square_of( tan(theta) ) ),
    _INI(d7out),
    eps6(1.0/(1.0+beta_s*(1.0+0.001*d7out))),
    eps7(1.0-eps6),
    get_h("h",vm),
    h0( get_h(0) ),
    Ua(t2/h0),
    _INI(mu7),
    mu6( sigma*mu7 ),
    _INI(kappa),
    _INI(phi7),
    phi6( kappa * phi7 ),
    _INI(tau_shift)
    {
        std::cerr << "Theta = " << Theta << std::endl;
        std::cerr << "sigma = " << sigma << std::endl;
        std::cerr << "theta = " << theta << " => c2   = " << c2   << "=> s2 = " << s2 << " => t2 = " << t2 << std::endl;
        std::cerr << "d7out = " << d7out << " => eps6 = " << eps6 << ", eps7=" << eps7 << std::endl;
        std::cerr << "mu7   = " << mu7   << " => mu6  = " << mu6  << std::endl;
        std::cerr << "kappa = " << kappa << std::endl;
        std::cerr << "phi7   = " << phi7   << " => phi6  = " << phi6  << std::endl;

    }

    inline ~LinSim() throw()
    {
    }

    inline void Compute(Array       &dY,
                        const double tau,
                        const Array &Y)
    {
        const double ac = Y[I_AC];
        const double b6 = Y[I_B6];
        const double b7 = Y[I_B7];
        const double h  = get_h(tau);
        const double U  = Ua*h;
        const double V  = U * ac;

        dY[I_AC] = 1.0 - ac * (1.0+U);
        dY[I_B6] = (Theta-b6)*mu6 + phi6 * V;
        dY[I_B7] = (Theta-b7)*mu7 + phi7 * V;

    }

    inline void setup(Array &Y)
    {
        Y[I_AC] = 1;
        Y[I_B6] = 0;
        Y[I_B7] = 0;

        std::cerr << "Initial Ratio0: " << (Theta*mu7+phi7*t2)/(Theta*mu6+phi6*t2) << std::endl;
    }

    inline double alpha0(double tau) const
    {
        return c2 + s2 * exp( -tau/c2 );
    }

    inline
    double Beta(double tau, double mu, double phi)
    {

        const double lam = 1.0/mu/c2;

        return
        Theta * Lithium::Grow(tau*mu)
        + phi * s2 *  ( Lithium::Grow(tau*mu)  + Lithium::Bump(lam, mu*tau) * t2) / mu;
    }

    inline
    double beta7(double tau)
    {
        return Beta(tau,mu7,phi7);
    }

    double beta6(double tau)
    {
        return Beta(tau,mu6,phi6);
    }

    double ratio(double tau)
    {
        return beta7(tau)/beta6(tau);
    }

    inline double computeDelta(const Array &Y) const
    {
        const double r = Y[I_B7]/Y[I_B6];
        return 1000.0 * ( (1+d7out/1000.0) * r - 1.0 );
    }

    inline double computeLiAll(const Array &Y ) const
    {
        const double b6 = Y[I_B6];
        const double b7 = Y[I_B7];
        return (eps6*b6+eps7*b7);
    }


    private:
    Y_DISABLE_COPY_AND_ASSIGN(LinSim);
};

namespace
{
    static inline
    void save( ios::ostream &fp, const double mark, const double r_pred, const double r_real, const double H, const double alpha0, const Array &Y)
    {
        fp("%.15g",mark);
        fp(" %.15g", r_pred);
        fp(" %.15g", r_real);
        fp(" %.15g", H);
        fp(" %.15g", alpha0);
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

    LinSim      lin(vm);
    ODEquation  diffeq( &lin, & LinSim::Compute );



    const double lt_min = -6;
    double       lt_max =  8;
    double       lt_amp = lt_max-lt_min;
    double       lt_stp = 0.002;
    double       lt_sav = 0.05;
    size_t       every  = 0;
    const size_t iters  = timings::setup(lt_amp, lt_stp, lt_sav, every);
    lt_max = lt_amp + lt_min;
    std::cerr << "#iters=" << iters  << ", saving every " << every << " lt_stp=" << lt_stp << " from " << lt_min << " to " << lt_max << std::endl;

    vector<double> Y( LinSim::NVAR );

    driver.start( Y.size() );
    progress bar;
    bar.start();

    lin.setup(Y);
    {
        ios::ocstream::overwrite("output.dat");
        ios::ocstream::overwrite("li.dat");
    }

    {
        vector<double> dY(3);
        lin.Compute(dY, 0, Y);
        std::cerr << "dY=" << dY << std::endl;
        std::cerr << "rho0=" << dY[LinSim::I_B7] / dY[LinSim::I_B6]  << std::endl;
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
            const double r_real = Y[LinSim::I_B7]/Y[LinSim::I_B6];
            {
                ios::ocstream fp("output.dat",true);
                const double r_pred = lin.ratio(t1);
                const double H      = lin.get_h(t1)/lin.h0;
                save(fp,lt1,r_pred,r_real,H,lin.alpha0(t1),Y);
            }
            {
                ios::ocstream fp("li.dat",true);
                const double tt = exp(lt1+lin.tau_shift);
                if(tt<=5000)
                {
                    fp("%.15g %.15g %.15g %.15g\n", lt1, lin.computeDelta(Y), lin.computeLiAll(Y), tt);
                }
            }
        }
        t0 = t1;
    }
    std::cerr << std::endl;
    std::cerr << "Final Ratio: " << Y[LinSim::I_B7]/Y[LinSim::I_B6] << std::endl;
}
Y_PROGRAM_END()



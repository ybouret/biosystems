
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

    const double rho;   //! h_end/h_ini

    const double beta7; //! final beta7
    const double xs7;   //!< beta7/Theta-1

    const double Omega; //! intake parameter
    const double cos2Omega;
    const double tan2Omega;

    const double phi7;
    const double cos2phi7;

    const double phi6;
    const double cos2phi6;

    const double f0;
    const double kappa;

    const double kappa_over_sigma;
    const double xs6;
    const double beta6;
    const double r_end;
    const double d7end;

    const double U7; // h0*Upsilon7
    const double U6; // kappa*U7
    const double Ua; // eps6*U6+eps7*U7

    const double eta;

    const double cos2psi;
    const double sin2psi;
    const double ac_end;  // final
    const double Q6;
    const double Q7;



    double check_r0( const double r0_guess )
    {
        std::cerr << "r0 guess=" << r0_guess << ", sigma=" << sigma << std::endl;
        if(r0_guess>=1/sigma)
        {
            throw exception("invalid r0");
        }
        return r0_guess;
    }

    inline double check_beta7( const double beta7_guess )
    {
        if(beta7_guess<=Theta)
        {
            throw exception("invalid beta7");
        }
        return beta7_guess;
    }

    inline LiSim( const Lua::VM & _vm ) :
    vm( _vm ),
    _INI(Theta),
    _INI(sigma),
    _INI(d7out),
    eps6(1.0/(1.0+beta_s*(1.0+0.001*d7out))),
    eps7(1.0-eps6),
    _INI(d7in),
    r0( check_r0( (1.0+d7in/1000.0)/(1.0+d7out/1000.0)) ),

    _INI(mu7),
    mu6( sigma*mu7 ),

    _INI(rho),

    beta7( check_beta7(vm->get<double>("beta7")) ),
    xs7( beta7/Theta-1.0 ),

    _INI( Omega ),
    cos2Omega( square_of( cos(Omega) ) ),
    tan2Omega( (1.0-cos2Omega)/cos2Omega ),

    _INI(phi7),
    cos2phi7( square_of(cos(phi7)) ),

    _INI(phi6),
    cos2phi6( square_of(cos(phi6)) ),

    f0( rho * cos2Omega * cos2phi7 / xs7 ),
    kappa( (1+(1-r0*sigma)/f0)/r0 ),

    kappa_over_sigma( kappa/sigma ),
    xs6( kappa_over_sigma * cos2phi6/cos2phi7 * xs7 ),
    beta6( Theta * (1+xs6) ),
    r_end( beta7/beta6 ),
    d7end( 1000 * ( (1+0.001*d7out) * r_end -1 ) ),
    U7( tan2Omega/rho/(eps6*kappa+eps7) ),
    U6( kappa*U7 ),
    Ua(eps6*U6+eps7*U7),

    eta( mu7*Theta*f0/U7 ),

    cos2psi( (eps6 * kappa * cos2phi6 + eps7 * cos2phi7 ) / ( eps6*kappa+eps7) ),
    sin2psi( 1.0 - cos2psi ),
    ac_end( cos2Omega * cos2psi + sin2psi ),

    Q6(0),
    Q7(0)
    {
        std::cerr << "Theta = " << Theta << std::endl;
        std::cerr << "sigma = " << sigma << std::endl;
        std::cerr << "mu7   = " << mu7   << std::endl;
        std::cerr << "my6   = " << mu6   << std::endl;
        std::cerr << "d7out = " << d7out << std::endl;
        std::cerr << "eps6  = " << eps6  << std::endl;
        std::cerr << "eps7  = " << eps7  << std::endl;
        std::cerr << "d7in  = " << d7in  << std::endl;
        std::cerr << "r0    = " << r0    << std::endl;
        std::cerr << "rho   = " << rho   << std::endl;
        std::cerr << "beta7 = " << beta7 << ", xs7=" << xs7 << std::endl;
        std::cerr << "Omega = " << Omega << ", cos2=" << cos2Omega << std::endl;
        std::cerr << "phi7  = " << phi7  << ", cos2=" << cos2phi7 << std::endl;
        std::cerr << "phi6  = " << phi6  << ", cos2=" << cos2phi6 << std::endl;

        std::cerr << "f0    = " << f0    << std::endl;
        std::cerr << "kappa = " << kappa << ", kappa/sigma=" << kappa/sigma << std::endl;
        std::cerr << "xs6   = " << xs6 << " (xs7=" << xs7 << ")" << std::endl;
        std::cerr << "beta6 = " << beta6 << std::endl;
        std::cerr << "r_end = " << r_end << " => d7end=" << d7end << "/d7out=" << d7out << std::endl;
        std::cerr << "U7    = " << U7 << std::endl;
        std::cerr << "U6    = " << U6 << std::endl;
        std::cerr << "Ua    = " << Ua << std::endl;

        std::cerr << "eta     = " << eta     << std::endl;
        std::cerr << "cos2psi = " << cos2psi << std::endl;
        std::cerr << "ac_end  = " << ac_end  << std::endl;
    }



    void Compute( array<double> &dY, double, const array<double> &Y )
    {
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

    return 0;

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

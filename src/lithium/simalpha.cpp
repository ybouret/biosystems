#include "y/program.hpp"
#include "y/lua++/state.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"
#include "y/sequence/vector.hpp"
#include "y/ios/ocstream.hpp"
#include "y/os/progress.hpp"
#include "y/math/timings.hpp"

using namespace upsylon;
using namespace math;

typedef array<double>                 Array;
typedef ODE::DriverCK<double>::Type   ODE_Driver;
typedef ODE::Field<double>::Equation  ODEquation;

static double lambda_s = 12.0192;
static double sigma    = 1.0/0.99772;





class Lithium
{
public:
    static const size_t I_AC = 1;
    static const size_t I_B6 = 2;
    static const size_t I_B7 = 3;
    static const size_t NVAR = 3;

    const double d7out;
    const double eps6;
    const double eps7;
    const double Theta;
    const double k7;
    const double k6;
    const double d7ini;
    const double r0;

    const double pH_ini;
    const double pH_end;
    const double h_ini;
    const double h_end;
    const double t_h;
    const double gamma_h;

    const double k0;
    const double mu;
    const double mup;
    const double kappa;
    const double r_mu;

    const double d7end;
    const double r_end;

    const double C2; //!< ac_end
    const double S2; //!< 1-C2;
    const double T2; //!< S2/C2
    const double eta_ini;
    const double eta_end;
    const double scaling;

    const double Lambda;

    Lithium(const double d7out_,
            const double Theta_,
            const double k7_,
            const double d7ini_,
            const double pHini_,
            const double pHend_,
            const double t_h_,
            const double k0_,
            const double mu_,
            const double d7end_,
            const double Lambda_):
    d7out(d7out_),
    eps6(1.0/(1.0+lambda_s*(1.0+1.0e-3*d7out))),
    eps7(1.0-eps6),
    Theta(Theta_),
    k7(k7_),
    k6(sigma*k7),
    d7ini(d7ini_),
    r0( check_r0() ),
    pH_ini(pHini_),
    pH_end(pHend_),
    h_ini( pow(10.0,-pH_ini) ),
    h_end( pow(10.0,-pH_end) ),
    t_h(  t_h_ ),
    gamma_h( compute_gamma_h() ),
    k0(k0_),
    mu(mu_),
    mup( (1.0+mu)/(r0*sigma) - 1.0 ),
    kappa( (1.0+mu)/mu * ( 1.0/r0 - sigma/(1.0+mu) ) ),
    r_mu( compute_r_mu() ),
    d7end( d7end_ ),
    r_end( check_r_end() ),
    C2( compute_C2() ),
    S2( 1.0 - C2 ),
    T2( S2/C2 ),
    eta_ini( get_eta(h_ini) ),
    eta_end( get_eta(h_end) ),
    scaling( (eta_end/eta_ini)*(h_ini/h_end)*T2 ),
    Lambda(Lambda_)
    {
        std::cerr << "d7out = " << d7out  << std::endl;
        std::cerr << "d7ini = " << d7ini  << std::endl;
        std::cerr << "r0    = " << r0     << std::endl;

        if( r0 >= 1.0/sigma )
        {
            throw exception("invalid d7ini");
        }
        std::cerr << "eps6   = " << eps6   << std::endl;
        std::cerr << "eps7   = " << eps7   << std::endl;
        std::cerr << "Theta  = " << Theta  << std::endl;
        std::cerr << "k7     = " << k7     << std::endl;

        std::cerr << "pH_ini = " << pH_ini << std::endl;
        std::cerr << "pH_end = " << pH_end << std::endl;
        std::cerr << "t_h    = " << t_h    << std::endl;
        std::cerr << "mu     = " << mu    << std::endl;
        std::cerr << "kappa  = " << kappa << std::endl;

        std::cerr << "r_mu   = " << r_mu  << std::endl;
        std::cerr << "d7end  = " << d7end << std::endl;
        std::cerr << "r_end  = " << r_end << std::endl;
        std::cerr << "C2     = " << C2    << std::endl;
        std::cerr << "eta_ini= " << eta_ini << std::endl;
        std::cerr << "eta_end= " << eta_end << std::endl;

        std::cerr << "Lambda = " << Lambda << std::endl;
    }


    inline double DeltaOf( const double ratio ) const throw()
    {
        return 1000.0 * ( (1.0+d7out/1000.0) * ratio - 1.0 );

    }

    inline double RatioOf( const double d ) const throw()
    {
        return (1.0+d/1000.0)/(1.0+d7out/1000.0);
    }


    inline double compute_gamma_h() const
    {
        if(h_end>h_ini)
        {
            throw exception("invalid [H] ratio");
        }
        return clamp<double>(0,h_end/h_ini,1);
    }

    inline double check_r0() const
    {
        const double ans = RatioOf(d7ini);
        if(ans >= 1.0/sigma)
        {
            throw exception("d7ini is too high!");
        }
        return ans;
    }

    inline double compute_r_mu() const
    {
        const double sr0 = sigma * r0;
        if(gamma_h>=1.0)
        {
            return sr0;
        }
        else
        {
            return (1.0+mu*gamma_h) / (1.0+ mup * gamma_h );
        }
    }

    inline double check_r_end() const
    {
        const double ans  = RatioOf(d7end);
        if(ans<r_mu)
        {
            throw exception("d7end is too small, should increase mu!");
        }

        if(ans>1.0)
        {
            throw exception("d7end is too high!!!");
        }

        return ans;
    }

    inline double compute_C2() const
    {
        if(r_end>1.0||r_end<r_mu) throw exception("final ratio is invalid, shouldn't happen at this point");
        const double omr = clamp<double>(0,1.0   - r_end,1);
        const double rmr = clamp<double>(0,r_end - r_mu,1);
        const double sr0 = sigma*r0;
        const double cof = ( (gamma_h>=1.0) ? 1.0/sr0 : 1.0 + ( (1.0+mu)/(sigma*r0) - 1.0) * gamma_h );
        return omr / ( omr + cof * rmr );
    }

    void setup( Array &Y )
    {
        Y[I_AC] = 1;
        Y[I_B6] = 0;
        Y[I_B7] = 0;
    }

    inline double get_h( double t )
    {
        if(t<=0)
        {
            return h_ini;
        }
        else
        {
            const double w_end = t / (t+t_h);
            const double w_ini = clamp<double>(0,1.0-w_end,1);
            return h_ini * w_ini + h_end * w_end;
        }
    }


    inline double get_eta( double h ) const
    {
        const double h_eta = 4.0e-7;
        const double p_eta = 1.70;
        const double U = pow(h/h_eta,p_eta);
        return U/(1.0+U);
    }

    void Compute(Array &dY, double t, const Array &Y )
    {
        //const Lithium &self = *this;

        const double ac    = Y[I_AC];
        const double beta6 = Y[I_B6];
        const double beta7 = Y[I_B7];
        const double h     = get_h(t);
        const double eta   = get_eta(h);

        const double phi    = ac * h / h_ini;

        dY[I_AC] = k0 * ( eta * (1.0-ac) - scaling * phi );
        dY[I_B7] = k7 * ( Theta * (1.0+mu*phi)   -   beta7);
        dY[I_B6] = k6 * ( Theta * (1.0+mup*phi)  -   beta6);

    }

    double get_total( const Array &Y ) const throw()
    {
        return Lambda* ( eps6 * Y[I_B6] + eps7 * Y[I_B7] );
    }

    void save(ios::ostream &fp,
              const double  lt,
              Array        &Y,
              const double *extra=0
              ) const
    {
        fp("%.15g",lt);
        for(size_t i=1;i<=Y.size();++i)
        {
            fp(" %.15g", Y[i]);
        }

        fp(" %.15g", get_total(Y));

        if(extra)
        {
            fp(" %.15g",*extra);
        }
        else
        {
            fp << " 0";
        }
        fp << '\n';
    }



};

#include "y/math/fcn/functions.hpp"


static inline Y_LUA_IMPL_CFUNCTION(erf,qerf)

#define INI(NAME) vm->get<double>(#NAME)
Y_PROGRAM_START()
{
    Lua::VM vm = new Lua::State();
    Y_LUA_LOAD_CFUNCTION(vm,erf);

    for( int i=1; i<argc; ++i)
    {
        vm->doFile(argv[i]);
    }


    Lithium  Li(INI(d7out),
                INI(Theta),
                INI(k7),
                INI(d7ini),
                INI(pH_ini),
                INI(pH_end),
                INI(t_h),
                INI(k0),
                INI(mu),
                INI(d7end),
                INI(Lambda));

    ODEquation  diffeq( &Li, & Lithium::Compute );
    ODE_Driver  driver;
    driver.eps = 1e-5;



    const double lt_min = 0;
    double       lt_max =  log(1000);
    double       lt_amp = lt_max-lt_min;
    double       lt_stp = 0.002;
    double       lt_sav = 0.05;
    size_t       every  = 0;
    const size_t iters  = timings::setup(lt_amp, lt_stp, lt_sav, every);
    lt_max = lt_amp + lt_min;
    std::cerr << "#iters=" << iters  << ", saving every " << every << " lt_stp=" << lt_stp << " from " << lt_min << " to " << lt_max << std::endl;

    vector<double> Y( Lithium::NVAR,0 );

    driver.start( Y.size() );
    progress bar;
    bar.start();

    vector<double> dY(Y.size());

    const string sim_name = "output.dat"; //vformat("output_mu%g_k%g_C%g.dat",Li.mu,Li.k0,Li.C2);

    std::cerr << "<saving into " << sim_name << ">" << std::endl;
    ios::ocstream::overwrite(sim_name);


    double t0   = 0;
    double ctrl = exp(lt_min)/1000;
    Li.setup(Y);
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
                ios::ocstream fp(sim_name,true);
                const double  d = Li.DeltaOf( Y[Lithium::I_B7] / Y[Lithium::I_B6] );
                Li.save(fp,lt1,Y,&d);
            }

        }
        t0 = t1;
    }
    std::cerr << std::endl;
    std::cerr << "<saved  into " << sim_name << ">" << std::endl;

}
Y_PROGRAM_END()


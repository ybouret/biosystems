#include "y/program.hpp"
#include "y/lua++/state.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"
#include "y/sequence/vector.hpp"
#include "y/ios/ocstream.hpp"
#include "y/os/progress.hpp"
#include "y/math/timings.hpp"
#include "y/math/fit/ls.hpp"
#include "y/math/fit/samples-io.hpp"


using namespace upsylon;
using namespace math;

////////////////////////////////////////////////////////////////////////////////
//
//
// types definitions
//
//
////////////////////////////////////////////////////////////////////////////////

typedef array<double>                 Array;
typedef ODE::DriverCK<double>::Type   ODE_Driver;
typedef ODE::Field<double>::Equation  ODEquation;

typedef Fit::Sample<double>       Sample;
typedef Fit::LeastSquares<double> LSF;
typedef LSF::Function             Function;
typedef Fit::Variables            Variables;


////////////////////////////////////////////////////////////////////////////////
//
//
// global definitions
//
//
////////////////////////////////////////////////////////////////////////////////
static double lambda_s = 12.0192;
static double sigma    = 1.0/0.99772;


////////////////////////////////////////////////////////////////////////////////
//
//
// lithium class
//
//
////////////////////////////////////////////////////////////////////////////////



class Lithium
{
public:
    static bool  Verbose;
    static const size_t I_AC = 1;
    static const size_t I_B6 = 2;
    static const size_t I_B7 = 3;
    static const size_t NVAR = 3;

    const double SCALING;
    const double d7out;
    const double eps6;
    const double eps7;
    const double sigmap;
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
    const double s0;
    const double s0m1;
    const double mu;
    const double mup;
    const double kappa;
    const double kappap; //!< eps6*kappa+eps7
    const double r_mu;

    const double d7end;
    const double r_end;

    const double C2; //!< ac_end
    const double S2; //!< 1-C2;
    const double T2; //!< S2/C2
    const double eta_ini;
    const double eta_end;
    const double scaling;

    const double Lambda; //!< total external lithium

    Lithium(const double SCALING_,
            const double d7out_,
            const double Theta_,
            const double k7_,
            const double d7ini_,
            const double pHini_,
            const double pHend_,
            const double t_h_,
            const double k0_,
            const double s0_,
            const double d7end_,
            const double Lambda_):
    SCALING(SCALING_),
    d7out(d7out_),
    eps6(1.0/(1.0+lambda_s*(1.0+1.0e-3*d7out))),
    eps7(1.0-eps6),
    sigmap( eps6*sigma + eps7 ),
    Theta(Theta_),
    k7(k7_*SCALING),
    k6(sigma*k7),
    d7ini(d7ini_),
    r0( check_r0() ),
    pH_ini(pHini_),
    pH_end(pHend_),
    h_ini( pow(10.0,-pH_ini) ),
    h_end( pow(10.0,-pH_end) ),
    t_h(  t_h_ ),
    gamma_h( compute_gamma_h() ),
    k0(k0_*SCALING),
    s0(s0_),
    s0m1(s0-1),
    mu( (r0*sigmap * s0m1 - eps6*(1-r0*sigma) ) / (eps7*r0+eps6) ),
    mup( (1.0+mu)/(r0*sigma) - 1.0 ),
    kappa( (1.0+mu)/mu * ( 1.0/r0 - sigma/(1.0+mu) ) ),
    kappap( eps6*kappa+eps7 ),
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
        if(Verbose)
        {
            std::cerr << "d7out   = " << d7out  << std::endl;
            std::cerr << "d7ini   = " << d7ini  << std::endl;
            std::cerr << "r0      = " << r0     << std::endl;

            if( r0 >= 1.0/sigma )
            {
                throw exception("invalid d7ini");
            }
            std::cerr << "eps6    = " << eps6   << std::endl;
            std::cerr << "eps7    = " << eps7   << std::endl;
            std::cerr << "sigma   = " << sigma  << std::endl;
            std::cerr << "sigmap  = " << sigmap << std::endl;
            std::cerr << "Theta   = " << Theta  << std::endl;
            std::cerr << "k7      = " << k7     << std::endl;

            std::cerr << "pH_ini  = " << pH_ini << std::endl;
            std::cerr << "pH_end  = " << pH_end << std::endl;
            std::cerr << "t_h     = " << t_h    << std::endl;

            std::cerr << "s0      = " << s0    << std::endl;
            std::cerr << "mu      = " << mu    << std::endl;
            std::cerr << "kappa   = " << kappa << std::endl;
            std::cerr << "kappap  = " << kappap << std::endl;

            std::cerr << "r_mu    = " << r_mu  << std::endl;
            std::cerr << "d7end   = " << d7end << std::endl;
            std::cerr << "r_end   = " << r_end << std::endl;
            std::cerr << "C2      = " << C2    << std::endl;
            std::cerr << "eta_ini = " << eta_ini << std::endl;
            std::cerr << "eta_end = " << eta_end << std::endl;

            std::cerr << "Lambda  = " << Lambda << std::endl;
        }
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

bool Lithium::Verbose = false;

////////////////////////////////////////////////////////////////////////////////
//
//
// perform one run based on lua settings
//
//
////////////////////////////////////////////////////////////////////////////////

#define INI(NAME) vm->get<double>(#NAME)

#define INI_LIST \
INI(SCALING),    \
INI(d7out),      \
INI(Theta),      \
INI(k7),         \
INI(d7ini),      \
INI(pH_ini),     \
INI(pH_end),     \
INI(t_h),        \
INI(k0),         \
INI(s0),         \
INI(d7end),      \
INI(Lambda)

static inline
void perform_simulation( Lua::VM &vm )
{
    // time parameters
    const double lt_min = 0;
    double       lt_max =  log(3*60*60);
    double       lt_amp = lt_max-lt_min;
    double       lt_stp = 0.002;
    double       lt_sav = 0.05;
    size_t       every  = 0;
    const size_t iters  = timings::setup(lt_amp, lt_stp, lt_sav, every);
    lt_max = lt_amp + lt_min;
    std::cerr << "#iters=" << iters  << ", saving every " << every << " lt_stp=" << lt_stp << " from " << lt_min << " to " << lt_max << std::endl;


    // computation
    Lithium  Li(INI_LIST);

    ODEquation  diffeq( &Li, & Lithium::Compute );
    ODE_Driver  driver;
    driver.eps = 1e-5;




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
#if 0
    std::cerr << "plot 'src/lithium/doc/nhe1_delta7_full_15mM_37_v2.txt' u (log($1)):2 w lp, 'output.dat' u 1:6 w lp, 'src/lithium/data/nhe1_intake_15mM.txt' u (log($1)):2 axis x1y2 w lp, 'output.dat' u 1:5 w l axis x1y2" << std::endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////
//
// fit function
//
////////////////////////////////////////////////////////////////////////////////
#undef INI
#define INI(NAME) vars(aorg,#NAME)
class LiFit
{
public:
    ODE_Driver     driver;
    vector<double> Y;
    inline  LiFit() :
    driver(), Y(Lithium::NVAR)
    {
        driver.start(Y.size());
        driver.eps = 1e-6;
    }

    inline ~LiFit() throw() {}

    inline
    double ComputeDelta7( double lt, const Array &aorg, const Variables &vars )
    {
        
        Lithium  Li(INI_LIST);

        ODEquation  diffeq( &Li, & Lithium::Compute );

        Li.setup(Y);

        double ctrl = lt/1000;
        driver( diffeq,Y,0,exp(lt), ctrl, NULL);

        const double r = Y[Lithium::I_B7]/Y[Lithium::I_B6];
        return Li.DeltaOf(r);
    }

    inline
    double ComputeTotal( double lt, const Array &aorg, const Variables &vars )
    {

        Lithium  Li(INI_LIST);

        ODEquation  diffeq( &Li, & Lithium::Compute );

        Li.setup(Y);

        double ctrl = lt/1000;
        driver( diffeq,Y,0,exp(lt), ctrl, NULL);

        return Li.get_total(Y);
    }



    void save_ln(const string &filename, const double tmax, const Array &aorg,const Variables &vars)
    {
        std::cerr << "<saving to '" << filename << "'>" << std::endl;
        const double lt_min = 0;
        double       lt_max =  log(tmax);
        double       lt_amp = lt_max-lt_min;
        double       lt_stp = 0.05;
        double       lt_sav = lt_stp;
        size_t       every  = 0;
        const size_t iters  = timings::setup(lt_amp, lt_stp, lt_sav, every);
        lt_max = lt_amp + lt_min;

        ios::ocstream::overwrite(filename);
        for(size_t i=1;i<=iters;++i)
        {
            const double lt1 = lt_min + ( (i-1)*lt_amp )/(iters-1);
            const double d7  = ComputeDelta7(lt1,aorg,vars);
            const double tot = ComputeTotal(lt1,aorg,vars);
            ios::ocstream::echo(filename, "%.15g %.15g %.15g\n",lt1,d7,tot);
        }
    }



private:
    Y_DISABLE_COPY_AND_ASSIGN(LiFit);
};

#include "y/math/fcn/functions.hpp"

static inline Y_LUA_IMPL_CFUNCTION(erf,qerf)

#if 0
static const double t_dose  = 60.0;
static const double lt_dose = log(t_dose);
#endif

Y_PROGRAM_START()
{
    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // create VM
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    Lua::VM vm = new Lua::State();
    Y_LUA_LOAD_CFUNCTION(vm,erf);

    if(argc<=2)
    {
        throw exception("usage: %s parameters.lua delta7.txt", program);
        //throw exception("need parameters.lua, delta7.txt, intake");
    }
    vm->doFile(argv[1]);

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    //
    // create delta7 sample
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    vector<double> lt;
    vector<double> delta7;
    vector<double> delta7fit;

    const size_t nd    = Fit::IO::Load(argv[2], 1, lt, 2, delta7, delta7fit);
    const double t_max = lt[nd]+30*60;


    for(size_t i=nd;i>0;--i)
    {
        lt[i] = log(lt[i]);
    }
    Sample       delta7Sample(lt,delta7,delta7fit);
    Variables   &vars = delta7Sample.variables;


    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // create delta7 variables
    //
    //
    ////////////////////////////////////////////////////////////////////////////
#undef  INI
#define INI(NAME) vars.create_global( #NAME )
    INI_LIST;

    const size_t   nvar = vars.size();
    vector<double> aorg(nvar,0);
    vector<double> aerr(nvar,0);
    vector<bool>   used(nvar,false);

#undef INI
#define INI(NAME) vars(aorg,#NAME) = vm->get<double>(#NAME)

    INI_LIST;
    vars.display(std::cerr, aorg, "\t" );

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // create fit functions and wrapper for rescaling
    //
    //
    ////////////////////////////////////////////////////////////////////////////


    LiFit     FitLithium;
    Function  delta7Fit( &FitLithium, & LiFit::ComputeDelta7); //!< fit delta

    LSF      ls;

    if(false)
    {
        perform_simulation(vm);
    }

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // initialize rescaling
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    int level = 0;
    string savename  = "savefit.dat";


    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Start Fitting
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    std::cerr << "Fitting..." << std::endl;

#define FIT_SESSION() do {\
++level;\
std::cerr << "level " << level << std::endl; \
if( ! ls.fit(delta7Sample, delta7Fit, aorg, aerr, used) ) throw exception("couldn't fit level-%d",level);\
vars.display(std::cerr, aorg, aerr, "\t");\
FitLithium.save_ln(savename,t_max, aorg, vars);\
std::cerr << std::endl; } while(false)


    if(true)
    {
        vars.on(used,"k7");
        FIT_SESSION();
    }

    if(true)
    {
        vars.on(used,"k0");
        FIT_SESSION();
    }
    
    if(true)
    {
        vars.on(used,"d7ini");
        FIT_SESSION();
    }
    
    if(false)
    {
        vars.on(used,"mu");
        FIT_SESSION();
    }

    

#if 0
    std::cerr << "plot 'src/lithium/doc/nhe1_delta7_full_15mM_37_v2.txt' u (log($1)):2 w lp, 'savefit.dat' u 1:2 w l, 'src/lithium/data/nhe1_intake_15mM.txt' u (log($1)):2 axis x1y2 w lp, 'savefit.dat' u 1:3 w l axis x1y2" << std::endl;
#endif

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Data exploitation
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    std::cerr << std::endl;
    std::cerr << "<DATA>" << std::endl;

    // create the fitted lithium simulator...
#undef INI
#define INI(NAME) vars(aorg,#NAME)
    Lithium::Verbose = true;
    Lithium  Li(INI_LIST);
    
    std::cerr << "<DATA/>" << std::endl;
    
}
Y_PROGRAM_END()


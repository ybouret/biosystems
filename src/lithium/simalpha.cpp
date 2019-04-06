#include "li-common.hpp"
#include "y/program.hpp"
#include "y/lua++/state.hpp"

static const double lambda_s = 12.0192;
static const double d7out    = 14.57;
static const double sigma    = 1.0/0.99772;
static const double eps6     = 1.0/(1.0+lambda_s*(1.0+1.0e-3*d7out));
static const double eps7     = 1.0-eps6;




class Lithium
{
public:
    static const size_t I_AC = 1;
    static const size_t I_B6 = 2;
    static const size_t I_B7 = 3;
    static const size_t NVAR = 3;

    Variables      vars;
    vector<double> aorg;


    explicit Lithium( Lua::VM &vm ) :
    vars(),
    aorg()
    {
        vars << "pH_ini" << "pH_end" << "t_h";

        vars << "k0" << "pH_eta" << "pw_eta";

        vars << "Omega";

        vars << "k7" << "Theta";

        aorg.make(vars.size(),0);
        std::cerr << "eps6=" << eps6 << std::endl;
        std::cerr << "eps7=" << eps7 << std::endl;
        
        initialize(vm);
    }


#define _LD(NAME) do { self[#NAME] = vm->get<double>(#NAME);\
std::cerr << #NAME " = " << self[#NAME] << std::endl;       \
} while(false)

    void initialize(Lua::VM &vm)
    {
        Lithium &self = *this;
        _LD(pH_ini);
        _LD(pH_end);
        _LD(t_h);
        _LD(k0);
        _LD(pH_eta);
        _LD(pw_eta);
        _LD(Omega);
        _LD(k7);
        _LD(Theta);
    }

    static inline double DeltaOf( const double ratio ) throw()
    {
        return 1000.0 * ( (1.0+d7out/1000.0) * ratio - 1.0 );

    }

    string make_sim_name() const
    {
        const Lithium &self = *this;
        return vformat("out_k0=%.1f_Omega=%.1f_ini=%.1f_end=%.1f_th=%g.dat",
                       self["k0"],
                       self["Omega"],
                       self["pH_ini"],
                       self["pH_end"],
                       self["t_h"]
                       );
    }


    inline double & operator[](const string &id)
    {
        return vars(aorg,id);
    }

    inline const double & operator[](const string &id) const
    {
        return vars(aorg,id);
    }

    inline double get_eta( const double h ) const
    {
        return __lithium::eta( h, aorg, vars);
    }

    inline double get_h_end() const
    {
        return __lithium::h_end(aorg,vars);
    }

    inline double get_h_ini() const
    {
        return __lithium::h_ini(aorg,vars);
    }

    void save_info(double pH_min, double pH_max) const
    {
        if(pH_min>=pH_max) cswap(pH_max,pH_min);

        {
            ios::ocstream fp("eta.dat");
            const size_t  N = 200;
            const double  h_inf   = pow(10.0,-pH_max);
            const double  eta_inf = get_eta(h_inf);
            for(size_t i=0;i<=N;++i)
            {
                const double pH  = pH_min + ( (pH_max-pH_min) * i)/N;
                const double h   = pow(10.0,-pH);
                const double eta = get_eta(h);
                fp("%g %g %g\n", pH, eta, (eta_inf/eta)*(h/h_inf));
            }
        }

    }


    void setup( Array &Y )
    {
        Y[I_AC] = 1;
        Y[I_B6] = 0;
        Y[I_B7] = 0;
    }

    double get_kappa( ) const
    {
        const Lithium &self = *this;
        const double   _end = __lithium::h_end(aorg,vars);
        return self["k0"] * __lithium::eta( _end,aorg,vars ) / _end * square_of( tan( self["Omega"] ) );
    }

    double get_h( double t ) const
    {
        return __lithium::h(t, aorg, vars);
    }




    void Compute( array<double> &dY, double t, const array<double> &Y )
    {
        const Lithium &self = *this;

        const double ac    = Y[I_AC];
        const double beta6 = Y[I_B6];
        const double beta7 = Y[I_B7];

        const double h     = get_h(t);
        const double eta   = get_eta(h);
        const double kh    = self["k0"] * eta;
        const double kappa = get_kappa();
        const double k7    = self["k7"];
        const double k6    = k7 * sigma;
        const double Theta = self["Theta"];

        const double ach = ac*h;
        dY[I_AC] = kh*(1.0-ac) - ach * kappa;
        dY[I_B6] = k6*(Theta-beta6);
        dY[I_B7] = k7*(Theta-beta7);
    }

    double get_gamma0() const
    {
        const Lithium &self = *this;
        return self["k0"] * get_eta( get_h_ini() );
    }

    double get_master(double t) const
    {
        const Lithium &self  = *this;
        const double   Omega = self["Omega"];
        const double   S2    = square_of( sin(Omega) );
        const double   C2    = square_of( cos(Omega) );
        const double   gam0  = get_gamma0();
        return 1.0 - S2 * (1.0 - exp( -t * gam0 / C2 ) );
    }

    double get_ac_end() const
    {
        const  Lithium &self = *this;
        return square_of( cos( self["Omega"] ) );
    }
};


Y_PROGRAM_START()
{
    Lua::VM vm = new Lua::State();
    for( int i=1; i<argc; ++i)
    {
        vm->doFile(argv[i]);
    }

    Lithium     Li(vm);
    ODEquation  diffeq( &Li, & Lithium::Compute );

    ODE_Driver driver;
    driver.eps = 1e-5;


    Li.save_info(5.8,6.8);

    const double lt_min = -6;
    double       lt_max =  log(60*60);
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

    const string sim_name = Li.make_sim_name();
    const double ac_end   = Li.get_ac_end();

    std::cerr << "ac_end  =" << ac_end << std::endl;
    std::cerr << "sim_name=[" << sim_name << "]" << std::endl;
    std::cerr << "gam0    =" << Li.get_gamma0()  << std::endl;

    ios::ocstream::overwrite(sim_name);


    double t0   = 0;
    double ctrl = exp(lt_min)/10;
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
                const double m = Li.get_master(t1);
                const double r = Y[Lithium::I_B7]/Y[Lithium::I_B6];
                const double d7 = Lithium::DeltaOf(r);
                __lithium::save(fp,lt1,Y,&m,&d7);
            }


        }
        t0 = t1;
    }
    std::cerr << std::endl;

}
Y_PROGRAM_END()


#include "li-common.hpp"
#include "y/program.hpp"

class Lithium
{
public:
    static const size_t I_AC = 1;
    static const size_t I_B6 = 2;
    static const size_t I_B7 = 3;
    static const size_t NVAR = 3;

    Variables      vars;
    vector<double> aorg;


    explicit Lithium() :
    vars(),
    aorg()
    {
        vars << "pH_ini" << "pH_end" << "t_h";

        vars << "k0" << "pH_eta" << "pw_eta";

        vars << "Omega";

        aorg.make(vars.size(),0);
        initialize();
    }



    void initialize()
    {
        Lithium &self = *this;

        self["pH_ini"] = 5.8;
        self["pH_end"] = 6.8;  //self["pH_end"] = self["pH_ini"];
        self["t_h"]    = 120.0;

        self["k0"]     = 0.5;
        self["pH_eta"] = 6.39;
        self["pw_eta"] = 1.70;

        self["Omega"]  = 0.5;
    }

    string make_sim_name() const
    {
        const Lithium &self = *this;
        return vformat("out_k0=%g_Omega=%g_ini=%g_end=%g_th=%g.dat",
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

    void save_info(double pH_min, double pH_max) const
    {
        if(pH_min>=pH_max) cswap(pH_max,pH_min);

        {
            ios::ocstream fp("eta.dat");
            const size_t  N = 200;
            for(size_t i=0;i<=N;++i)
            {
                const double pH = pH_min + ( (pH_max-pH_min) * i)/N;
                fp("%g %g\n", pH, __lithium::eta( pow(10.0,-pH), aorg, vars) );
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
        //const Lithium &self = *this;
        return __lithium::h(t, aorg, vars);
    }

    double get_h_ini() const
    {
        return __lithium::h_ini(aorg,vars);
    }


    double get_gamma0() const
    {
        const Lithium &self = *this;
        return self["k0"] * __lithium::eta(__lithium::h_ini(aorg,vars),aorg,vars);
    }

    void Compute( array<double> &dY, double t, const array<double> &Y )
    {
        const Lithium &self = *this;

        const double ac    = Y[I_AC];
        const double h     =  get_h(t);
        const double eta   = __lithium::eta(h, aorg, vars);
        const double kh    = self["k0"] * eta;
        const double kappa = get_kappa();

        dY[I_AC] = kh*(1.0-ac) - ac * h * kappa;
        dY[I_B6] = 0;
        dY[I_B7] = 0;
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
    Lithium     Li;
    ODEquation  diffeq( &Li, & Lithium::Compute );

    ODE_Driver driver;
    driver.eps = 1e-5;


    Li.save_info(5,8);


    const double lt_min = -6;
    double       lt_max =  log(30*60);
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
                //const double rho = Y[Lithium::I_AC] / Li.get_master(t1);
                const double phi = Y[Lithium::I_AC] * Li.get_h(t1)/Li.get_h_ini();
                __lithium::save(fp,lt1,Y,&phi);
            }

        }
        t0 = t1;
    }
    std::cerr << std::endl;

}
Y_PROGRAM_END()


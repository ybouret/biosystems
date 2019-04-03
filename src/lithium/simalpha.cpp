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

        //vars << "kappa";

        aorg.make(vars.size(),0);
        initialize();
    }



    void initialize()
    {
        Lithium &self = *this;

        self["pH_ini"] = 5.8;
        self["pH_end"] = 6.8;  self["pH_end"] = self["pH_ini"];
        self["t_h"]    = 30.0;

        self["k0"]     = 0.5;
        self["pH_eta"] = 6.39;
        self["pw_eta"] = 1.70;

        //self["kappa"]  = 1;
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

    void Compute( array<double> &dY, double t, const array<double> &Y )
    {
        const Lithium &self = *this;

        const double ac  = Y[I_AC];
        const double h   = __lithium::h(t, aorg, vars);
        const double eta = __lithium::eta(h, aorg, vars);
        const double kh  = self["k0"] * eta;

        dY[I_AC] = kh*(1.0-ac);
        dY[I_B6] = 0;
        dY[I_B7] = 0;
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
    
    const string sim_name = "output.dat";
    
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
                __lithium::save(fp,lt1,Y);
            }

        }
        t0 = t1;
    }
    std::cerr << std::endl;

}
Y_PROGRAM_END()


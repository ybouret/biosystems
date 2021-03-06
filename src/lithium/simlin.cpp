
#include "y/program.hpp"
#include "li-common.hpp"


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
        
        vars << "kappa";
        
        aorg.make(vars.size(),0);
        initialize();
    }

    string make_sim_name() const
    {
        const  Lithium &self = *this;
        return vformat("output_k0=%g_kappa=%g_ini=%g_end=%g.dat",
                       self["k0"],
                       self["kappa"],
                       self["pH_ini"],
                       self["pH_end"]
                       );
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
        
        self["kappa"]  = 1;
    }
    
    inline double & operator[](const string &id)
    {
        return vars(aorg,id);
    }
    
    inline const double & operator[](const string &id) const
    {
        return vars(aorg,id);
    }

    inline double get_h_ini() const
    {
        const Lithium &self  = *this;
        return pow(10.0,-self["pH_ini"]);
    }


    inline double get_h_end() const
    {
        const Lithium &self  = *this;
        return pow(10.0,-self["pH_end"]);
    }


    inline double get_h(double t) const
    {
        const Lithium &self  = *this;
        const double   h_ini = get_h_ini();
        if(t<=0)
        {
            return h_ini;
        }
        else
        {
            const double   h_end = get_h_end();
            return h_ini + (h_end-h_ini) * t/(t+self["t_h"]);
        }
    }

    inline double get_pH(double t) const
    {
        return -log10( get_h(t) );
    }
    
    inline double get_eta( double h ) const
    {
        const Lithium &self = *this;
        const double   pH   = -log10(h);
        const double   u    = self["pw_eta"] * ( self["pH_eta"] - pH );
        const double   v    = pow(10.0,u);
        return v/(1.0+v);
    }

    inline double get_ac_inf() const
    {
        const Lithium &self  = *this;
        const double f = self["k0"]    * get_eta( get_h_end() );
        const double g = self["kappa"] * get_h_end()/get_h_ini();
        return f/(f+g);
    }
    
    virtual ~Lithium() throw() {}
    
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
    
    void Compute( array<double> &dY, double t, const array<double> &Y )
    {
        const Lithium &self = *this;
        
        const double ac  = Y[I_AC];
        const double h   = get_h(t);
        const double eta = get_eta(h);
        const double kh  = self["k0"] * eta;
        const double rho = h/get_h_ini();
        
        dY[I_AC] = kh*(1.0-ac) - self["kappa"] * rho * ac;
        dY[I_B6] = 0;
        dY[I_B7] = 0;
    }
    
    void setup( Array &Y )
    {
        Y[I_AC] = 1;
        Y[I_B6] = 0;
        Y[I_B7] = 0;
    }

    
private:
    Y_DISABLE_COPY_AND_ASSIGN(Lithium);
    
};

static inline
void save(ios::ostream &fp, const double lt, const Array &Y, const double *extra=NULL)
{
    fp("%.15g", lt);
    for(size_t i=1;i<=Y.size();++i)
    {
        fp(" %.15g", Y[i]);
    }
    fp(" %.15g", exp(lt));
    if( extra )
    {
        fp(" %.15g", *extra);
    }
    fp << "\n";
}

Y_PROGRAM_START()
{
    Lithium    Li;
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
    const double ac_inf = Li.get_ac_inf();
    std::cerr << "saving into [" << sim_name << "]" << std::endl;
    std::cerr << "ac_inf=" << ac_inf << std::endl;

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
                save(fp,lt1,Y,&ac_inf);
            }

        }
        t0 = t1;
    }
    std::cerr << std::endl;
    
}
Y_PROGRAM_END()


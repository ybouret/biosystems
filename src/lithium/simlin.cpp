
#include "y/program.hpp"
#include "y/math/fit/ls.hpp"
#include "y/ios/ocstream.hpp"

using namespace upsylon;
using namespace math;

typedef Fit::Variables Variables;

class Lithium
{
public:
    Variables      vars;
    vector<double> aorg;
    
    explicit Lithium() :
    vars(),
    aorg()
    {
        vars << "pH_ini" << "pH_end" << "t_h";
        vars << "pH_eta" << "pw_eta";
        
        aorg.make(vars.size(),0);
        setup();
    }
    
    void setup()
    {
        Lithium &self = *this;
        self["pH_eta"] = 6.39;
        self["pw_eta"] = 1.70;
    }
    
    inline double & operator[](const string &id)
    {
        return vars(aorg,id);
    }
    
    inline const double & operator[](const string &id) const
    {
        return vars(aorg,id);
    }
    
    inline double get_h(double t) const
    {
        const Lithium &self = *this;
        const double   h_ini = pow(10.0,-self["pH_ini"]);
        const double   h_end = pow(10.0,-self["pH_end"]);
        if(t<=0)
        {
            return h_ini;
        }
        else
        {
            return h_ini + (h_end-h_ini) * t/(t+self["t_h"]);
        }
    }

    inline double get_pH(double t) const
    {
        return -log10( get_h(t) );
    }
    
    inline double eta( double h ) const
    {
        const Lithium &self = *this;
        const double   pH   = -log10(h);
        const double   u    = self["pw_eta"] * ( self["pH_eta"] - pH );
        const double   v    = pow(10.0,u);
        return v/(1.0+v);
    }
    
    virtual ~Lithium() throw() {}
    
    void save_info(double pH_min, double pH_max) const
    {
        //const Lithium &self  = *this;
        if(pH_min>=pH_max)
            cswap(pH_max,pH_min);
        
        {
            ios::ocstream fp("eta.dat");
            const size_t  N = 100;
            for(size_t i=0;i<=N;++i)
            {
                const double pH = pH_min + ( (pH_max-pH_min) * i)/N;
                fp("%g %g\n", pH, eta( pow(10.0,-pH) ) );
            }
        }
        
    }
    

    
private:
    Y_DISABLE_COPY_AND_ASSIGN(Lithium);
    
};


Y_PROGRAM_START()
{
    Lithium Li;
    
    Li["pH_ini"] = 5.8;
    Li["pH_end"] = 6.2;
    Li["t_h"]    = 30.0;
    
    Li.save_info(5,7);
    
}
Y_PROGRAM_END()


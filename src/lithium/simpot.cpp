#include "y/program.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"
#include "y/math/ode/explicit/explode.hpp"
#include "y/ios/ocstream.hpp"
#include "y/math/fit/ls.hpp"
#include "y/math/fit/samples-io.hpp"
#include "y/type/physics.hpp"
#include "y/sort/heap.hpp"

using namespace upsylon;
using namespace math;

typedef array<double>                Array;

typedef ODE::ExplODE<double>         IODE;
typedef IODE::Solver                 Solver;
typedef IODE::ProblemType            ProblemType;
typedef Fit::Variables               Variables;

static double lambda_s = 12.0192;
static double sigma    = 1.0/0.99772;
static double d7out    = 14.57;
static double Lambda   = 15.0;

static double Texp     = 37.0;

static inline double V2Z( const double V )
{
    return (Y_FARADAY*V)/(Y_R*(Y_ZERO+Texp));
}

static inline double V2Theta( const double V )
{
    return exp( -(Y_FARADAY*V)/(Y_R*(Y_ZERO+Texp)) );
}


static inline double Theta2V( const double th )
{
    return -(Y_R*(Y_ZERO+Texp)) * log(th) / Y_FARADAY;
}

class SimPot : public ProblemType
{
public:
    typedef IODE::Embedded<SimPot>::Type Alias;
    static const size_t I_B7   = 1;
    static const size_t I_B6   = 2;
    static const size_t I_V    = 3;

    const double eps6;
    const double eps7;
    double k7;
    double k6;
    double V0;
    double km;
    double ym;

    explicit SimPot() :
    eps6(1.0/(1.0+lambda_s*(1.0+1.0e-3*d7out))),
    eps7(1.0-eps6),
    k7(0), k6(sigma*k7), V0(0), km(0), ym(0)
    {
    }

    void load( const Array &aorg, const Variables &vars )
    {
        k7     = vars( aorg, "k7" ) * 1e-4;
        k6     = sigma * k7;
        V0     = vars(aorg,"V0");
        km     = vars(aorg,"km");
        ym     = vars(aorg,"ym");
    }

    virtual ~SimPot() throw()
    {
    }

    virtual size_t dimension() const throw() { return 3; }

    virtual void   setup( Array &Y ) const throw()
    {
        Y[I_B6] = 0;
        Y[I_B7] = 0;
        Y[I_V ] = V0;
    }

    virtual double begin() const throw() { return 0;   }

    virtual double delta() const throw() { return 1.0; }

    double  beta_of(const Array &Y ) const
    {
        return (eps6*Y[I_B6]+eps7*Y[I_B7]);
    }

    double  LiTot(const Array &Y ) const
    {
        return Lambda *  beta_of(Y);
    }

    virtual void compute( Array &dYdt, double, const Array &Y )
    {
        const double beta7 = Y[I_B7];
        const double beta6 = Y[I_B6];
        const double V     = Y[I_V];
        const double Theta = V2Theta(V);
        const double beta  = eps6*beta6 + eps7*beta7;

        const double db7 = dYdt[I_B7] = k7*(Theta - beta7);
        const double db6 = dYdt[I_B6] = k6*(Theta - beta6);
        const double dbeta = eps6*db6+eps7*db7;
        dYdt[I_V]  = km*(V0-V)+ym*dbeta;
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(SimPot);
};


class Leak : public SimPot::Alias
{
public:
    IODE iode;

    explicit Leak(const Solver &solver) :
    SimPot::Alias(),
    iode( solver, pointer )
    {
    }

    virtual ~Leak() throw()
    {
    }

    double Compute( double t, const Array &aorg, const Variables &vars )
    {
        SimPot &self = **this;
        self.load(aorg,vars);
        return self.LiTot( iode.at(t) );
    }


    void save(const double tmax)
    {
        SimPot &self = **this;
        iode.reset();

        {
            ios::ocstream fp("beta.dat");
            for(double x=0;x<=tmax; x += 5 )
            {
                const Array &Y = iode.update(x);
                fp("%g %g %g\n", x, self.LiTot(Y), Y[self.I_V]*1000.0 );
            }
        }
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Leak);
};


#define FIT_SESSION() do {\
++level; std::cerr << "Fit Session #" << level << std::endl;\
vars.display(std::cerr,used,"\tusing ");\
if(!ls.fit(sample, F, aorg, aerr, used) ) throw exception("error in fit level-%d", level);\
vars.display(std::cerr, aorg, aerr, "\t");\
leak.save(t_sim);\
} while(false)

Y_PROGRAM_START()
{
    if(argc<=1)
    {
        throw exception("usage: %s data",program);
    }

    vector<double> t;
    vector<double> C;
    vector<double> Cf;

    Fit::Sample<double> sample(t,C,Cf);
    Fit::IO::Load(argv[1], 1,t, 2,C, Cf);
    const size_t N = t.size();
    std::cerr << "Loaded #" << N << std::endl;
    hsort(t,C,comparison::increasing<double>);

    const double t_sim = 4*3600;

    Variables & vars = sample.variables;
    vars << "k7" << "V0" << "km" << "ym";

    Solver        solver   = ODE::DriverCK<double>::New();
    solver->eps = 1e-5;
    Leak          leak(solver);

    Fit::Type<double>::Function F( &leak, & Leak::Compute );


    const size_t   nvar = vars.size();
    vector<double> aorg( nvar, 0 );
    vector<double> aerr( nvar, 0 );
    vector<bool>   used( nvar, false );

    vars(aorg,"k7") = 3;
    vars(aorg,"V0") = -40e-3;
    vars(aorg,"km") = 0.001;
    vars(aorg,"ym") = 0.0;

    Fit::LeastSquares<double> ls;


    leak->load(aorg,vars);
    (void)sample.computeD2(F,aorg);


    std::cerr << "Saving points..." << std::endl;
    {
        leak.save(t_sim);
    }


    int level = 0;
    vars.only_on(used,"k7");
    FIT_SESSION();

    vars.only_on(used,"ym");
    FIT_SESSION();


    vars.only_on(used,"k7");
    FIT_SESSION();

    vars.only_on(used,"k7:ym");
    FIT_SESSION();

    vars.only_on(used,"km");
    FIT_SESSION();

    vars.only_on(used,"km:ym");
    FIT_SESSION();

    vars.only_on(used,"k7");
    FIT_SESSION();

    vars.only_on(used,"k7:ym");
    FIT_SESSION();

    vars.only_off(used,"V0");
    FIT_SESSION();


#if 0
    std::cerr << "Fitting k7" << std::endl;
    vars.on(used,"k7:coeff");
    if(!ls.fit(sample, F, aorg, aerr, used) )
    {
        throw exception("Couldn't fit k7");
    }
    vars.display(std::cerr, aorg, aerr, "\t");

    if(false)
    {
        std::cerr << "Fitting k7/u" << std::endl;

        vars.on(used,"k7:u");
        if(!ls.fit(sample, F, aorg, aerr, used) )
        {
            throw exception("Couldn't fit k7:u");
        }
        vars.display(std::cerr, aorg, aerr, "\t");

    }

    std::cerr << "Saving points..." << std::endl;
    {
        {
            ios::ocstream fp("beta.dat");
            for(size_t i=1;i<=N;++i)
            {
                fp("%g %g %g\n", t[i], C[i], F(t[i],aorg,vars) );
            }
        }


    }

    std::cerr << "Saving fit and potential for unit coefficient..." << std::endl;
    {
        ios::ocstream fp("beta_fit.dat");
        leak->load(aorg,vars);
        
        leak.iode.reset();
        
        for(double x=0;x<=2*t[N]; x += 10 )
        {
            const array<double> &Y    = leak.iode.update(x);
            const double         beta = leak->beta_of(Y);
            const double         y    = beta*Lambda*vars(aorg,"coeff");
            const double         Theta = vars(aorg,"Theta0") * exp( - vars(aorg,"u") * beta );
            const double         V     = -1000.0 * Y_R * (37+Y_ZERO) * log(Theta) / Y_FARADAY;
            fp("%g %g %g %g\n",x,y,Theta,V);
        }
    }



    if(false)
    {
        const string cf = "coeff.dat";
        ios::ocstream::overwrite(cf);
        ios::ocstream::echo(cf,"%g %g %g\n", vars(aorg,"coeff"), vars(aorg,"k7"), vars(aorg,"u"));
        for( double c=0.9;c>=0.49; c -= 0.05 )
        {
            std::cerr << "Fitting with coeff=" << c << std::endl;
            vars(aorg,"coeff") = c;
            if(!ls.fit(sample, F, aorg, aerr, used) )
            {
                throw exception("Couldn't fit k7:u");
            }
            vars.display(std::cerr, aorg, aerr, "\t");
            ios::ocstream::echo(cf,"%g %g %g\n", vars(aorg,"coeff"), vars(aorg,"k7"), vars(aorg,"u"));
        }
    }
#endif


}
Y_PROGRAM_END()


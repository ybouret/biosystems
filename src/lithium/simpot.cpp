#include "y/program.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"
#include "y/math/ode/explicit/explode.hpp"
#include "y/ios/ocstream.hpp"
#include "y/math/fit/ls.hpp"
#include "y/math/fit/samples-io.hpp"
#include "y/type/physics.hpp"

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

class SimPot : public ProblemType
{
public:
    typedef IODE::Embedded<SimPot>::Type Alias;
    static const size_t I_B7 = 1;
    static const size_t I_B6 = 2;

    const double eps6;
    const double eps7;
    double k7;
    double k6;
    double Theta0;
    double u;

    explicit SimPot() :
    eps6(1.0/(1.0+lambda_s*(1.0+1.0e-3*d7out))),
    eps7(1.0-eps6),
    k7(1e-2), k6(sigma*k7), Theta0(4.47), u(0)
    {
    }

    void load( const Array &aorg, const Variables &vars )
    {
        k7     = vars( aorg, "k7" ) * 1e-4;
        k6     = sigma * k7;
        Theta0 = vars(aorg,"Theta0");
        u      = vars(aorg,"u");
    }

    virtual ~SimPot() throw()
    {
    }

    virtual size_t dimension() const throw() { return 2; }
    virtual void   setup( Array &Y ) const throw()
    {
        Y[1] = 0;
        Y[2] = 0;
    }

    virtual double begin() const throw() { return 0; }
    virtual double delta() const throw() { return 1.0; }

    double  beta_of(const Array &Y ) const
    {
        return (eps6*Y[I_B6]+eps7*Y[I_B7]);
    }

    virtual double query( const Array &Y, double ) const
    {
        return Lambda * beta_of(Y);
    }

    virtual void compute( Array &dYdt, double, const Array &Y )
    {
        const double beta7 = Y[I_B7];
        const double beta6 = Y[I_B6];
        const double beta  = eps6*beta6 + eps7*beta7;
        const double Theta = Theta0*exp(-u*beta);

        dYdt[I_B7] = k7*(Theta - beta7);
        dYdt[I_B6] = k6*(Theta - beta6);

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
        return iode.at(t);
    }


private:
    Y_DISABLE_COPY_AND_ASSIGN(Leak);
};



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

    Variables & vars = sample.variables;
    vars << "k7" << "Theta0" << "u";

    Solver        solver   = ODE::DriverCK<double>::New();
    solver->eps = 1e-5;
    Leak          leak(solver);

    Fit::Type<double>::Function F( &leak, & Leak::Compute );


    const size_t   nvar = vars.size();
    vector<double> aorg( nvar, 0 );
    vector<double> aerr( nvar, 0 );
    vector<bool>   used( nvar, false );

    vars(aorg,"k7")     = 3;
    vars(aorg,"Theta0") = 4.47;
    vars(aorg,"u")      = 0.7;


    Fit::LeastSquares<double> ls;
    vars.on(used,"k7");
    if(!ls.fit(sample, F, aorg, aerr, used) )
    {
        throw exception("Couldn't fit k7");
    }
    vars.on(used,"k7:u");
    if(!ls.fit(sample, F, aorg, aerr, used) )
    {
        throw exception("Couldn't fit k7:u");
    }

    vars.display(std::cerr, aorg, aerr, "\t");

    {
        ios::ocstream fp("beta.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g\n", t[i], C[i], F(t[i],aorg,vars) );
        }
    }

    {
        ios::ocstream fp("beta_fit.dat");
        leak->load(aorg,vars);

        for(double x=0;x<=2*t[N]; x += 10 )
        {
            const array<double> &Y    = leak.iode.state_at(x);
            const double         y    = leak->query(Y,x);
            const double         beta = leak->beta_of(Y);
            const double         Theta = vars(aorg,"Theta0") * exp( - vars(aorg,"u") * beta );
            const double         V     = -1000.0 * Y_R * (37+Y_ZERO) * log(Theta) / Y_FARADAY;
            fp("%g %g %g %g\n",x,y,Theta,V);
        }
    }




}
Y_PROGRAM_END()


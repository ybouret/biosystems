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


    inline double DeltaOf( const double ratio ) const throw()
    {
        return 1000.0 * ( (1.0+d7out/1000.0) * ratio - 1.0 );

    }

    Lithium( const double d7out_ ) :
    d7out(d7out_),
    eps6(1.0/(1.0+lambda_s*(1.0+1.0e-3*d7out))),
    eps7(1.0-eps6)
    {
        std::cerr << "d7out = " << d7out << std::endl;
        std::cerr << "eps6  = " << eps6  << std::endl;
        std::cerr << "eps7  = " << eps7  << std::endl;
    }


    void setup( Array &Y )
    {
        Y[I_AC] = 1;
        Y[I_B6] = 0;
        Y[I_B7] = 0;
    }




    void Compute(Array &dY, double t, const Array &Y )
    {
        const Lithium &self = *this;

        const double ac    = Y[I_AC];
        const double beta6 = Y[I_B6];
        const double beta7 = Y[I_B7];

        dY.ld(0);
    }

    void save(ios::ostream &fp,
              Array        &Y,
              const double  mark) const
    {
        fp("%.15g",mark);
        for(size_t i=1;i<=Y.size();++i)
        {
            fp(" %.15g", Y[i]);
        }
        fp << '\n';
    }



};


#define INI(NAME) vm->get<double>(#NAME)
Y_PROGRAM_START()
{
    Lua::VM vm = new Lua::State();
    for( int i=1; i<argc; ++i)
    {
        vm->doFile(argv[i]);
    }


    Lithium  Li( INI(d7out) );
    ODEquation  diffeq( &Li, & Lithium::Compute );
    ODE_Driver  driver;
    driver.eps = 1e-5;



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

    const string sim_name = "output.dat";

    ios::ocstream::overwrite(sim_name);


    double t0   = 0;
    double ctrl = exp(lt_min)/10;
    Li.setup(Y);
    for(size_t i=1;i<=iters;++i)
    {
        const double lt1 = lt_min + ( (i-1)*lt_amp )/(iters-1);
        const double t1  = exp(lt1);
        //driver( diffeq, Y, t0, t1, ctrl, NULL);
        if(1==i||0==(i%every))
        {
            bar.update(i,iters);
            bar.display(std::cerr) << '\r';
#if 0
            {
                ios::ocstream fp(sim_name,true);
                const double m = Li.get_master(t1);
                const double r = Y[Lithium::I_B7]/Y[Lithium::I_B6];
                const double d7 = Lithium::DeltaOf(r);
                __lithium::save(fp,lt1,Y,&m,&d7);
            }
#endif

        }
        t0 = t1;
    }
    std::cerr << std::endl;

}
Y_PROGRAM_END()


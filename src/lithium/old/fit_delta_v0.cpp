#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fit/fit.hpp"

using namespace yocto;
using namespace math;

static const double rho_s = 12.0192;
static const double omp_max = 0.01;

static inline
double Xi(const double u, const double p)
{
    const double omp = 1.0-p;
    if(fabs(omp)<0.01)
    {
        const double central = u*exp(-u);
        double       sum     = 1;
        double       mulby   = u*omp;
        size_t       divby   = 2;
        double       factor  = mulby/divby;
        while(true)
        {
            const double new_sum = sum + factor;
            if(new_sum<=sum)
                break;
            factor *= mulby;
            factor /= ++divby;
            factor = -factor;
            sum    = new_sum;
        }
        //std::cerr << "(" << sum << ")";
        return central * sum;
    }
    else
    {
        return (exp( -p*u ) - exp(-u))/omp;
    }
}

static inline
void testXi()
{
    static const double p[] = { 0.01, 0.1, 1-omp_max/2, 1, 1+omp_max/2, 10, 100 };
    ios::wcstream fp("xi_test.dat");
    for(double u=0;u<10;u += 0.01)
    {
        fp("%.15g", u);
        for(size_t j=0;j<sizeof(p)/sizeof(p[0]);++j)
        {
            fp(" %.15g", Xi(u,p[j]) );
        }
        fp("\n");
    }
}

static inline
double Growth(const double u) throw()
{
    return 1.0 - exp(-u);
}


typedef Fit::Samples<double> Samples;
typedef Fit::Sample<double>  Sample;
typedef Fit::Variables       Variables;
typedef vector<double>       Vector;

class DeltaFit
{
public:
    explicit DeltaFit() : t_reponse(60.0) {}
    virtual ~DeltaFit() throw() {}

    double t_reponse;

    double DeltaOfTime( const double t, const array<double> &aorg, const Variables &vars )
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double sigma  = aorg[ vars["sigma"]  ];
        const double psi    = aorg[ vars["psi"]    ];
        const double d7out  = aorg[ vars["d7out"]  ];

        const double tau7    = k7*t;
        const double tau6    = lambda*tau7;
        const double beta7   = Growth(tau7) + psi * Xi(tau7,sigma);
        const double beta6   = Growth(tau6) + psi * Xi(tau6,sigma/lambda);

        return 1000.0 * ( (1.0+d7out/1000.0) * (beta7/beta6) - 1.0 );
    }



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(DeltaFit);
};

YOCTO_PROGRAM_START()
{

    testXi();


    Samples samples;

    //__________________________________________________________________________
    //
    // Loading Data
    //__________________________________________________________________________

    const string workdir = "src/lithium/doc/";

    vector<double> t_full;      //! col 1
    vector<double> delta_full;  //! col 3
    vector<double> delta_full_fit;


    {
        const string filename = workdir + "nhe1_delta7_full_15mM_37.txt";
        std::cerr << "Loading " << filename << std::endl;
        data_set<double> ds;
        ds.use(1,t_full);
        ds.use(3,delta_full);
        ios::icstream fp(filename);
        ds.load(fp);
        std::cerr << "...done" << std::endl;
    }
    delta_full_fit.make( t_full.size() );
    Sample &sample_full = samples.add(t_full,delta_full,delta_full_fit);

    vector<double> LiExtDR;
    vector<double> LiIntDR;
    vector<double> deltaDR;
    vector<double> stddevDR;
    {
        const string filename = workdir + "nhe1_reponse_60s_37.txt";
        std::cerr << "Loading " << filename << std::endl;
        data_set<double> ds;
        ds.use(1,LiExtDR);
        ds.use(2,LiIntDR);
        ds.use(3,deltaDR);
        ds.use(4,stddevDR);
        ios::icstream fp(filename);
        ds.load(fp);
        std::cerr << "...done" << std::endl;
    }

    vector<double> t_short;
    vector<double> Li_short;
    {
        const string filename = workdir + "nhe1_total_short_times_15mM_37.txt";
        std::cerr << "Loading " << filename << std::endl;
        data_set<double> ds;
        ds.use(1,t_short);
        ds.use(2,Li_short);
        ios::icstream fp(filename);
        ds.load(fp);
        std::cerr << "...done" << std::endl;
    }

    //__________________________________________________________________________
    //
    // global variables
    //__________________________________________________________________________
    Variables &vars = samples.variables;
    vars << "k7" << "lambda" << "sigma" << "psi" << "d7out";
    for(Samples::iterator i=samples.begin();i!=samples.end();++i)
    {
        (**i).variables = vars;
    }
    samples.link();

    const size_t nvar = vars.size();
    Vector       aorg(nvar);
    Vector       aerr(nvar);
    vector<bool> used(nvar,true);

    used[ vars["d7out"] ] = false;
    aorg[ vars["d7out"] ] = 14.98;



    Fit::LS<double> lsf;

}
YOCTO_PROGRAM_END()


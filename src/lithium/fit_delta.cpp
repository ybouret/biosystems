#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fit/fit.hpp"
#include "yocto/container/utils.hpp"
#include "yocto/sort/quick.hpp"
using namespace yocto;
using namespace math;

//static const double rho_s = 12.0192;
static const double omp_max = 0.001;

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
    double t_short;
    double t_long;

    explicit DeltaFit() :
    t_short(0.0),
    t_long(0.0)
    {}

    virtual ~DeltaFit() throw() {}

    double t_reponse;


    double Compute(const double t, const array<double> &aorg, const Variables &vars )
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];
        const double sigma  = aorg[ vars["sigma"]  ];
        const double tau7    = k7*t;
        const double tau6    = lambda*tau7;
        if(t>=t_long)
        {
            return 1000.0 * ( (1.0+d7out/1000.0) * (Growth(tau7)/Growth(tau6)) - 1.0 );
        }

        if(t<=t_short)
        {
            return 1000.0 * ( (1.0+d7out/1000.0) * (Xi(tau7,sigma)/Xi(tau6,sigma/lambda)) - 1.0 );
        }

        return 0;
    }


    double d7ini( const array<double> &aorg, const Variables &vars )
    {
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];
        return (d7out + 1000.0 * (1.0-lambda))/lambda;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(DeltaFit);
};

static inline void save_data(ios::ocstream       &fp,
                             const array<double> &t,
                             const array<double> &f,
                             const array<double> &g,
                             const array<double> *h=NULL)
{
    for(size_t i=1;i<=t.size();++i)
    {
        fp("%.15g %.15g %.15g", log(t[i]), f[i], g[i] );
        if(h)
        {
            fp(" %.15g", (*h)[i] );
        }
        fp << "\n";
    }
}

YOCTO_PROGRAM_START()
{
    DeltaFit        dfn;

    if(false)
    {
        testXi();
    }

    if(argc<=2)
    {
        throw exception("usage: %s t_short t_long",program);
    }
    double &t_short = dfn.t_short;
    double &t_long  = dfn.t_long;
    t_short = strconv::to<double>(argv[1],"t_short");
    t_long  = strconv::to<double>(argv[2],"t_long");


    //__________________________________________________________________________
    //
    // Loading Data
    //__________________________________________________________________________

    const string workdir = "src/lithium/doc/";

    Vector t;      //! col 1
    Vector delta;  //! col 3


    {
        const string filename = workdir + "nhe1_delta7_full_15mM_37.txt";
        std::cerr << "Loading " << filename << std::endl;
        data_set<double> ds;
        ds.use(1,t);
        ds.use(3,delta);
        ios::icstream fp(filename);
        ds.load(fp);
        std::cerr << "...done" << std::endl;
        co_qsort(t,delta,__compare<double>);
    }

    const size_t   N = t.size();
    Vector         deltaFit(N);
    Vector         ratio(N);
    Sample         sample(t,delta,deltaFit);

    //__________________________________________________________________________
    //
    // create sub vectors
    //__________________________________________________________________________
    Vector tIni(N,as_capacity), deltaIni(N,as_capacity);
    Vector tEnd(N,as_capacity), deltaEnd(N,as_capacity);
    for(size_t i=1;i<=N;++i)
    {
        if(t[i]<=t_short)
        {
            tIni.push_back(t[i]);
            deltaIni.push_back(delta[i]);
        }
        if(t[i]>=t_long)
        {
            tEnd.push_back(t[i]);
            deltaEnd.push_back(delta[i]);
        }
    }

    {
        ios::wcstream fp("delta0.dat");
        save_data(fp,t,delta,delta);
    }
    //__________________________________________________________________________
    //
    // create sub samples
    //__________________________________________________________________________
    Fit::Samples<double> multiple(3);

    const size_t N_Ini = tIni.size();
    Vector       deltaIniFit(N_Ini);
    Fit::Sample<double> &sampleIni = multiple.add(tIni,deltaIni,deltaIniFit);

    const size_t N_End = tEnd.size();
    Vector       deltaEndFit(N_End);
    Fit::Sample<double> &sampleEnd = multiple.add(tEnd,deltaEnd,deltaEndFit);

    //__________________________________________________________________________
    //
    // global variables
    //__________________________________________________________________________
    Fit::Variables &vars = multiple.variables;
    vars << "k7" << "lambda" << "d7out" << "sigma";
    sampleIni.variables = vars;
    sampleEnd.variables = vars;


    const size_t nvar = vars.size();
    Vector       aorg(nvar);
    Vector       aerr(nvar);
    vector<bool> used(nvar,false);
    std::cerr << "nvar=" << nvar << std::endl;
    double &k7     = aorg[ vars["k7"]     ];
    double &lambda = aorg[ vars["lambda"] ];
    double &d7out  = aorg[ vars["d7out"]  ];
    double &sigma  = aorg[ vars["sigma"]  ];

    const double dmin = find_min_of(delta);

    k7     = 0.003;
    d7out  = 15.00;
    lambda = (1000.0+d7out)/(1000.0+dmin);
    sigma  = 10.0;

    Fit::LS<double> lsf;
    Fit::Type<double>::Function F(  &dfn, & DeltaFit::Compute  );

    used[ vars["k7"]     ] = true;
    used[ vars["lambda"] ] = false;
    used[ vars["d7out"]  ] = false;
    used[ vars["sigma"]  ] = false;

    //__________________________________________________________________________
    //
    // fit the end to get time scale
    //__________________________________________________________________________
    std::cerr << "Fit relaxing scale..." << std::endl;
    if(!lsf.run(sampleEnd,F,aorg,used,aerr))
    {
        throw exception("cannnot fit end");
    }
    sampleEnd.display(std::cerr,aorg,aerr);
    {
        ios::wcstream fp("delta_end.dat");
        save_data(fp,tEnd,deltaEnd,deltaEndFit);
    }

    used[ vars["k7"]     ] = false;
    used[ vars["sigma"]  ] = true;
    
    if(!lsf.run(sampleIni,F,aorg,used,aerr))
    {
        throw exception("cannnot fit ini");
    }
    sampleEnd.display(std::cerr,aorg,aerr);
    {
        ios::wcstream fp("delta_ini.dat");
        save_data(fp,tIni,deltaIni,deltaIniFit);
    }


}
YOCTO_PROGRAM_END()


#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fit/fit.hpp"
#include "yocto/container/utils.hpp"
#include "yocto/sort/quick.hpp"
#include "yocto/math/core/tao.hpp"
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

//static const bool USE_LOG = true; //! USE_LOG => t is log(t) at load time!

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


        if(t<=t_short)
        {
            return 1000.0 * ( (1.0+d7out/1000.0) * (Xi(tau7,sigma)/Xi(tau6,sigma/lambda)) - 1.0 );
        }

        if(t>=t_long)
        {
            return 1000.0 * ( (1.0+d7out/1000.0) * (Growth(tau7)/Growth(tau6)) - 1.0 );
        }

        return 0;
    }

    double ComputeIni(const double t, const array<double> &aorg, const Variables &vars )
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];
        const double sigma  = aorg[ vars["sigma"]  ];
        const double tau7    = k7*t;
        const double tau6    = lambda*tau7;

        return 1000.0 * ( (1.0+d7out/1000.0) * (Xi(tau7,sigma)/Xi(tau6,sigma/lambda)) - 1.0 );
    }


    double ComputeEnd(const double t, const array<double> &aorg, const Variables &vars )
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];
        const double sigma  = aorg[ vars["sigma"]  ];
        const double tau7    = k7*t;
        const double tau6    = lambda*tau7;

        return 1000.0 * ( (1.0+d7out/1000.0) * (Growth(tau7)/Growth(tau6)) - 1.0 );
    }

    double ComputeFull(const double t, const array<double> &aorg, const Variables &vars )
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];
        const double sigma  = aorg[ vars["sigma"]  ];

        const double tau7    = k7*t;
        const double tau6    = lambda*tau7;

        const double lo = clamp<double>(0,CouplingRamp(t,aorg,vars),1);
        const double hi = 1.0 - lo;

#if 1
        const double beta7    = lo * Xi(tau7,sigma)        + hi * Growth(tau7);
        const double beta6    = lo * Xi(tau6,sigma/lambda) + hi * Growth(tau6);


        return 1000.0 * ( (1.0+d7out/1000.0) * (beta7/beta6) - 1.0 );
#else

        const double ratio_lo = Xi(tau7,sigma)/Xi(tau6,sigma/lambda);
        const double ratio_hi = Growth(tau7)/Growth(tau6);
        const double ratio    = lo * ratio_lo + hi * ratio_hi;
        return 1000.0 * ( (1.0+d7out/1000.0) * ratio - 1.0 );
#endif

    }

    double d7ini( const array<double> &aorg, const Variables &vars )
    {
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];
        return (d7out + 1000.0 * (1.0-lambda))/lambda;
    }

    double CouplingRamp(const double t, const array<double> &aorg, const Variables &vars )
    {
        const double t0    = aorg[ vars[ "t0" ] ];
        const double scale = aorg[ vars[ "scale" ] ];
#if 0
        if(t<=t0)
        {
            return 1.0;
        }
        else if( t<=t0+scale )
        {
            return 1.0-(t-t0)/scale;
        }
        else
        {
            return 0.0;
        }
#endif
        return 0.5*(1.0-tanh( (t-t0)/scale ) );
    }

    double ComputePsi(const double t,
                      const double delta7,
                      const array<double> &aorg,
                      const Variables &vars )
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double sigma  = aorg[ vars["sigma"]  ];
        const double d7out  = aorg[ vars["d7out"]  ];

        const double tau7    = k7*t;
        const double tau6    = lambda*tau7;

        const double Xi7 = Xi(tau7,sigma);
        const double Xi6 = Xi(tau6,sigma/lambda);
        const double G7  = Growth(tau7);
        const double G6  = Growth(tau6);

        const double rho = (1.0+delta7/1000.0)/(1.0+d7out/1000.0);
        const double num = Xi7 - rho * Xi6;
        const double den = rho*G6 - G7;

        return num/den;
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
        //const double tt = USE_LOG ? t[i] : log(t[i]);
        const double tt = log(t[i]);
        fp("%.15g %.15g %.15g", tt, f[i], g[i] );
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
    vars << "k7" << "lambda" << "d7out" << "sigma" << "t0" << "scale";
    sampleIni.variables = vars;
    sampleEnd.variables = vars;
    sample.variables    = vars;


    const size_t nvar = vars.size();
    Vector       aorg(nvar);
    Vector       aerr(nvar);
    vector<bool> used(nvar,false);
    std::cerr << "nvar=" << nvar << std::endl;
    double &k7     = aorg[ vars["k7"]     ];
    double &lambda = aorg[ vars["lambda"] ];
    double &d7out  = aorg[ vars["d7out"]  ];
    double &sigma  = aorg[ vars["sigma"]  ];
    double &t0     = aorg[ vars["t0"]     ];
    double &scale  = aorg[ vars["scale"]  ];
    const double dmin = find_min_of(delta);

    k7     = 0.003;
    d7out  = 14.98;
    lambda = (1000.0+d7out)/(1000.0+dmin);
    sigma  = 8;

    Fit::LS<double> lsf;
    Fit::Type<double>::Function F(     &dfn, & DeltaFit::Compute   );
    Fit::Type<double>::Function fIni(  &dfn, & DeltaFit::ComputeIni);
    Fit::Type<double>::Function fEnd(  &dfn, & DeltaFit::ComputeEnd);

    const size_t NP = 500;


    if(true)
    {
        ios::wcstream fp("ramp.dat");
        for(size_t i=0;i<=NP;++i)
        {
            const double tt = t[1] + (i*(t[N]-t[1]))/NP;
            fp("%.15g %.15g %.15g\n", tt, dfn.CouplingRamp(tt,aorg,vars), log(tt) );
        }
    }

    //__________________________________________________________________________
    //
    // fit the end to get time scale
    //__________________________________________________________________________
    std::cerr << "Fit relaxing scale..." << std::endl;
    tao::ld(used,false);
    used[ vars["k7"]     ] = true;

    if(!lsf.run(sampleEnd,fEnd,aorg,used,aerr))
    {
        throw exception("cannnot fit end/k7");
    }
    sampleEnd.display(std::cerr,aorg,aerr);
    {
        ios::wcstream fp("delta_end.dat");
        save_data(fp,tEnd,deltaEnd,deltaEndFit);
    }

    //__________________________________________________________________________
    //
    // using the initial scaling, find initial sigma
    //__________________________________________________________________________
    std::cerr << "Fit initial catalytic..." << std::endl;
    tao::ld(used,false);
    used[ vars["sigma"]  ] = true;
    used[ vars["k7"]  ] = true;

    if(!lsf.run(sampleIni,fIni,aorg,used,aerr))
    {
        throw exception("cannnot fit ini/sigma");
    }
    sampleEnd.display(std::cerr,aorg,aerr);
    {
        ios::wcstream fp("delta_ini.dat");
        save_data(fp,tIni,deltaIni,deltaIniFit);
    }

    return 0;

    //__________________________________________________________________________
    //
    // OK, try both...
    //__________________________________________________________________________
    std::cerr << "Fit common k7" << std::endl;
    tao::ld(used,false);
    used[ vars["k7"]     ] = true;
    used[ vars["lambda"] ] = true;
    used[ vars["sigma"]  ] = true;
    //used[ vars["d7out"]  ] = true;

    if(!lsf.run(multiple,F,aorg,used,aerr))
    {
        throw exception("cannnot fit common");
    }

    multiple.display(std::cerr,aorg,aerr);
    {
        ios::wcstream fp("delta_com.dat");
        save_data(fp,tIni,deltaIni,deltaIniFit);
        save_data(fp,tEnd,deltaEnd,deltaEndFit);
    }
    {
        ios::wcstream fp("delta_dfn.dat");
        for(size_t i=0;i<=NP;++i)
        {
            const double tmp  = tIni[1] + (i*(tIni[N_Ini]-tIni[1]))/NP;
            const double tt = log(tmp);
            fp("%.15g %.15g\n", tt, F(tmp,aorg,vars) );
        }
        fp << "\n";
        for(size_t i=0;i<=NP;++i)
        {
            const double tmp  = tEnd[1] + (i*(tEnd[N_End]-tEnd[1]))/NP;
            const double tt   = log(tmp);
            fp("%.15g %.15g\n", tt, F(tmp,aorg,vars) );
        }
    }

    {
        ios::wcstream fp("psi.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g\n", log(t[i]), dfn.ComputePsi(t[i],delta[i],aorg,vars) );
        }
    }


    //__________________________________________________________________________
    //
    // OK, try full
    //__________________________________________________________________________
    return 0;
    std::cerr << "Fit full" << std::endl;
    Fit::Type<double>::Function Full(  &dfn, & DeltaFit::ComputeFull);


    t0    = 170;
    scale = 190;

    tao::ld(used,false);
    used[ vars["k7"]     ] = true;
    used[ vars["sigma"]     ] = true;

    used[ vars["t0"]    ] = true;
    used[ vars["scale"] ] = true;

    if(!lsf.run(sample,Full,aorg,used,aerr))
    {
        throw exception("cannnot fit full");
    }
    sample.display(std::cerr,aorg,aerr);

    {
        ios::wcstream fp("delta_full.dat");
        for(size_t i=0;i<=NP;++i)
        {
            const double tmp = t[1] + (i*(t[N]-t[1]))/NP;
            const double tt  = log(tmp);
            fp("%.15g %.15g\n", tt, Full(tmp,aorg,vars) );
        }
    }

}
YOCTO_PROGRAM_END()


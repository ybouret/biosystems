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

typedef Fit::Samples<double>        Samples;
typedef Fit::Sample<double>         Sample;
typedef Fit::Variables              Variables;
typedef vector<double>              Vector;
typedef Fit::Type<double>::Function FitFunction;
typedef numeric<double>::function   Function;

class DeltaFit
{
public:
    size_t               i0;
    size_t               i1;
    double               t0;
    double               t1;
    const array<double> *pAorg;
    const Variables     *pVars;
    Function             f;
    Function             g;
    derivatives<double>  drvs;
    
    explicit DeltaFit() :
    i0(0),
    i1(0),
    t0(0),
    t1(0),
    pAorg(0),
    pVars(0),
    f( this, & DeltaFit::FunctionIni ),
    g( this, & DeltaFit::FunctionEnd )
    {}

    virtual ~DeltaFit() throw() {}

    double t_reponse;


    double Compute(const double t, const array<double> &aorg, const Variables &vars )
    {
        if(t<=t0)
        {
            return ComputeIni(t,aorg,vars);
        }
        else if(t>=t1)
        {
            return ComputeEnd(t,aorg,vars);
        }
        else
        {
            pAorg = &aorg;
            pVars = &vars;
            const double dd  = t1-t0; assert(dd>0);
            const double dd2 = dd * dd;
            const double dd3 = dd2 * dd;
            const double dd4 = dd2*dd2;
            const double f0  = f(t0);
            const double df0 = drvs.compute(f,t0,1e-4);
            const double g1  = g(t1);
            const double dg1 = drvs.compute(g,t1,1e-4);
            const double A   = g1-f0 - dd * df0;
            const double B   = dg1 - df0;
            const double __beta  = ((3*dd2)*A - dd3*B)/dd4;
            const double __gamma = ((-2*dd)*A + dd2*B)/dd4;
            const double dt = t-t0;
            return f0 + df0 * dt + __beta * dt*dt + __gamma * dt *dt * dt;
        }
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

    double FunctionIni(const double t)
    {
        assert(pAorg);
        assert(pVars);
        return ComputeIni(t,*pAorg,*pVars);
    }


    double ComputeEnd(const double t, const array<double> &aorg, const Variables &vars )
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];
        const double tau7    = k7*t;
        const double tau6    = lambda*tau7;

        return 1000.0 * ( (1.0+d7out/1000.0) * (Growth(tau7)/Growth(tau6)) - 1.0 );
    }

    double FunctionEnd(const double t)
    {
        assert(pAorg);
        assert(pVars);
        return ComputeEnd(t,*pAorg,*pVars);
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

    const double t_short = strconv::to<double>(argv[1],"t_short");
    const double t_long  = strconv::to<double>(argv[2],"t_long");
    if(t_long<t_short)
    {
        throw exception("t_long=%g<t_short=%g", t_long, t_short);
    }


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
    size_t &i0 = dfn.i0;
    size_t &i1 = dfn.i1;
    i0=1;
    i1=N;

    for(size_t i=1;i<=N;++i)
    {
        if(t[i]<=t_short)
        {
            tIni.push_back(t[i]);
            deltaIni.push_back(delta[i]);
            i0=max_of(i0,i);
        }
        if(t[i]>=t_long)
        {
            tEnd.push_back(t[i]);
            deltaEnd.push_back(delta[i]);
            i1 = min_of(i1,i);
        }
    }
    if(i0>=i1)
    {
        throw exception("overlap of short time and long time points");
    }

    dfn.t0 = t[i0];
    dfn.t1 = t[i1];


    const size_t nIni = tIni.size(); if(nIni<=0) throw exception("no ini point...");
    const size_t nEnd = tEnd.size(); if(nEnd<=0) throw exception("no end point...");



    {
        ios::wcstream fp("delta0.dat");
        save_data(fp,t,delta,delta);
    }
    //__________________________________________________________________________
    //
    // create sub samples
    //__________________________________________________________________________
    Fit::Samples<double> multiple(3);

    Vector       deltaIniFit(nIni);
    Fit::Sample<double> &sampleIni = multiple.add(tIni,deltaIni,deltaIniFit);

    Vector       deltaEndFit(nEnd);
    Fit::Sample<double> &sampleEnd = multiple.add(tEnd,deltaEnd,deltaEndFit);

    //__________________________________________________________________________
    //
    // global variables
    //__________________________________________________________________________
    Fit::Variables &vars = multiple.variables;
    vars << "k7" << "lambda" << "d7out" << "sigma";
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
    const double dmin = find_min_of(delta);

    k7     = 0.003;
    d7out  = 14.98;
    lambda = (1000.0+d7out)/(1000.0+dmin);
    sigma  = 8;

    Fit::LS<double> lsf;
    Fit::Type<double>::Function f(  &dfn, & DeltaFit::ComputeIni);
    Fit::Type<double>::Function g(  &dfn, & DeltaFit::ComputeEnd);
    Fit::Type<double>::Function F(  &dfn, & DeltaFit::Compute);

    const size_t NP = 1000;


    //__________________________________________________________________________
    //
    // fit the end to get time scale
    //__________________________________________________________________________
    std::cerr << "Fit relaxing scale..." << std::endl;
    tao::ld(used,false);
    used[ vars["k7"]     ] = true;

    if(!lsf.run(sampleEnd,g,aorg,used,aerr))
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

    if(!lsf.run(sampleIni,f,aorg,used,aerr))
    {
        throw exception("cannnot fit ini/sigma");
    }
    sampleEnd.display(std::cerr,aorg,aerr);
    {
        ios::wcstream fp("delta_ini.dat");
        save_data(fp,tIni,deltaIni,deltaIniFit);
    }

    //__________________________________________________________________________
    //
    // Fit both
    //__________________________________________________________________________
    std::cerr << "Fit Both" << std::endl;
    tao::ld(used,false);
    //used[ vars["lambda"]  ] = true;
    used[ vars["k7"]      ] = true;
    used[ vars["lambda"]      ] = true;
    //used[ vars["d7out"]      ] = true;
    used[ vars["sigma"]  ] = true;

    if(!lsf.run(sample,F,aorg,used,aerr))
    {
        throw exception("cannnot fit both");
    }

    sample.display(std::cerr,aorg,aerr);

    (void)sampleIni.computeD2(f,aorg);
    (void)sampleEnd.computeD2(g,aorg);
    {
        ios::wcstream fp("delta_ini.dat");
        save_data(fp,tIni,deltaIni,deltaIniFit);
    }

    {
        ios::wcstream fp("delta_end.dat");
        save_data(fp,tEnd,deltaEnd,deltaEndFit);
    }

    {
        ios::wcstream fp("delta_dfn.dat");

        for(size_t i=0;i<=NP;++i)
        {
            const double tt  = t[1] + ( i * (t[N]-t[1]) ) / NP;
            const double f   = dfn.Compute(tt,aorg,vars);
            fp("%g %g\n", log(tt), f);
        }

    }

}
YOCTO_PROGRAM_END()


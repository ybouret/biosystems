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
    const array<double> *pAorg;
    const Variables     *pVars;
    derivatives<double>  drvs;
    
    explicit DeltaFit() :
    pAorg(0),
    pVars(0)
    {}

    virtual ~DeltaFit() throw() {}


    double ComputeIni(const double t, const array<double> &aorg, const Variables &vars )
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];
        const double sigma  = aorg[ vars["sigma"]  ];
        const double tau7   = k7*t;
        const double tau6   = lambda*tau7;

        return 1000.0 * ( (1.0+d7out/1000.0) * (Xi(tau7,sigma)/Xi(tau6,sigma/lambda)) - 1.0 );
    }

    double Compute(const double t, const array<double> &aorg, const Variables &vars )
    {

        const double t_c    = aorg[ vars["t_c"]    ];
        if(t<=t_c)
        {
            return ComputeIni(t,aorg,vars);
        }
        else
        {
            const double k7     = aorg[ vars["k7"]     ];
            const double lambda = aorg[ vars["lambda"] ];
            const double d7out  = aorg[ vars["d7out"]  ];
            const double sigma  = aorg[ vars["sigma"]  ];
            const double Beta   = aorg[ vars["Beta"] ];
            const double tp     = t-t_c;
            const double tau7p  = k7*tp;
            const double tau6p  = lambda * tau7p;
            const double tau7c  = k7*t_c;
            const double tau6c  = lambda*tau7c;
            const double beta7c = Xi(tau7c,sigma);
            const double beta6c = Xi(tau6c,sigma/lambda);
            const double e7     = exp( -tau7p );
            const double e6     = exp( -tau6p );
            const double beta7  = beta7c * e7 + Beta * (1.0-e7);
            const double beta6  = beta6c * e6 + Beta * (1.0-e6);
            return 1000.0 * ( (1.0+d7out/1000.0) * beta7/beta6 - 1.0);
        }
    }

    double __beta7c(const array<double> &aorg, const Variables &vars ) const
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double sigma  = aorg[ vars["sigma"]  ];
        const double t_c    = aorg[ vars["t_c"]    ];
        const double tau7c  = k7*t_c;
        return Xi(tau7c,sigma);
    }

    double __beta6c(const array<double> &aorg, const Variables &vars ) const
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double sigma  = aorg[ vars["sigma"]  ];
        const double t_c    = aorg[ vars["t_c"]    ];
        const double lambda = aorg[ vars["lambda"] ];
        const double tau7c  = k7*t_c;
        return Xi(tau7c*lambda,sigma/lambda);
    }



#if 0
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
#endif

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

    

    {
        ios::wcstream fp("delta0.dat");
        save_data(fp,t,delta,delta);
    }

    const double dmin = find_min_of(delta);
    const double dmax = find_max_of(delta);
    const double dlim = dmin + 0.25 * (dmax-dmin);
    Vector t_ini(N,as_capacity);
    Vector delta_ini(N,as_capacity);
    for(size_t i=1;i<=N;++i)
    {
        if(delta[i]<=dlim)
        {
            t_ini.push_back(t[i]);
            delta_ini.push_back(delta[i]);
        }
        else
            break;
    }
    const size_t N_ini = t_ini.size();
    Vector deltaFit_ini(N_ini);
    Sample sample_ini(t_ini,delta_ini,deltaFit_ini);




    //__________________________________________________________________________
    //
    // global variables
    //__________________________________________________________________________
    Fit::Variables &vars = sample.variables;
    vars << "k7" << "lambda" << "t_c" << "d7out" << "sigma" << "Beta";
    sample_ini.variables = vars;


    const size_t nvar = vars.size();
    Vector       aorg(nvar);
    Vector       aerr(nvar);
    vector<bool> used(nvar,false);
    std::cerr << "nvar=" << nvar << std::endl;
    double &k7     = aorg[ vars["k7"]     ];
    double &lambda = aorg[ vars["lambda"] ];
    double &d7out  = aorg[ vars["d7out"]  ];
    double &sigma  = aorg[ vars["sigma"]  ];
    double &t_c    = aorg[ vars["t_c"] ];
    double &Beta   = aorg[ vars["Beta"] ];

    k7     = 0.003;
    d7out  = 14.98;
    lambda = (1000.0+d7out)/(1000.0+dmin);
    sigma  = 100;
    t_c    = 90;

    Fit::LS<double>             lsf;
    Fit::Type<double>::Function F(  &dfn, & DeltaFit::Compute);
    Fit::Type<double>::Function f(  &dfn, & DeltaFit::ComputeIni);

    used[ vars["k7"]    ] = true;
    used[ vars["sigma"]    ] = true;
    used[ vars["lambda"]    ] = true;

    if(!lsf.run(sample_ini,f,aorg,used,aerr))
    {
        throw exception("couldn't fit init");
    }

    sample_ini.display(std::cerr, aorg, aerr);
    {
        ios::wcstream fp("delta1.dat");
        save_data(fp,t_ini,delta_ini,deltaFit_ini);
    }

    tao::ld(used,false);
    Beta = 0.5 * (dfn.__beta6c(aorg,vars)+dfn.__beta7c(aorg,vars));
    used[ vars["Beta"]    ] = true;
    std::cerr << "Starting from Beta=" << Beta << std::endl;
    if(!lsf.run(sample,F,aorg,used,aerr))
    {
        throw exception("couldn't fit Beta");
    }

    tao::ld(used,false);

    used[ vars["Beta"]   ] = true;
    used[ vars["k7"]     ] = true;
    used[ vars["lambda"] ] = true;
    //used[ vars["d7out"] ] = true;

    if(!lsf.run(sample,F,aorg,used,aerr))
    {
        throw exception("couldn't fit Beta/k7/lambda");
    }

    tao::ld(used,false);
    used[ vars["sigma"]     ] = true;
    if(!lsf.run(sample,F,aorg,used,aerr)) throw ("no sigma level2");


    sample.display(std::cerr, aorg, aerr);
    {
        ios::wcstream fp("delta2.dat");
        save_data(fp,t,delta,deltaFit);
    }

    const size_t NP = 1000;
    {
        ios::wcstream fp("dfn.dat");
        for(size_t i=0;i<=NP;++i)
        {
            const double tt = t[1] + 2.0 * (i*(t[N]-t[1])/double(NP));
            fp("%g %g\n", log(tt), F(tt,aorg,vars));
        }
    }



}
YOCTO_PROGRAM_END()


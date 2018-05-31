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
    explicit DeltaFit() : t_reponse(60.0) {}
    virtual ~DeltaFit() throw() {}

    double t_reponse;

    double Compute( const double t, const array<double> &aorg, const Variables &vars )
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

    double Compute0( const double t, const array<double> &aorg, const Variables &vars )
    {
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];

        const double tau7    = k7*t;
        const double tau6    = lambda*tau7;
        const double beta7   = Growth(tau7);
        const double beta6   = Growth(tau6);

        return 1000.0 * ( (1.0+d7out/1000.0) * (beta7/beta6) - 1.0 );
    }

    double Ratio(const double t, const array<double> &aorg, const Variables &vars )
    {
        return Compute(t,aorg,vars)/Compute0(t,aorg,vars);
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

    if(false)
    {
        testXi();
    }

    if(argc<=1)
    {
        throw exception("usage: %s t_cut",program);
    }
    const double t_cut = strconv::to<double>(argv[1],"t_cut");


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
    // cut long time data
    //__________________________________________________________________________
    Vector tEnd(N,as_capacity);
    Vector deltaEnd(N,as_capacity);
    for(size_t i=1;i<=N;++i)
    {
        if(t[i]>=t_cut)
        {
            tEnd.push_back(t[i]);
            deltaEnd.push_back(delta[i]);
        }
    }
    const size_t Ncut = tEnd.size();
    Vector       deltaEndFit(Ncut);
    Sample       sampleEnd(tEnd,deltaEnd,deltaEndFit);

    //__________________________________________________________________________
    //
    // global variables
    //__________________________________________________________________________
    Variables &vars = sample.variables;
    vars << "k7" << "lambda" << "sigma" << "psi" << "d7out";
    sampleEnd.variables = vars;

    sampleEnd.link();
    sample.link();

    const size_t nvar = vars.size();
    Vector       aorg(nvar);
    Vector       aerr(nvar);
    vector<bool> used(nvar,false);

    double &k7     = aorg[ vars["k7"]     ];
    double &lambda = aorg[ vars["lambda"] ];
    double &sigma  = aorg[ vars["sigma"]  ];
    double &psi    = aorg[ vars["psi"]    ];
    double &d7out  = aorg[ vars["d7out"]  ];

    // initialize variables
    const double dmin = find_min_of(delta);

    k7     = 0.003;
    lambda = 1.0;
    psi    = 0.0;
    sigma  = 0.01;
    d7out  = 15.00;

    lambda = (1000.0+d7out)/(1000.0+dmin);

    used[ vars["k7"]     ] = true;
    used[ vars["lambda"] ] = true;
    //used[ vars["d7out"] ]  = false;

    Fit::LS<double> lsf;
    DeltaFit        dfn;
    Fit::Type<double>::Function F(  &dfn, & DeltaFit::Compute  );
    Fit::Type<double>::Function F0( &dfn, & DeltaFit::Compute0 );

    if( !lsf.run(sampleEnd,F,aorg, used, aerr) )
    {
        throw exception("couldn't fit k7/lam");
    }
    else
    {
        sampleEnd.display(std::cerr,aorg,aerr);

        {
            ios::wcstream fp("delta_fit_end.dat");
            save_data(fp,tEnd,deltaEnd,deltaEndFit);
        }

        {
            ios::wcstream fp("delta_fit0.dat");
            (void)sample.computeD2(F,aorg);
            save_data(fp,t,delta,deltaFit);
        }


        {
            ios::wcstream fp("delta_dfn0.dat");
            const size_t M = 1000;
            for(size_t i=0;i<=M;++i)
            {
                const double tt = t[1] + (i*(t[N]-t[1]))/M;
                fp("%.15g %.15g\n", log(tt), F(tt,aorg,vars));
            }
        }
    }


#if 0
    if( !lsf.run(sample,F,aorg, used, aerr) )
    {
        throw exception("couldn't fit k7/lam");
    }
    else
    {
        sample.display(std::cerr,aorg,aerr);
        for(size_t i=N;i>0;--i)
        {
            ratio[i] = delta[i]/F0(t[i],aorg,vars);
        }
        {
            ios::wcstream fp("delta_fit0.dat");
            save_data(fp,t,delta,deltaFit,&ratio);
        }




        std::cerr << "d7ini=" << dfn.d7ini(aorg,vars) << std::endl;
        {
            ios::wcstream fp("delta_dfn0.dat");
            const size_t M = 1000;
            for(size_t i=0;i<=M;++i)
            {
                const double tt = t[1] + (i*(t[N]-t[1]))/M;
                fp("%.15g %.15g\n", log(tt), F(tt,aorg,vars));
            }
        }
    }
#endif

}
YOCTO_PROGRAM_END()


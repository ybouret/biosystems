#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fit/fit.hpp"

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

static inline void BuildHull( Vector &X, Vector &Y, const array<double> &x, const array<double> &y )
{
    assert(x.size()==y.size());
    const size_t N = x.size();
    X.free();
    Y.free();
    X.ensure(N);
    Y.ensure(N);

    // take minimal point
    double xcurr = x[1];
    double ycurr = y[1];
    size_t icurr = 1;
    for(size_t i=N;i>1;--i)
    {
        const double ytmp = y[i];
        if(ytmp<ycurr)
        {
            xcurr = x[i];
            ycurr = ytmp;
            icurr = i;
        }
    }

    X.push_back(xcurr);
    Y.push_back(ycurr);

    for(size_t i=icurr;i<N;)
    {
        std::cerr << "@" << i << ": xcurr=" << xcurr << ", ycurr=" << ycurr << std::endl;
        size_t jhull = i+1;
        double xhull = x[jhull];
        double yhull = y[jhull];
        double angle = math::Atan2(yhull-ycurr,xhull-xcurr);
        std::cerr << "\t@" << jhull <<": xhull=" << xhull << ", yhull=" << yhull << ", angle=" << angle << " --" << std::endl;
        assert(xhull>xcurr);
        for(size_t j=jhull+1;j<=N;++j)
        {
            const double xj = x[j];
            const double yj = y[j];
            const double aj = math::Atan2(yj-ycurr,xj-xcurr);
            std::cerr << "\t@" << j <<": xhull=" << xj << ", yhull=" << yj << ", angle=" << aj << std::endl;
            if(aj<=angle)
            {
                angle=aj;
                xhull=xj;
                yhull=yj;
            }
        }
        X.push_back(xhull);
        Y.push_back(yhull);
        i=jhull;
        xcurr=xhull;
        ycurr=yhull;
    }


}

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

    testXi();


    //__________________________________________________________________________
    //
    // Loading Data
    //__________________________________________________________________________

    const string workdir = "src/lithium/doc/";

    vector<double> t;      //! col 1
    vector<double> delta;  //! col 3
    vector<double> deltaFit;
    vector<double> ratio;

    {
        const string filename = workdir + "nhe1_delta7_full_15mM_37.txt";
        std::cerr << "Loading " << filename << std::endl;
        data_set<double> ds;
        ds.use(1,t);
        ds.use(3,delta);
        ios::icstream fp(filename);
        ds.load(fp);
        std::cerr << "...done" << std::endl;
    }
    const size_t N = t.size();
    deltaFit.make( N );
    ratio.make(N);
    Sample sample(t,delta,deltaFit);

#if 0
    vector<double> t_hull;
    vector<double> d_hull;
    BuildHull(t_hull,d_hull,t,delta);
    {
        ios::wcstream fp("hull.dat");
        save_data(fp,t_hull,d_hull,d_hull);
    }
#endif
    
    //__________________________________________________________________________
    //
    // global variables
    //__________________________________________________________________________
    Variables &vars = sample.variables;
    vars << "k7" << "lambda" << "sigma" << "psi" << "d7out";


    const size_t nvar = vars.size();
    Vector       aorg(nvar);
    Vector       aerr(nvar);
    vector<bool> used(nvar,false);

    double &k7     = aorg[ vars["k7"]     ];
    double &lambda = aorg[ vars["lambda"] ];
    double &sigma  = aorg[ vars["sigma"]  ];
    double &psi    = aorg[ vars["psi"]    ];
    double &d7out  = aorg[ vars["d7out"]  ];

    k7     = 0.003;
    lambda = 1.01;
    psi    = 0.0;
    sigma  = 0.01;
    d7out  = 15.00;

    used[ vars["k7"]     ] = true;
    used[ vars["lambda"] ] = true;
    used[ vars["d7out"] ]  = true;

    Fit::LS<double> lsf;
    DeltaFit        dfn;
    Fit::Type<double>::Function F(  &dfn, & DeltaFit::Compute  );
    Fit::Type<double>::Function F0( &dfn, & DeltaFit::Compute0 );



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
}
YOCTO_PROGRAM_END()


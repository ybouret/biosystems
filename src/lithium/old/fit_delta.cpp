#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/lua/lua-funcs.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/math/fit/fit.hpp"

using namespace yocto;
using namespace math;

typedef vector<double> Vector;
typedef array<double>  Array;

typedef ode::driverCK<double>::type  Solver;
typedef ode::Field<double>::Equation DiffEq;
typedef Fit::Type<double>::Function  FitFunction;
typedef Fit::Variables               Variables;

#define NVAR 3
static const size_t IA = 1;
static const size_t I6 = 2;
static const size_t I7 = 3;

class System
{
public:

    Lua::Function<double> h;
    double                gamma;
    double                eta;
    double                mu;
    double                lambda;
    double                Theta;
    double                k7;
    DiffEq                diffeq;
    const double          h0;

    virtual ~System() throw()
    {
    }

#define __LUA_GET(NAME) NAME( vm->Get<double>(#NAME) )
    explicit System( Lua::State::Pointer &vm ) :
    h(vm,"h"),
    __LUA_GET(gamma),
    __LUA_GET(eta),
    __LUA_GET(mu),
    __LUA_GET(lambda),
    __LUA_GET(Theta),
    __LUA_GET(k7),
    diffeq(this,&System::compute),
    h0( h(0) )
    {
    }

    void compute( Array &dYdt, double t, const Array &Y )
    {
        const double alpha   = Y[IA];
        const double hh      = h(t)/h0;
        {
            const double scaling = gamma * hh;
            dYdt[IA] = k7 * (scaling - (scaling+eta) * alpha );
        }
        const double Beta = Theta + mu*hh*(1.0-alpha);

        dYdt[I7] =          k7 * ( Beta - Y[I7] );
        dYdt[I6] = lambda * k7 * ( Beta - Y[I6] );

    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(System);
};



class Generator
{
public:
    Solver odeint;

    virtual ~Generator() throw() {}

    explicit Generator() :
    odeint(1e-4)
    {
        odeint.start(NVAR);
    }

    void compute( Array &Y, const double t, DiffEq &diffeq )
    {
        static const double dt_max = 0.5;
        double h       = dt_max;
        tao::ld(Y,0);
        //odeint( diffeq, Y, 0, t, h, NULL);
        //return;
        double t_start = 0;
        for(;;)
        {
            double t_final = t_start + dt_max;
            bool   done    = false;
            if(t_final>=t)
            {
                t_final = t;
                done    = true;
            }
            h = min_of(h,dt_max);
            odeint( diffeq, Y, t_start, t_final, h, NULL);
            t_start = t_final;
            if(done) break;
        }
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Generator);
};

class DeltaLi
{
public:
    System    sys;
    Generator gen;
    Vector    Y;

    explicit DeltaLi( Lua::State::Pointer &vm ) :
    sys(vm),
    gen(),
    Y(NVAR)
    {

    }

    virtual ~DeltaLi() throw() {}

    double compute( const double t, const Array &aorg, const Variables &vars )
    {
        const double d7out = vars(aorg,"d7out");
        gen.compute(Y,t,sys.diffeq);

        return 1000.0 * ( (1.0+1e-3*d7out) * Y[I7]/Y[I6] - 1.0 );
    }



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(DeltaLi);
};


YOCTO_PROGRAM_START()
{
    Lua::State::Pointer VM( new Lua::State() );
    if(argc>1)
    {
        VM->DoFile(argv[1]);
    }

    DeltaLi   dfn(VM);
    Variables vars;
    vars << "d7out";
    const size_t nfit = vars.size();
    Vector aorg(nfit);

    double &d7out = vars(aorg,"d7out");

    d7out = 14.98;

    FitFunction F( &dfn, & DeltaLi::compute);

    ios::ocstream::overwrite("intg.dat");
    for(double lt=1.5;lt<=8.2;lt+=0.1)
    {
        const double t = exp(lt);
        std::cerr << "t=" << t << std::endl;
        const double d7 = F(t,aorg,vars);
        {
            ios::acstream fp("intg.dat");
            fp("%.15g", lt);
            fp(" %.15g",dfn.Y[IA]);
            fp(" %.15g",d7);
            fp(" %.15g",dfn.Y[I6]);
            fp(" %.15g",dfn.Y[I7]);
            fp(" %.15g",d7);

            fp("\n");
        }
    }
}
YOCTO_PROGRAM_END()




#if 0
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
    Function             fIni;
    Function             fEnd;
    double               tIni;
    double               tEnd;
    explicit DeltaFit() :
    pAorg(0),
    pVars(0),
    fIni( this, & DeltaFit:: FunctionIni ),
    fEnd( this, & DeltaFit:: FunctionEnd ),
    tIni(0),
    tEnd(0)
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

    double ComputeEnd(const double t, const array<double> &aorg, const Variables &vars )
    {
        const double t_c    = aorg[ vars["t_c"] ];
        const double k7     = aorg[ vars["k7"]     ];
        const double lambda = aorg[ vars["lambda"] ];
        const double d7out  = aorg[ vars["d7out"]  ];
        const double gamma  = aorg[ vars["gamma"]  ];
        const double tp     = t-t_c;
        const double tau7p  = k7*tp;
        const double tau6p  = lambda * tau7p;
        const double tau7c  = k7*t_c;
        const double tau6c  = lambda*tau7c;
        const double beta7c = 1.0-exp(-tau7c);
        const double beta6c = 1.0-exp(-tau6c);
        const double e7     = exp( -tau7p );
        const double e6     = exp( -tau6p );
        const double beta7  = beta7c * e7 + gamma * (1.0-e7);
        const double beta6  = beta6c * e6 + gamma * (1.0-e6);
        return 1000.0 * ( (1.0+d7out/1000.0) * beta7/beta6 - 1.0);
    }

    double Compute(const double t, const array<double> &aorg, const Variables &vars )
    {

        if(t<=tIni)
        {
            return ComputeIni(t,aorg,vars);
        }
        else if(t>=tEnd)
        {
            return ComputeEnd(t,aorg,vars);
        }
        else
        {
            pAorg = &aorg;
            pVars = &vars;
            const double f1   = fIni(tIni);
            const double f2   = fEnd(tEnd);
            const double df1  = drvs.compute(fIni,tIni,1e-4);
            const double df2  = drvs.compute(fEnd,tEnd,1e-4);
            polynomial<double> P;
            __poly::gap(P, tIni, f1, df1, tEnd, f2, df2);
            return P(t-tIni);
        }
    }


    double FunctionIni(const double t)
    {
        assert(pAorg);
        assert(pVars);
        return ComputeIni(t,*pAorg,*pVars);
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

    if( argc <= 2 )
    {
        throw exception("usage: %s tIni tEnd", program);
    }
    {
        size_t iarg=0;
        dfn.tIni = strconv::to<double>(argv[++iarg],"tIni");
        //dfn.tCut = strconv::to<double>(argv[++iarg],"tCut");
        dfn.tEnd = strconv::to<double>(argv[++iarg],"tEnd");
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

    Samples        multi;
    

    {
        ios::wcstream fp("delta0.dat");
        save_data(fp,t,delta,delta);
    }

    const double dmin = find_min_of(delta);
    Vector t_ini(N,as_capacity);
    Vector delta_ini(N,as_capacity);
    Vector t_end(N,as_capacity);
    Vector delta_end(N,as_capacity);
    for(size_t i=1;i<=N;++i)
    {
        if(t[i]<=dfn.tIni)
        {
            t_ini.push_back(t[i]);
            delta_ini.push_back(delta[i]);
        }

        if(t[i]>=dfn.tEnd)
        {
            t_end.push_back(t[i]);
            delta_end.push_back(delta[i]);
        }
    }
    const size_t N_ini = t_ini.size();
    Vector  deltaFit_ini(N_ini);
    Sample &sample_ini = multi.add(t_ini,delta_ini,deltaFit_ini);

    const size_t N_end = t_end.size();
    Vector  deltaFit_end(N_end);
    Sample &sample_end = multi.add(t_end,delta_end,deltaFit_end);



    //__________________________________________________________________________
    //
    // global variables
    //__________________________________________________________________________
    Fit::Variables &vars = sample.variables;
    vars << "k7" << "lambda" << "d7out" << "sigma" << "gamma" << "t_c";
    multi.variables       = vars;
    sample_ini.variables  = vars;
    sample_end.variables  = vars;

    const size_t nvar = vars.size();
    Vector       aorg(nvar);
    Vector       aerr(nvar);
    vector<bool> used(nvar,false);
    std::cerr << "nvar=" << nvar << std::endl;
    double &k7     = aorg[ vars["k7"]     ];
    double &lambda = aorg[ vars["lambda"] ];
    double &d7out  = aorg[ vars["d7out"]  ];
    double &sigma  = aorg[ vars["sigma"]  ];
    double &t_c    = aorg[ vars["t_c"]    ];
    double &gamma  = aorg[ vars["gamma"]  ];


    k7     = 0.003;
    d7out  = 14.98;
    lambda = (1000.0+d7out)/(1000.0+dmin);
    sigma  = 100;
    t_c    = 120;
    gamma  = 0.1;

    Fit::LS<double>             lsf;
    Fit::Type<double>::Function F(  &dfn, & DeltaFit::Compute);
    Fit::Type<double>::Function F_ini(  &dfn, & DeltaFit::ComputeIni);
    Fit::Type<double>::Function F_end(  &dfn, & DeltaFit::ComputeEnd);

    tao::ld(used,false);
    used[ vars["k7"]     ] = true;
    used[ vars["sigma"]  ] = true;
    used[ vars["lambda"] ] = true;

    if(!lsf.run(sample_ini,F_ini,aorg,used,aerr))
    {
        throw exception("couldn't fit init");
    }

    sample_ini.display(std::cerr, aorg, aerr) << std::endl;
    {
        ios::wcstream fp("delta1.dat");
        save_data(fp,t_ini,delta_ini,deltaFit_ini);
    }


    const size_t NP = 1000;
    {
        ios::wcstream fp("dfn1.dat");
        for(size_t i=0;i<=NP;++i)
        {
            const double tt = t_ini.front()/2 +  (i*(dfn.tIni-t_ini.front()/2)/double(NP));
            fp("%g %g\n", log(tt), F_ini(tt,aorg,vars));
        }
    }

    tao::ld(used,false);
    used[ vars["gamma"] ] = true;
    if(!lsf.run(sample_end,F_end,aorg,used,aerr))
    {
        throw exception("couldn't fit end");
    }
    sample_end.display(std::cerr, aorg, aerr) << std::endl;
    {
        ios::wcstream fp("delta2.dat");
        save_data(fp,t_end,delta_end,deltaFit_end);
    }

    {
        ios::wcstream fp("dfn2.dat");
        for(size_t i=0;i<=NP;++i)
        {
            const double tt = dfn.tEnd +  (2*i*(t_end.back()-dfn.tEnd)/double(NP));
            fp("%g %g\n", log(tt), F_end(tt,aorg,vars));
        }
    }

    tao::ld(used,false);
    vars(used,"k7")     = true;
    //vars(used,"lambda") = true;
    //vars(used,"sigma")  = true;
    //vars(used,"d7out")  = true;

    if(!lsf.run(sample,F,aorg,used,aerr)) throw exception("couldn't fit both");
    sample.display(std::cerr, aorg, aerr) << std::endl;

    {
        ios::wcstream fp("dfn0.dat");
        const double tt0 = t[1]/2;
        const double tt1 = t[1] + 1.5*(t[N]-t[1]);
        const double dtt = tt1-tt0;
        for(size_t i=0;i<=NP;++i)
        {
            const double tt = tt0 + (i*dtt)/NP;
            fp("%g %g\n", log(tt), F(tt,aorg,vars));
        }
    }



}
YOCTO_PROGRAM_END()
#endif


#include "yocto/program.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/fit/glsf.hpp"
#include "yocto/container/utils.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/math/dat/linear.hpp"

using namespace yocto;
using namespace math;

static double dLi7_out = 0;


class Lithium
{
public:

    static double BumpFull(const double tau, const double sigma)
    {
        return ( exp(-tau) - exp( -sigma * tau ) )/(sigma-1.0);
    }

    static double BumpZero(const double tau, const double sigma)
    {
        const double st = (sigma-1.0) * tau;
        return tau * exp( -tau ) * ( 1.0 - st*(1.0-0.5*st) );
    }

    static double Bump(const double tau, const double sigma)
    {
        if( Fabs(sigma-1.0)<1e-6 )
        {
            return BumpZero(tau,sigma);
        }
        else
        {
            return BumpFull(tau,sigma);
        }
    }

    static void TestBump()
    {
        const double sigma[] = { 0.9, 0.99, 0.999, 1.001, 1.01, 1.1 };
        const double dtau = 0.001;
        ios::wcstream fp("bump.dat");
        for(size_t i=0;i<sizeof(sigma)/sizeof(sigma[0]);++i)
        {
            const double sig = sigma[i];
            for(double tau=dtau;tau<=5.0;tau+=dtau)
            {
                fp("%g %g %g\n", tau, BumpFull(tau,sig), BumpZero(tau,sig));
            }
            fp("\n");
        }
    }

    static double Grow(const double tau)
    {
        return (1.0-exp(-tau));
    }

    static double Core(const double tau, const double phi6, const double sigma )
    {
        const double B = Bump(tau,sigma);
        return (1.0+phi6)*B/(Grow(tau)+phi6*B);
    }

    Lithium()
    {
    }

    ~Lithium()
    {
    }

    inline double Fit(const double u, const array<double> &a )
    {
        const double du     = a[1];
        const double gamma7 = a[2];
        const double phi6   = a[3];
        const double sigma  = a[4];

        //const double fac6 = phi6 / (1.0+phi6);
        const double tau  = exp(u-du);
        return gamma7 * phi6/(1.0+phi6) * Core(tau,phi6,sigma);
    }

    void saveFit(const char *fn, const array<double> &U, const array<double> &a)
    {
        ios::wcstream fp(fn);
        const double umin = U[1];
        const double umax = U[U.size()];
        const size_t M = 1000;
        for(size_t i=0;i<=M;++i)
        {
            const double u = umin + i * (umax-umin) / double(M);
            fp("%g %g\n", u, Fit(u,a) );
        }
    }

};



YOCTO_PROGRAM_START()
{
    Lithium::TestBump();
    if(argc<=2) throw exception("usage: %s dLi7_out dLi7.dat",program);

    dLi7_out = strconv::to_double(argv[1],"dLi7_out");

    std::cerr << "dLi7_out=" << dLi7_out << std::endl;

    vector<double> t;
    vector<double> dLi7;

    {
        data_set<double> ds;
        ds.use(1,t);
        ds.use(2,dLi7);
        ios::icstream fp(argv[2]);
        ds.load(fp);
    }
    const size_t N = t.size();
    std::cerr << "#data=" << N << std::endl;
    std::cerr << "t   =" << t << std::endl;
    std::cerr << "dLi7=" << dLi7 << std::endl;

    vector<double> u(N);
    vector<double> Omega(N);
    vector<double> OmegaFit(N);

    for(size_t i=1;i<=N;++i)
    {
        u[i]     = Log(t[i]);
        Omega[i] = (1e-3*(dLi7[i]-dLi7_out))/(1.0+1e-3*dLi7_out);
    }

    {
        ios::wcstream fp("omega.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%.15g %.15g %.15g\n", u[i], Omega[i], dLi7[i]);
        }
    }


    Lithium Li;
    GLS<double>::Samples samples(1);
    GLS<double>::Sample &sample = samples.append(u,Omega,OmegaFit);
    GLS<double>::Function F( &Li, & Lithium::Fit );
    vector<double> aorg(4);
    vector<bool>   used( aorg.size(), false );
    vector<double> aerr( aorg.size(), 0 );
    double &du     = aorg[1];
    double &gamma7 = aorg[2];
    double &phi6   = aorg[3];
    double &sigma  = aorg[4];

    samples.prepare( aorg.size() );

    const double Omega0 = Omega[1];
    sigma  = 10;
    phi6   = 1;
    const double fac6 = phi6/(1.0+phi6);
    gamma7 = Omega0/fac6;
    {
        vector<double> uh;
        linear_find(0.5*Omega0,uh,u, Omega);
        std::cerr << "uh=" << uh << std::endl;
        if(uh.size()<=0)
        {
            throw exception("couldn't find half value for Omega");
        }
        du = tao::sum(uh)/uh.size();
    }
    std::cerr << "gamma7=" << gamma7 << std::endl;
    std::cerr << "du    =" << du     << std::endl;

    sample.computeD2(F,aorg);
    {
        ios::wcstream fp("fit0.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%.15g %.15g %.15g\n", u[i], Omega[i], OmegaFit[i]);
        }
    }
    std::cerr << "du standalone" << std::endl;
    used[1] = true;
    if( !samples.fit_with(F,aorg,used,aerr) )
    {
        throw exception("couldn't find du");
    }
    GLS<double>::display(std::cerr,aorg,aerr);
    Li.saveFit("omfit0.dat",u,aorg);



    {
        ios::wcstream fp("fit1.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%.15g %.15g %.15g\n", u[i], Omega[i], OmegaFit[i]);
        }
    }
    Li.saveFit("omfit1.dat",u,aorg);


    std::cerr << "sigma/phi6" << std::endl;
    used[1] = false;
    used[4] = true;
    used[3] = true;
    used[2] = true;
    if( !samples.fit_with(F,aorg,used,aerr) )
    {
        throw exception("couldn't find sigma/phi6");
    }
    GLS<double>::display(std::cerr,aorg,aerr);
    {
        ios::wcstream fp("fit2.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%.15g %.15g %.15g\n", u[i], Omega[i], OmegaFit[i]);
        }
    }
    Li.saveFit("omfit2.dat",u,aorg);


    std::cerr << "all" << std::endl;
    used[1] = true;
    used[4] = true;
    used[3] = true;
    used[2] = true;
    if( !samples.fit_with(F,aorg,used,aerr) )
    {
        throw exception("couldn't find all");
    }
    GLS<double>::display(std::cerr,aorg,aerr);
    {
        ios::wcstream fp("fit3.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%.15g %.15g %.15g\n", u[i], Omega[i], OmegaFit[i]);
        }
    }
    Li.saveFit("omfit3.dat",u,aorg);



    return 0;


}
YOCTO_PROGRAM_END()



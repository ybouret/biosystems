
#include "yocto/program.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fit/glsf.hpp"
#include "yocto/container/utils.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/math/dat/linear.hpp"

using namespace yocto;
using namespace math;

class Intra
{
public:
    inline Intra() throw()
    {
    }

    inline ~Intra() throw()
    {
    }

    inline double Fit(const double t, const array<double> &a )
    {
        const double A = a[1];
        const double B = a[2];
        return A * (1.0 - exp( -B*t ) );
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Intra);
};

YOCTO_PROGRAM_START()
{
    if(argc<=1)
    {
        throw exception("usage: %s intra.txt",program);
    }

    vector<double> t;
    vector<double> C;

    {
        data_set<double> ds;
        ds.use(1,t);
        ds.use(2,C);
        ios::icstream fp(argv[1]);
        ds.load(fp);
    }
    const size_t N = t.size();
    std::cerr << "#data=" << N << std::endl;
    std::cerr << "t   ="  << t << std::endl;
    std::cerr << "C   ="  << C << std::endl;

    vector<double> Cfit(N);

    Intra intake;
    GLS<double>::Function F( &intake, &Intra::Fit );
    GLS<double>::Samples samples(1);
    GLS<double>::Sample &sample = samples.append(t,C,Cfit);
    (void)sample;

    const size_t    nvar = 2;
    samples.prepare(nvar);
    vector<double>  aorg(nvar);
    vector<double>  aerr(nvar);
    vector<bool>    used(nvar,true);
    aorg[1] = 2*C[N];
    aorg[2] = log(2.0)/t[N];
    sample.computeD2(F,aorg);
    {
        ios::wcstream fp("intra0.txt");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g\n", t[i], C[i], Cfit[i]);
        }
    }

    used[1] = false;

    std::cerr << "initialize time scale" << std::endl;
    if( !samples.fit_with(F, aorg, used, aerr) )
    {
        throw exception("couldn't adjust time scale");
    }
    GLS<double>::display(std::cerr, aorg, aerr);
    {
        ios::wcstream fp("intra1.txt");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g\n", t[i], C[i], Cfit[i]);
        }
    }

    used[1] = true;
    if( !samples.fit_with(F, aorg, used, aerr) )
    {
        throw exception("couldn't adjust time scale");
    }
    GLS<double>::display(std::cerr, aorg, aerr);
    {
        ios::wcstream fp("intra2.txt");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g\n", t[i], C[i], Cfit[i]);
        }
    }
    std::cerr << "1/kl=" << aorg[2] << std::endl;
    std::cerr << "kl  =" << 1.0/aorg[2] << "/s" << std::endl;
}
YOCTO_PROGRAM_END()


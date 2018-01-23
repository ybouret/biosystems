#include "yocto/program.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;

static double dLi7_out = 0;



YOCTO_PROGRAM_START()
{
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

    vector<double> Omega(N);
    for(size_t i=1;i<=N;++i)
    {
        Omega[i] = (1e-3*(dLi7[i]-dLi7_out))/(1.0+1e-3*dLi7_out);
    }

    {
        ios::wcstream fp("omega.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%.15g %.15g %.15g\n", t[i], Omega[i], dLi7[i]);
        }
    }

}
YOCTO_PROGRAM_END()



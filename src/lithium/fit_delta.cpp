#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;

YOCTO_PROGRAM_START()
{

    if(argc<=3)
    {
        throw exception("usage: %s file colIndex longTime", program);
    }

    const string   fileName = argv[1];
    const size_t   colIndex = strconv::to_size(  argv[2],"colIndex");
    const double   longTime = strconv::to_double(argv[3],"longTime");

    vector<double> t;
    vector<double> delta7;

    {
        data_set<double> ds;
        ds.use(1,t);
        ds.use(colIndex,delta7);
        ios::icstream fp(fileName);
        ds.load(fp);
    }
    const size_t N = t.size();
    std::cerr << "-- loaded " << N << " data" << std::endl;

    vector<double> tL(N,as_capacity);
    vector<double> d7L(N,as_capacity);
    for(size_t i=1;i<=N;++i)
    {
        if(t[i]>=longTime)
        {
            tL.push_back(t[i]);
            d7L.push_back(delta7[i]);
        }
    }
    const size_t NL = tL.size();
    std::cerr << "-- #longTime = " << NL << std::endl;

    {
        ios::wcstream fp("lt.dat");
        for(size_t i=1;i<=NL;++i)
        {
            fp("%.15g %.15g\n", tL[i], d7L[i]);
        }
    }


}
YOCTO_PROGRAM_END()


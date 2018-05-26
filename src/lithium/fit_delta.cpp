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

    if(argc<=2)
    {
        throw exception("usage: %s file #deltaColumn", program);
    }

    const string   fileName = argv[1];
    const size_t   colIndex = strconv::to_size(argv[2],"#deltaColumn");
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


}
YOCTO_PROGRAM_END()


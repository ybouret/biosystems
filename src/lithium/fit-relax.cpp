
#include "y/program.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/math/io/data-set.hpp"
#include "y/math/fit/ls.hpp"
#include "y/ptr/shared.hpp"
#include "y/sort/unique.hpp"
#include "y/math/fit/vectors.hpp"
#include "y/math/stat/metrics.hpp"
#include "y/math/fcn/zfind.hpp"
#include "y/sort/sorted-sum.hpp"

using namespace upsylon;
using namespace math;

typedef Fit::Sample<double>     Sample;
typedef Fit::Samples<double>    Samples;
typedef vector<double>          Vector;
typedef array<double>           Array;


Y_PROGRAM_START()
{

    if(argc<=1)
    {
        return 0;
    }

    const string file_name = argv[1];
    Vector t;
    Vector d7;
    {
        ios::icstream fp(file_name);
        data_set<double> ds;
        ds.use(1,t);
        ds.use(2,d7);
        ds.load(fp);
    }
    const size_t np = t.size();
    Vector d7fit(np,0);

    Sample sample(t,d7,d7fit);
    



}
Y_PROGRAM_END()


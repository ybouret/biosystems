#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;

//static const double rho_s = 12.0192;

YOCTO_PROGRAM_START()
{

    const string workdir = "src/lithium/doc/";

    vector<double> t_full;      //! col 1
    vector<double> delta_full;  //! col 3


    {
        const string filename = workdir + "nhe1_delta7_full_15mM_37.txt";
        std::cerr << "Loading " << filename << std::endl;
        data_set<double> ds;
        ds.use(1,t_full);
        ds.use(3,delta_full);
        ios::icstream fp(filename);
        ds.load(fp);
        std::cerr << "...done" << std::endl;
    }

    vector<double> LiExtDR;
    vector<double> LiIntDR;
    vector<double> deltaDR;
    vector<double> stddevDR;
    {
        const string filename = workdir + "nhe1_reponse_60s_37.txt";
        std::cerr << "Loading " << filename << std::endl;
        data_set<double> ds;
        ds.use(1,LiExtDR);
        ds.use(2,LiIntDR);
        ds.use(3,deltaDR);
        ds.use(4,stddevDR);
        ios::icstream fp(filename);
        ds.load(fp);
        std::cerr << "...done" << std::endl;
    }

    vector<double> t_short;
    vector<double> Li_short;
    {
        const string filename = workdir + "nhe1_total_short_times_15mM_37.txt";
        std::cerr << "Loading " << filename << std::endl;
        data_set<double> ds;
        ds.use(1,t_short);
        ds.use(2,Li_short);
        ios::icstream fp(filename);
        ds.load(fp);
        std::cerr << "...done" << std::endl;
    }

}
YOCTO_PROGRAM_END()


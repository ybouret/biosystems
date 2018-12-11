#include "y/program.hpp"
#include "y/math/io/data-set.hpp"
#include "y/sequence/vector.hpp"
#include "y/ios/ocstream.hpp"

using namespace upsylon;
using namespace math;

Y_PROGRAM_START()
{
    const double Vcell = 2600 * ipower(1e-6,3) * 1000; //!< in liters
    const double Ncell = 3.5e6;
    const double Vtot  = Vcell*Ncell;
    const double Vwork = 2e-3; //! in L

    std::cerr << "Vtot=" << Vtot << " L" << std::endl;

    size_t iCol = 2;

    vector<double> t;  // in secs
    vector<double> Li; // in muM

    if(argc>1)
    {
        const string filename = argv[1];
        {
            math::data_set<double> ds;
            ds.use(1, t);
            ds.use(iCol, Li);
            ds.load(filename);
        }

        const size_t n=t.size();
        for(size_t i=1;i<=n;++i)
        {
            Li[i] = (Li[i]*1e-6*Vwork)/Vtot;
            std::cerr << "t=" << t[i] << " => " << Li[i] << " M" << std::endl;
        }

        {
            const string  outname = "li-all-out.txt";
            ios::ocstream fp(outname);
            fp << "#t Li(M)\n";
            for(size_t i=1;i<=n;++i)
            {
                fp("%.15g %.15g\n",t[i],Li[i]);
            }
        }
    }


}
Y_PROGRAM_END()

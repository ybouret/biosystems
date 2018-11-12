#include "y/program.hpp"
#include "y/math/io/data-set.hpp"
#include "y/sequence/vector.hpp"
#include "y/string/convert.hpp"
#include "y/string/tokenizer.hpp"
#include "y/ios/icstream.hpp"
#include "y/ios/ocstream.hpp"

using namespace upsylon;
using namespace math;

static inline bool is_sep(const char C)
{
    return ' '==C || '\t' == C;
}

Y_PROGRAM_START()
{
    vector<string,memory::pooled> words;
    for(int iarg=1;iarg<argc;++iarg)
    {
        // loading data
        const string  filename = argv[iarg];
        std::cerr << std::endl << filename << std::endl;
        vector<double> t;
        vector<double> pH;
        vector<double> sem;
        {
            ios::icstream fp(filename);
            string line;
            if(!fp.gets(line))
            {
                throw exception("missing header in '%s'", *filename);
            }
            words.free();
            if(3!=tokenizer<char>::split(words,line,is_sep))
            {
                throw exception("invalid header in '%s'", *filename);
            }
            std::cerr << "header=" << words << std::endl;
            if("#t"!=words[1]) throw exception("invalid first header field in '%s'", *filename);
            const double pH_asymp = string_convert::to<double>(words[2],"pH");
            std::cerr << "|_pH_asymp=" << pH_asymp << std::endl;
            const double pH_asymp_sem = string_convert::to<double>(words[2],"pH");
            std::cerr << "|_  |_sem =" << pH_asymp_sem << std::endl;

            data_set<double> ds;
            ds.use(1, t);
            ds.use(2, pH);
            ds.use(3, sem);
            ds.load(fp);
        }
        size_t n = t.size();
        std::cerr << "#data="  << n << std::endl;

        {
            string outname = filename + ".dat";
            ios::ocstream fp(outname);
            for(size_t i=1;i<=n;++i)
            {
                fp("%.15g %.15g %.15g\n", t[i], pH[i], sem[i] );
            }
        }

    }
}
Y_PROGRAM_END()


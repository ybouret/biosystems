#include "y/program.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/sequence/vector.hpp"
#include "y/math/io/data-set.hpp"
#include "y/math/fit/ls.hpp"
#include "y/ptr/shared.hpp"

using namespace upsylon;
using namespace math;

typedef Fit::Sample<double>     Sample;
typedef Fit::Samples<double>    Samples;
typedef vector<double>          Vector;
typedef shared_ptr< Vector  >   VectorPtr;

Y_PROGRAM_START()
{

    vector<string>    files;
    vector<double>    concs;
    vector<VectorPtr> vectors;

    // finding files parameters
    std::cerr << "-- Parsing Arguments" << std::endl;
    {
        Lang::MatchString     match( "[:digit:]+" );
        vector<Lang::Token>   tokens(4,as_capacity);

        for(int i=1;i<argc;++i)
        {
            const string fn = argv[i];
            if( match(tokens,fn) <= 0 )
            {
                throw exception("no concentration found in '%s'", *fn );
            }
            const string first_C = tokens.front().to_string();
            std::cerr << " |_found '" << first_C << "'" << std::endl;

            files.push_back( fn );
            concs.push_back( string_convert::to<double>(first_C,"concentration") );
        }
    }

    std::cerr << " |_Processing " << files << ", with C=" << concs << std::endl;

    std::cerr << "-- Loading Files and Building Samples" << std::endl;
    const size_t         ns = files.size();
    if(ns<=0)
    {
        return 0;
    }

    vectors.ensure(3*ns);
    Samples samples(4,ns);
    for(size_t i=1;i<=ns;++i)
    {
        // prepare X/Y
        VectorPtr pX = new Vector();
        VectorPtr pY = new Vector();
        data_set<double> ds;
        ds.use(1, *pX);
        ds.use(2, *pY);
        {
            ios::icstream fp( files[i] );
            ds.load(fp);
        }
        // prepare Z with same size
        const size_t nd = pX->size();
        VectorPtr pZ = new Vector(nd);
        vectors.push_back(pX);
        vectors.push_back(pY);
        vectors.push_back(pZ);
        std::cerr << " |_Loaded '" << files[i] << "', #=" << nd << std::endl;

        // new sample
        samples.add(*pX, *pY, *pZ);
    }

    std::cerr << "-- Preparing fit functions" << std::endl;
    for(size_t i=1;i<=ns;++i)
    {
        std::cerr << " |_samples[" << i << "] : #=" << samples[i]->count() << std::endl;
    }
    


}
Y_PROGRAM_END()


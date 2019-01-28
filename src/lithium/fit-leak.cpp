#include "y/program.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/sequence/vector.hpp"
#include "y/math/io/data-set.hpp"
#include "y/math/fit/ls.hpp"
#include "y/ptr/shared.hpp"
#include "y/sort/unique.hpp"

using namespace upsylon;
using namespace math;

typedef Fit::Sample<double>     Sample;
typedef Fit::Samples<double>    Samples;
typedef vector<double>          Vector;
typedef shared_ptr< Vector  >   VectorPtr;
typedef array<double>           Array;
typedef Fit::Variables          Variables;


static const double beta_s = 12.0192;
static const double d7out  = 14.57;
static const double eps6   = 1.0/(1.0+beta_s*(1.0+0.001*d7out));
static const double eps7   = 1.0 - eps6;
static const double sigma  = 1.0/0.99772;

static  double grow6(const double tau) { return 1.0-exp(-sigma*tau); }
static  double grow7(const double tau) { return 1.0-exp(tau);        }
static  double grow(const double tau)  { return eps6 * grow6(tau) + eps7 * grow7(tau); }

class Leak
{
public:

    explicit Leak() 
    {
    }

    virtual ~Leak() throw()
    {
    }

    double Compute(double           t,
                   const Array     &aorg,
                   const Variables &vars)
    {
        const double Theta = vars(aorg,"Theta");
        const double k7    = vars(aorg,"k7");
        const double LiOut = vars(aorg,"Li");

        return LiOut * Theta * grow( k7 * t );

    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Leak);
};

Y_PROGRAM_START()
{

    vector<string>    files;
    vector<string>    labels;
    vector<double>    concs;
    vector<VectorPtr> vectors;

    // finding files parameters
    std::cerr << "-- Parsing Arguments" << std::endl;
    {
        Lang::MatchString     match_C( "[:digit:]+" );
        Lang::MatchString     match_L( "_[:word:]+" );
        vector<Lang::Token>   tokens(4,as_capacity);

        for(int i=1;i<argc;++i)
        {
            const string fn = argv[i];
            if( match_C(tokens,fn) <= 0 )
            {
                throw exception("no concentration found in '%s'", *fn );
            }
            const string first_C = tokens.front().to_string();
            std::cerr << " |_found '" << first_C << "'" << std::endl;

            if( match_L(tokens,fn) <= 0 )
            {
                throw exception("no acceptable label found in '%s'", *fn );
            }
            const string first_L = tokens.front().to_string(1,0);
            std::cerr << " |_found '" << first_L << "'" << std::endl;

            files.push_back( fn );
            labels.push_back( first_L );
            concs.push_back( string_convert::to<double>(first_C,"concentration") );
        }
    }

    vector<string> unique_labels( labels );
    unique(unique_labels);


    std::cerr << " |_Processing " << files << std::endl;
    std::cerr << "  |_with C=" << concs << std::endl;
    std::cerr << "  |_labels=" << unique_labels << std::endl;

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

    Variables &vars = samples.variables;
    vars << "k7"; //! global k7

    std::cerr << "vars=" << vars << std::endl;


    for(size_t i=1;i<=ns;++i)
    {
        std::cerr << " |_samples[" << i << "] : #=" << samples[i]->count() << std::endl;

        // individual fit

    }




}
Y_PROGRAM_END()


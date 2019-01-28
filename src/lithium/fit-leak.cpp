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


#if 0
static const double beta_s = 12.0192;
static const double d7out  = 14.57;
static const double eps6   = 1.0/(1.0+beta_s*(1.0+0.001*d7out));
#else
static const double eps6 = 0;
#endif

static const double eps7   = 1.0 - eps6;
static const double sigma  = 1.0/0.99772;

static  double grow6(const double tau) { return 1.0-exp(-sigma*tau); }
static  double grow7(const double tau) { return 1.0-exp(-tau);        }
static  double grow( const double tau) { return eps6 * grow6(tau) + eps7 * grow7(tau); }

static double grow6_linear(const double tau) { return sigma*tau; }
static double grow7_linear(const double tau) { return       tau; }
static double grow_linear(const double tau)  { return eps6 * grow6_linear(tau) + eps7 * grow7_linear(tau); }

static double grow6_sq(const double tau) { return sigma*tau-square_of(sigma*tau)/2; }
static double grow7_sq(const double tau) { return       tau-square_of(tau)/2; }
static double grow_sq(const double tau)  { return eps6 * grow6_sq(tau) + eps7 * grow7_sq(tau); }


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

    double ComputeLinear(double           t,
                         const Array     &aorg,
                         const Variables &vars)
    {
        const double Theta = vars(aorg,"Theta");
        const double k7    = vars(aorg,"k7");
        const double LiOut = vars(aorg,"Li");

        return LiOut * Theta * grow_linear( k7 * t );
    }

    double ComputeSq(double           t,
                     const Array     &aorg,
                     const Variables &vars)
    {
        const double Theta = vars(aorg,"Theta");
        const double k7    = vars(aorg,"k7");
        const double LiOut = vars(aorg,"Li");

        return LiOut * Theta * grow_sq( k7 * t );
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


    std::cerr << "-- Preparing fit variables" << std::endl;

    Variables &vars = samples.variables;
    vars << "k7"; // global k7

    bool useSameTheta      = false;
    bool useDifferentTheta = !useSameTheta;

    // different Theta's
    if(useSameTheta)
    {
        for(size_t i=1;i<=unique_labels.size();++i)
        {
            {
                const string th = "Theta_" + unique_labels[i];
            vars << th;
            }
        }
    }


    // different Li
    for(size_t i=1;i<=ns;++i)
    {
        const string Li = vformat("Li#%u", unsigned(i) );
        vars << Li;
        if(useDifferentTheta)
        {
            const string Th = vformat("Theta#%u",unsigned(i));
            vars << Th;
        }
    }



    std::cerr << "vars=" << vars << std::endl;

    const size_t nv = vars.size();
    Vector       aorg(nv,0);
    Vector       aerr(nv,0);
    vector<bool> used(nv,false);


    for(size_t i=1;i<=ns;++i)
    {
        std::cerr << " |_samples[" << i << "] : #=" << samples[i]->count() << std::endl;

        Variables &local = samples[i]->variables;
        local("k7", vars["k7"]);

        if(useSameTheta)
        {
            const string Th = "Theta_" + labels[i];
            local("Theta",vars[Th]);
            local(aorg,"Theta") = i;
        }
        else
        {
            const string Th = vformat("Theta#%u",unsigned(i));
            local("Theta",vars[Th]);
            local(aorg,"Theta") = i;
        }
        local(used,"Theta") = true;

        {
            const string Li = vformat("Li#%u", unsigned(i) );
            local("Li",vars[Li]);
            local(aorg,"Li") = concs[i];
        }

        std::cerr << "  |_local[" << i << "]=" << local << std::endl;
    }

    vars(aorg,"k7") = 0.002;

    std::cerr << "aorg=" << aorg << std::endl;


    // preparing function
    Leak                                leak;
    Fit::LeastSquares<double>           LS;

    vars(used,"k7")=true;

    {
        std::cerr << "== Fit Linear ==" << std::endl;
        Fit::LeastSquares<double>::Function F( &leak, & Leak::ComputeLinear );
        if( !LS.fit(samples, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit");
        }
        std::cerr << "-- Results: " << std::endl;
        const double R2 = samples.computeR2();
        std::cerr << "R2=" << R2 << std::endl;

        vars.display(std::cerr,aorg,aerr);
    }

    {
        std::cerr << "== Fit Quadratic ==" << std::endl;
        Fit::LeastSquares<double>::Function F( &leak, & Leak::ComputeSq );
        if( !LS.fit(samples, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit");
        }
        std::cerr << "-- Results: " << std::endl;
        const double R2 = samples.computeR2();
        std::cerr << "R2=" << R2 << std::endl;

        vars.display(std::cerr,aorg,aerr);
    }

    {
        std::cerr << "== Fit Full ==" << std::endl;
        Fit::LeastSquares<double>::Function F( &leak, & Leak::Compute );
        if( !LS.fit(samples, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit");
        }
        std::cerr << "-- Results: " << std::endl;
        const double R2 = samples.computeR2();
        std::cerr << "R2=" << R2 << std::endl;

        vars.display(std::cerr,aorg,aerr);
    }




    for(size_t i=1;i<=ns;++i)
    {
        std::cerr << "  |_Saving from " << files[i] << std::endl;
        string fitname = vfs::base_name_from(files[i]);
        vfs::change_extension(fitname, "fit.dat");
        std::cerr << "   |_Saving [" << fitname << "]" << std::endl;
        ios::ocstream fp(fitname);
        const Sample &s = *samples[i];
        for(size_t i=1;i<=s.X.size();++i)
        {
            fp("%.15g %.15g %.15g\n",s.X[i],s.Y[i],s.Yf[i]);
        }
    }



}
Y_PROGRAM_END()


#include "y/program.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/math/io/data-set.hpp"
#include "y/math/fit/ls.hpp"
#include "y/ptr/shared.hpp"
#include "y/sort/unique.hpp"
#include "y/math/fit/vectors.hpp"

using namespace upsylon;
using namespace math;

typedef Fit::Sample<double>     Sample;
typedef Fit::Samples<double>    Samples;
typedef vector<double>          Vector;
typedef array<double>           Array;
typedef Fit::Variables          Variables;
typedef Fit::VectorsDB<double>  VecDB;
typedef Fit::Vectors<double>    Vectors;


static const double lambda_s = 12.0192;
static const double d7out  = 14.57;
static const double sigma  = 1.0/0.99772;
static const double eps6   = 1.0/(1.0+lambda_s*(1.0+1.0e-3*d7out));
static const double eps7   = 1.0-eps6;

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
        const double k7      = vars(aorg,"k7");
        const double Theta   = vars(aorg,"Theta");
        const double Li      = vars(aorg,"Li");
        const double tau     = t*k7;

        return Li * Theta * ( eps6*(1.0-exp(-sigma*tau)) + eps7*(1.0-exp(-tau)));
    }




private:
    Y_DISABLE_COPY_AND_ASSIGN(Leak);
};

Y_PROGRAM_START()
{

    vector<string>         files;
    vector<string>         labels;
    vector<string>         concsID;
    vector<double>         concs;
    VecDB                  vdb( argc );
    Samples                samples(16,argc);


    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Processing and Loading files
    //
    //
    ////////////////////////////////////////////////////////////////////////////

    std::cerr << "-- Parsing Arguments" << std::endl;
    {
        Lang::MatchString     match_C( "[:digit:]+" );
        Lang::MatchString     match_L( "_[:word:]+" );
        vector<Lang::Token>   tokens(4,as_capacity);

        for(int i=1;i<argc;++i)
        {
            ////////////////////////////////////////////////////////////////////
            // Parsing file names
            ////////////////////////////////////////////////////////////////////
            const string fn = argv[i];
            std::cerr << " |_using '" << fn << "'" << std::endl;
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
            concsID.push_back(first_C);
            concs.push_back( string_convert::to<double>(first_C,"concentration") );

            ////////////////////////////////////////////////////////////////////
            // Loading data and making a sample out of it
            ////////////////////////////////////////////////////////////////////
            std::cerr << " |_loading '" << fn << "'" << std::endl;
            {
                Vectors &vecs = vdb.create(fn);
                {
                    data_set<double> ds;
                    ds.use(1, vecs.X);
                    ds.use(2, vecs.Y);
                    ios::icstream fp( fn);
                    ds.load(fp);
                }
                (void) vecs.add_to(samples);
            }
            std::cerr << "  \\_done" << std::endl;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Building global parameters
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    std::cerr << "concsID=" << concsID << std::endl;
    std::cerr << "labels =" << labels  << std::endl;



    vector<string> cell_types( labels );
    unique(cell_types);
    std::cerr << "cell_types=" << cell_types << std::endl;

    vector<string> li_values( concsID );
    unique(li_values);
    std::cerr << "li_values=" << li_values << std::endl;


    const size_t ns   = samples.size();
    Variables   &vars = samples.variables;

    // preparing variables
    vars << "Theta";
    for(size_t i=1;i<=ns;++i)
    {
        Variables &local = samples[i]->variables;
        local("Theta",vars);

        {
            const string Li = vformat("Li#%u", unsigned(i));
            vars << Li;
            local("Li",vars[Li]);
        }

        {
            const string k7 = vformat("k7#%u", unsigned(i));
            vars << k7;
            local("k7",vars[k7]);
        }
    }

    {
        std::cerr << "vars=" << vars << std::endl;
    }

    const size_t nv = vars.size();
    Vector       aorg(nv,0);
    Vector       aerr(nv,0);
    vector<bool> used(nv,false);

    // setting initial variables
    vars(aorg,"Theta") = 4.47;


#if 0
    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Building global parameters
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    const size_t ns   = samples.size();
    Variables   &vars = samples.variables;
    vars << "k";  // common leak rate


    vector<string> species( labels );
    unique(species);


    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Building parameters for samples
    //
    //
    ////////////////////////////////////////////////////////////////////////////

    for(size_t i=1;i<=ns;++i)
    {
        const string Ai = "A" + vformat("%u",unsigned(i));
        vars << Ai;
    }

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Allocating parameters
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    std::cerr << "vars=" << vars << std::endl;

    const size_t nv = vars.size();
    Vector       aorg(nv,0);
    Vector       aerr(nv,0);
    vector<bool> used(nv,false);


    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // initializing
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    const double k0 = 0.001;
    vars(aorg,"k") = k0;
    vars(used,"k") = true;

    for(size_t i=1;i<=ns;++i)
    {
        std::cerr << "-- Building sample#" << i << "=" << files[i] << std::endl;
        Sample    &s     = *samples[i];
        Variables &local = s.variables;

        // link to master coefficient
        local("k",vars);

        // setting amplitudes
        {
            const string Ai = "A" + vformat("%u",unsigned(i));
            local("A",vars[Ai]);
            local(aorg,"A") = s.slope()/k0;
            local(used,"A") = true;
        }


        std::cerr << "local=" << local << std::endl;
    }

    std::cerr << std::endl;
    vars.display(std::cerr,aorg,"\t(set) ");
    std::cerr << std::endl;
    vars.display(std::cerr,used, "\t(use) ");


    // preparing function
    Leak                                leak;
    Fit::LeastSquares<double>           LS;

    LS.verbose = true;

    Fit::LeastSquares<double>::Function F( &leak, & Leak::Compute );




    {
        std::cerr << "== Fit Linear ==" << std::endl;
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
        std::cerr << "  |_Saving [" << fitname << "]" << std::endl;
        ios::ocstream fp(fitname);
        const Sample &s = *samples[i];
        for(size_t i=1;i<=s.X.size();++i)
        {
            fp("%.15g %.15g %.15g\n",s.X[i],s.Y[i],s.Yf[i]);
        }
    }
#endif

}
Y_PROGRAM_END()


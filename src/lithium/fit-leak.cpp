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
        const double k      = vars(aorg,"k");
        const double Theta0 = vars(aorg,"Theta0");
        const double u2     = square_of(vars(aorg,"u"));
        const double LiOut  = vars(aorg,"LiOut");

        const double ThLi    = Theta0 * LiOut;
        const double u2ThLi   = u2 * ThLi;
        const double u2ThLip1 = 1.0 + u2ThLi;

        return (ThLi/u2ThLip1) * ( 1.0 - exp( -k * u2ThLip1 * t ) );
    }




private:
    Y_DISABLE_COPY_AND_ASSIGN(Leak);
};

Y_PROGRAM_START()
{

    vector<string>         files;
    vector<string>         labels;
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
    const size_t ns   = samples.size();
    Variables   &vars = samples.variables;
    vars << "k";       // common leak rate
    vars << "Theta0";  // common steady-state GHK level


    vector<string> species( labels );
    unique(species);

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Building parameters for species
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    vector<string> u_per_species(species.size(),as_capacity);
    for(size_t i=1;i<=species.size();++i)
    {
        const string u = "u_" + species[i];
        vars << u;
        u_per_species.push_back(u);
    }

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Building parameters for samples
    //
    //
    ////////////////////////////////////////////////////////////////////////////

    for(size_t i=1;i<=ns;++i)
    {
        const string Li = "Li#" + vformat("%u", unsigned(i));
        vars << Li;
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

    vars(aorg,"k")      = 0.001;
    vars(aorg,"Theta0") = 4.47;
    for(size_t i=1;i<=species.size();++i)
    {
        const string u = "u_" + species[i];
        vars(used,u) = false;
    }

    for(size_t i=1;i<=ns;++i)
    {
        std::cerr << "-- Building sample#" << i << "=" << files[i] << std::endl;
        Sample    &s     = *samples[i];
        Variables &local = s.variables;
        local("k",vars);
        local("Theta0",vars);
        // setting outside concentration
        {
            const string Li = "Li#" + vformat("%u", unsigned(i));
            local("LiOut", vars[Li]);
            local(aorg,"LiOut") = concs[i];
            local(used,"LiOut") = false;
        }

        // setting shielding parameters
        {
            const string u = "u_" + labels[i];
            local("u",vars[u]);
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



    vars(used,"k") = true;

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

    vars(used,u_per_species[1]) = true; vars(aorg,u_per_species[1]) = 0.01;
    vars(used,u_per_species[2]) = true; vars(aorg,u_per_species[2]) = 0.01;

    {
        std::cerr << "== Fit Linear 1 ==" << std::endl;
        if( !LS.fit(samples, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit");
        }
        std::cerr << "-- Results: " << std::endl;
        const double R2 = samples.computeR2();
        std::cerr << "R2=" << R2 << std::endl;

        vars.display(std::cerr,aorg,aerr);
    }

    vars(used,"Theta0") = true;
    {
        std::cerr << "== Fit Linear 2 ==" << std::endl;
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


}
Y_PROGRAM_END()


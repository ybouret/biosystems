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
static const double d7out    = 14.57;
static const double sigma0   = 1.0/0.99772;
static const double eps6     = 1.0/(1.0+lambda_s*(1.0+1.0e-3*d7out));
static const double eps7     = 1.0-eps6;

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
        const double sigma   = vars(aorg,"sigma");
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
    // Building global and local parameters
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
    vars << "sigma";

    if(true)
    {
        for(size_t i=1;i<=ns;++i)
        {
            Variables &local = samples[i]->variables;
            local("Theta",vars);
            local("sigma",vars);
            local("Li",vars,i);
            local("k7",vars,i);
        }
    }
    else
    {
        for(size_t i=1;i<=cell_types.size();++i)
        {
            vars << ( "k7" + cell_types[i] );
        }

        for(size_t i=1;i<=ns;++i)
        {
            Variables &local = samples[i]->variables;
            local("Theta",vars);
            local("Li",vars,i);
            const string id = "k7" + labels[i];
            local("k7",vars[id]);
        }
    }

    {
        std::cerr << "vars=" << vars << std::endl;
        for(size_t i=1;i<=ns;++i)
        {
            std::cerr << "|_sub#" << i << "=" << samples[i]->variables << std::endl;
        }
    }

    const size_t nv = vars.size();
    Vector       aorg(nv,0);
    Vector       aerr(nv,0);
    vector<bool> used(nv,false);

    // setting initial variables
    vars(aorg,"Theta") = 4.47;
    vars(aorg,"sigma") = sigma0;
    
    for(size_t i=1;i<=ns;++i)
    {
        const Variables &local = samples[i]->variables;
        local(aorg,"Li") = concs[i];
        local(aorg,"k7") = 0.001;
        local.on(used,"k7");
    }

    std::cerr << "values: " << std::endl;
    vars.display(std::cerr,aorg);
    std::cerr << "used:   " << std::endl;
    vars.display(std::cerr,used);

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // Building LeastSquares context
    //
    //
    ////////////////////////////////////////////////////////////////////////////

    Leak                                leak;
    Fit::LeastSquares<double>           LS;

    //LS.verbose = true;

    Fit::LeastSquares<double>::Function F( &leak, & Leak::Compute );

    int level = 0;
    {
        std::cerr << "== Fit == Level " << ++level << std::endl;
        if( !LS.fit(samples, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit");
        }
        std::cerr << "-- Results: " << std::endl;
        const double R2 = samples.computeR2();
        std::cerr << "R2=" << R2 << std::endl;

        vars.display(std::cerr,aorg,aerr);
    }

    // relax Theta
    if(false)
    {
        vars.on(used, "Theta");
        std::cerr << "== Fit == Level " << ++level << std::endl;
        if( !LS.fit(samples, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit");
        }
        std::cerr << "-- Results: " << std::endl;
        const double R2 = samples.computeR2();
        std::cerr << "R2=" << R2 << std::endl;

        vars.display(std::cerr,aorg,aerr);
    }


    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // save results
    //
    //
    ////////////////////////////////////////////////////////////////////////////
    ios::ocstream lf("fit-leaks.dat");
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
        const Variables &V = samples[i]->variables;
        lf << fitname << '\n';
        lf("\tLi    = %.15g\n", V(aorg,"Li") );
        lf("\tTheta = %.15g\n", V(aorg,"Theta"));
        lf("\tk7    = %.15g +/- %.15g\n", V(aorg,"k7"), V(aerr,"k7") );
        lf("\t#x60:\n");
        lf("\tk7    = %.15g +/- %.15g\n", V(aorg,"k7")*60, V(aerr,"k7")*60 );
        //lf("\t#/60:\n");
        //lf("\tk7    = %.15e +/- %.15e\n", V(aorg,"k7")/60, V(aerr,"k7")/60 );

    }

    ////////////////////////////////////////////////////////////////////////////
    //
    //
    // save log_file
    //
    //
    ////////////////////////////////////////////////////////////////////////////


}
Y_PROGRAM_END()



#include "y/program.hpp"
#include "y/math/io/data-set.hpp"
#include "y/sequence/vector.hpp"
#include "y/string/convert.hpp"
#include "y/string/tokenizer.hpp"
#include "y/ios/icstream.hpp"
#include "y/ios/ocstream.hpp"
#include "y/math/signal/linear.hpp"
#include "y/math/utils.hpp"
#include "y/core/locate.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/math/fit/ls.hpp"


using namespace upsylon;
using namespace math;


class Proton
{
public:
    const double vmax;
    Proton(const double pH_asymp) : vmax( pH_asymp )
    {
    }

    ~Proton() throw()
    {
    }

    double Compute(double                t,
                   const array<double>  &aorg,
                   const Fit::Variables &vars)
    {
        const double vmin = vars(aorg,"vmin");
        const double t1   = vars(aorg,"t1");
        const double a    = vars(aorg,"a");
        const double t0   = vars(aorg,"t0");

        const double tt  = (t-t0)/t1;
        const double a2  = square_of(a/t1);
        const double tt2 = tt*tt;
        const double tt3 = tt2*tt;
        if(t<=t0)
        {
            return vmin;
        }
        else
        {
            return vmax + (vmin-vmax) * exp( - tt3 / (a2+tt2) );
        }

    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Proton);
};

static inline bool is_sep(const char C)
{
    return ' '==C || '\t' == C;
}

Y_PROGRAM_START()
{
    vector<string,memory::pooled> words;
    vector<double>                tz(4,as_capacity);
    Lang::Matching                match = "[:digit:]+([.][:digit:]*)?mM";
    list<Lang::Token>             concs;
    
    for(int iarg=1;iarg<argc;++iarg)
    {
        ////////////////////////////////////////////////////////////////////////
        //
        // loading data
        //
        ////////////////////////////////////////////////////////////////////////
        const string  filename = argv[iarg];
        std::cerr << std::endl << filename << std::endl;
        match.find(concs,filename);
        if( 1 != concs.size() )
        {
            throw exception("Need To find 1 concentration in file name '%s'", *filename );
        }
        const string conc_str = concs.front().to_string(0,2);
        const double conc     = string_convert::to<double>( conc_str, "[Li]" );
        std::cerr << "[Li]=" << conc << "mM" << std::endl;


        vector<double> t;
        vector<double> pH;
        vector<double> sem;
        double pH_asymp = 0;
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
            pH_asymp = string_convert::to<double>(words[2],"pH");
            std::cerr << "|_pH_asymp=" << pH_asymp << std::endl;
            const double pH_asymp_sem = string_convert::to<double>(words[3],"sem");
            std::cerr << "|_  |_sem =" << pH_asymp_sem << std::endl;

            data_set<double> ds;
            ds.use(1, t);
            ds.use(2, pH);
            ds.use(3, sem);
            ds.load(fp);
        }
        size_t n = t.size();
        std::cerr << "#data="  << n << std::endl;

        ////////////////////////////////////////////////////////////////////////
        //
        // find half rise
        //
        ////////////////////////////////////////////////////////////////////////
        double pH_min = pH[1];
        double t_min  = t[1];
        for(size_t i=2;i<n;++i)
        {
            const double tmp = pH[i];
            if(tmp<=pH_min)
            {
                pH_min = tmp;
                t_min  = t[i];
            }
        }

        std::cerr << "min=" << pH_min << std::endl;
        const double pH_mid = 0.5*(pH_min+pH_asymp);
        linear::zfind(tz, pH_mid, t, pH);
        std::cerr << "tz=" << tz << std::endl;
        if(tz.size()!=1)
        {
            throw exception("cannot find just 1 half-time for '%s'", *filename);
        }
        const double thalf = tz[1];
        size_t       ihalf = 0;
        (void)core::locate(thalf, *t, t.size(), comparison::increasing<double>, ihalf);
        std::cerr << "@ihalf=" << ihalf << std::endl;

        ////////////////////////////////////////////////////////////////////////
        //
        // build regular curve
        //
        ////////////////////////////////////////////////////////////////////////
        vector<double> Y( n );
        {
            double tmp = pH_mid;
            for(size_t i=ihalf;i>0;--i)
            {
                Y[i] = (tmp=min_of( pH[i], tmp ));
            }
        }

        {
            double tmp = pH_mid;
            for(size_t i=ihalf+1;i<=n;++i)
            {
                Y[i] = (tmp=max_of( pH[i], tmp ));
            }
        }

        {
            const string  outname = vformat("out%s.dat", *conc_str);
            ios::ocstream fp(outname);
            const double h_ini = pow(10.0,-pH_min);
            const double h_end = pow(10.0,-pH_asymp);
            for(size_t i=1;i<=n;++i)
            {
                //fp("%g %g %g\n", t[i], log( (pH_asymp-Y[i])/(pH_asymp-pH_min)), pH[i] );
                fp("%g %g\n",
                   t[i]-t_min,
                   ((pow(10.0,-Y[i])-h_ini)/(h_end-h_ini))
                    );
            }
        }

        ////////////////////////////////////////////////////////////////////////
        //
        // preparing fit
        //
        ////////////////////////////////////////////////////////////////////////
        Fit::LeastSquares<double>   ls;
        Proton                      proton(pH_asymp);
        Fit::Type<double>::Function F( & proton, & Proton::Compute );

        vector<double>       Z(n);
        Fit::Sample<double>  sample(t,Y,Z);
        Fit::Variables       &vars  = sample.variables;
        vars << "t1" << "a" << "vmin" << "t0";

        const size_t   nvar = vars.size();
        vector<double> aorg( nvar );
        vector<bool>   used( nvar, true );
        vector<double> aerr( nvar, 0 );

        vars(aorg,"vmin") = pH_min;

        double &t0 = vars(aorg,"t0");
        t0 = t_min;
        std::cerr << "t0=" << t0 << std::endl;
        vars(aorg,"t1")   = thalf-t0;
        vars(aorg,"a")    = 0.66535 * vars(aorg,"t1");


        (void)sample.computeD2(F,aorg);
        {
            const string  outname = vformat("fit%s.dat", *conc_str);
            ios::ocstream fp(outname);
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", t[i], Y[i], Z[i]);
            }
        }

        tao::ld(used,false);
        vars(used,"t1") = true;
       // vars(used, "a") = true;
       // vars(used, "vmin") = true;
        //vars(used, "t0")= true;
        //exit(1);
        ls.verbose = true;
        if(!ls.fit(sample, F, aorg, aerr, used) )
        {
            throw exception("cannot fit");
        }

        vars.diplay(std::cerr, aorg, aerr);
        {
            const string  outname = vformat("fit%s.dat", *conc_str);
            ios::ocstream fp(outname);
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", t[i], Y[i], Z[i]);
            }
        }


    }

}
Y_PROGRAM_END()


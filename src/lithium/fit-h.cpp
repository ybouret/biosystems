
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
#include "y/sort/index.hpp"

using namespace upsylon;
using namespace math;


class Proton
{
public:
    Proton()
    {
    }

    ~Proton() throw()
    {
    }

    double Compute(double                t,
                   const array<double>  &aorg,
                   const Fit::Variables &vars)
    {
        const double t0   = vars(aorg,"t0");
        const double vmin = vars(aorg,"vmin");
        const double vmax = vars(aorg,"vmax");

        if(t<=t0)
        {
            return vmin;
        }
        else
        {
            const double q = vars(aorg,"q");
            const double p = vars(aorg,"p");
            const double tt = (t-t0)/q;
            const double U  = pow(tt,p);
            return vmin + (vmax-vmin) * U/(1.0+U);
        }



    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Proton);
};

static inline bool is_sep(const char C)
{
    return ' '==C || '\t' == C;
}

static inline
void save_fit(const Fit::Sample<double> &sample,
              const string              &outname,
              const array<double>       &aorg,
              const array<double>       &aerr)
{
    std::cerr << std::endl;
    std::cerr << "----" << std::endl;
    sample.variables.diplay(std::cerr,aorg,aerr);
    ios::ocstream fp(outname);
    for(size_t i=1;i<=sample.count();++i)
    {
        fp("%.15g %.15f %.15g\n", sample.X[i], sample.Y[i], sample.Yf[i] );
    }
}


Y_PROGRAM_START()
{
    vector<string,memory::pooled> words;
    vector<double>                tz(4,as_capacity);
    Lang::Matching                match = "[:digit:]+([.][:digit:]*)?mM";
    list<Lang::Token>             concs;

    vector<double> C;
    vector<double> P;
    vector<double> Q;
    vector<double> pHmax;
    vector<double> pHmin;

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
        vector<double> H(n);
        for(size_t i=1;i<=n;++i)
        {
            H[i] =pow(10.0,-Y[i]);
        }
        const double Hini = pow(10.0,-pH_min);
        const double Hend = pow(10.0,-pH_asymp);

        {
            const string  outname = vformat("out%s.dat", *conc_str);
            ios::ocstream fp(outname);
            for(size_t i=1;i<=n;++i)
            {
                //fp("%g %g %g\n", t[i], log( (pH_asymp-Y[i])/(pH_asymp-pH_min)), pH[i] );
                fp("%g %g\n",
                   t[i]-t_min,
                   (H[i]-Hini)/(Hend-Hini)
                    );
            }
        }

        ////////////////////////////////////////////////////////////////////////
        //
        // preparing fit
        //
        ////////////////////////////////////////////////////////////////////////
        Fit::LeastSquares<double>   ls; //ls.verbose = true;
        Proton                      proton;
        Fit::Type<double>::Function F( & proton, & Proton::Compute );

        vector<double>       Z(n);
        Fit::Sample<double>  sample(t,H,Z);
        Fit::Variables       &vars  = sample.variables;
        vars << "p" << "q" << "vmin" << "vmax" << "t0";

        const size_t   nvar = vars.size();
        vector<double> aorg( nvar );
        vector<bool>   used( nvar, false );
        vector<double> aerr( nvar, 0 );

        vars(aorg,"vmin") = Hini;
        vars(aorg,"vmax") = Hend;

        vars(aorg,"vmax") = pow(10.0,-7.0);

        double &t0 = vars(aorg,"t0");
        t0 = t_min;
        std::cerr << "t0=" << t0 << std::endl;
        vars(aorg,"q")   = thalf-t0;
        vars(aorg,"p")   = 1.2;

        const string  outname = vformat("fit%s.dat", *conc_str);

        (void)sample.computeD2(F,aorg);
        save_fit(sample, outname, aorg, aerr);



        vars.on(used,"p");
        vars.on(used,"q");

        if(!ls.fit(sample, F, aorg, aerr, used) )
        {
            throw exception("cannot fit");
        }

        save_fit(sample, outname, aorg, aerr);

        vars.on(used,"vmin");
        if(!ls.fit(sample, F, aorg, aerr, used) )
        {
            throw exception("cannot fit");
        }
        save_fit(sample, outname, aorg, aerr);


#if 1
        vars.on(used,"t0");
        if(!ls.fit(sample, F, aorg, aerr, used) )
        {
            throw exception("cannot fit");
        }
        save_fit(sample, outname, aorg, aerr);
#endif

        if(!ls.fit(sample, F, aorg, aerr, used) )
        {
            throw exception("cannot fit");
        }
        save_fit(sample, outname, aorg, aerr);


        {
            const string resname = vformat("res%s.dat",*conc_str);
            ios::ocstream fp(resname);
            for(size_t i=1;i<=n;++i)
            {
                fp("%.15g %.15g %.15g\n",t[i],pH[i],-log10(sample.Yf[i]));
            }
        }

        C.push_back( conc );
        P.push_back( vars(aorg, "p" ) );
        Q.push_back( vars(aorg, "q" ) );
        pHmax.push_back( -log10( vars(aorg,"vmax" ) ) );
        pHmin.push_back( -log10( vars(aorg,"vmin" ) ) );
    }





    vector<size_t> idx( C.size() );
    indexing::make(idx, comparison::increasing<double>, C);
    indexing::rank(C,idx);
    indexing::rank(P,idx);
    indexing::rank(Q,idx);
    indexing::rank(pHmax,idx);
    indexing::rank(pHmin,idx);

    std::cerr << "C=" << C << std::endl;
    std::cerr << "P=" << P << std::endl;
    std::cerr << "Q=" << Q << std::endl;
    std::cerr << "pHmin=" << pHmin << std::endl;
    std::cerr << "pHmax=" << pHmax << std::endl;

    {
        ios::ocstream fp("fit_proton.log");
        fp << "#C P Q pHmin pHmax\n";
        for(size_t i=1;i<=C.size();++i)
        {
            fp("%.15g %.15g %.15g %.15g %.15g\n", C[i], P[i], Q[i], pHmin[i], pHmax[i]);
        }
    }

}
Y_PROGRAM_END()


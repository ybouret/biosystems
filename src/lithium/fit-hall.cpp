
#include "y/program.hpp"
#include "y/math/io/data-set.hpp"
#include "y/sequence/vector.hpp"
#include "y/string/convert.hpp"
#include "y/string/tokenizer.hpp"
#include "y/ios/icstream.hpp"
#include "y/ios/ocstream.hpp"
#include "y/ptr/counted.hpp"
#include "y/ptr/intr.hpp"
#include "y/associative/set.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/math/signal/linear.hpp"
#include "y/core/locate.hpp"
#include "y/math/fit/ls.hpp"

using namespace upsylon;
using namespace math;

typedef vector<double>       Vector;
typedef array<double>        Array;
typedef Fit::Sample<double>  Sample;
typedef Fit::Samples<double> Samples;
typedef Fit::Variables       Variables;



class Proton
{
public:
    Proton()
    {
    }

    ~Proton() throw()
    {
    }

    double Compute(double           t,
                   const Array     &aorg,
                   const Variables &vars)
    {
        const double t0 = vars(aorg,"t0");
        const double vi = vars(aorg,"vi");
        const double ve = vars(aorg,"ve");
        const double hi = pow(10.0,-vi);

        if(t<=t0)
        {
            return hi;
        }
        else
        {
            const double q = vars(aorg,"q");
            const double p = vars(aorg,"p");
            const double tt = (t-t0)/q;
            const double U  = pow(tt,p);
            const double he = pow(10.0,-ve);
            return hi + (he-hi) * U/(1.0+U);
        }



    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Proton);
};



class Record : public counted_object
{
public:
    typedef intr_ptr<string,Record>     Pointer;
    typedef set<string,Pointer>         DataBase;

    const string filename;
    const string name; //!< from concentration
    const double conc; //!< value
    double       pH_end;
    double       pH_SE_end;

    Vector         t;
    Vector         pH;
    Vector         pH_SE;
    Vector         Y;
    Vector         H;
    Vector         Z; //!< fitted, whatever
    size_t         n;

    double         pH_min;
    double         t_min;
    double         thalf;

    double         Hini;
    double         Hend;

    Sample::Pointer sample;
    Variables      &vars;

    inline const string & key() const throw() { return name; }

    static inline bool is_sep(const char C)
    {
        return ' '==C || '\t' == C;
    }

    inline Record(const char *id ) :
    filename(id),
    name( GetNameFromFilename() ),
    conc( string_convert::to<double>(name,"conc") ),
    pH_end(0),
    pH_SE_end(0),
    t(),
    pH(),
    pH_SE(),
    n(0),
    pH_min(0),
    t_min(0),
    thalf(0),
    Hini(0),
    Hend(0),
    sample( new Sample(t,H,Z) ),
    vars( sample->variables )
    {
        vector<string,memory::pooled> words;
        vector<double>                tz(4,as_capacity);


        ios::icstream fp(filename);
        string        line;

        ////////////////////////////////////////////////////////////////////////
        //
        // loading data from first line
        //
        ////////////////////////////////////////////////////////////////////////
        if(!fp.gets(line))
        {
            throw exception("missing header in '%s'", *filename);
        }
        if(3!=tokenizer<char>::split(words,line,is_sep))
        {
            throw exception("invalid header in '%s'", *filename);
        }
        std::cerr << "header=" << words << std::endl;
        if("#t"!=words[1]) throw exception("invalid first header field in '%s'", *filename);
        pH_end     = string_convert::to<double>(words[2],"pH_end"   );
        pH_SE_end  = string_convert::to<double>(words[3],"pH_SE_end");
        std::cerr << "pH_end=" << pH_end << " +/- " << pH_SE_end << std::endl;

        ////////////////////////////////////////////////////////////////////////
        //
        // loading pH
        //
        ////////////////////////////////////////////////////////////////////////
        data_set<double> ds;
        ds.use(1, t);
        ds.use(2, pH);
        ds.use(3, pH_SE);
        ds.load(fp);
        n = t.size();
        std::cerr << "..loaded #data=" << n << std::endl;
        Y.make(n,0);
        H.make(n,0);
        Z.make(n,0);

        ////////////////////////////////////////////////////////////////////////
        //
        // find half rise
        //
        ////////////////////////////////////////////////////////////////////////
        pH_min = pH[1];
        t_min  = t[1];
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
        const double pH_mid = 0.5*(pH_min+pH_end);
        linear::zfind(tz, pH_mid, t, pH);
        std::cerr << "tz=" << tz << std::endl;
        if(tz.size()!=1)
        {
            throw exception("cannot find just 1 half-time for '%s'", *filename);
        }
        thalf = tz[1];
        size_t ihalf = 0;
        (void)core::locate(thalf, *t, t.size(), comparison::increasing<double>, ihalf);
        std::cerr << "@ihalf=" << ihalf << std::endl;

        ////////////////////////////////////////////////////////////////////////
        //
        // build regular curve
        //
        ////////////////////////////////////////////////////////////////////////
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
                Y[i] = (tmp=max_of( pH[i], tmp )); // TODO: check this part...
            }
        }

        for(size_t i=1;i<=n;++i)
        {
            H[i] =pow(10.0,-Y[i]);
        }
        Hini = pow(10.0,-pH_min);
        Hend = pow(10.0,-pH_end);
    }

    inline string GetNameFromFilename() const
    {
        Lang::Matching    match = "[:digit:]+([.][:digit:]*)?mM";
        list<Lang::Token> concs;
        match.find(concs,filename);
        if( 1 != concs.size() )
        {
            throw exception("Need To find 1 concentration in file name '%s'", *filename );
        }
        return concs.front().to_string(0,2);
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Record);
};

typedef Record::DataBase::iterator Iterator;

static inline int compare_natural( const string &lhs, const string &rhs )
{
    const double l = atof(*lhs);
    const double r = atof(*rhs);
    return comparison::increasing(l,r);
}

Y_PROGRAM_START()
{
    Record::DataBase db;

    for(int i=1;i<argc;++i)
    {
        const Record::Pointer r = new Record(argv[i]);
        if(!db.insert(r))
        {
            throw exception("multiple file/conc '%s'",*(r->name));
        }
    }
    const size_t nrec = db.size();
    std::cerr << "Loaded #" << nrec << " files" << std::endl;
    db.sort_keys(compare_natural);

    Samples samples;
    for( Iterator i=db.begin();i!=db.end();++i)
    {
        Record &r = **i;
        std::cerr << r.filename << std::endl;
        samples.push_back(r.sample);
    }
    assert(samples.size()==nrec);

    ////////////////////////////////////////////////////////////////////////////
    //
    // prepare fit
    //
    ////////////////////////////////////////////////////////////////////////////
    Fit::LeastSquares<double>   ls;
    Proton                      proton;
    Fit::Type<double>::Function F( &proton, & Proton::Compute );

    ////////////////////////////////////////////////////////////////////////////
    //
    // First pass, individual fits
    //
    ////////////////////////////////////////////////////////////////////////////
    Variables &gvars = samples.variables;
    Vector C;
    Vector P;
    Vector Q;
    Vector pHi;
    Vector pHe;
    {
        gvars.free();
        gvars << "p" << "q" << "vi" << "ve" << "t0";
        const size_t nvar = gvars.size();
        Vector aorg( nvar );
        Vector aerr( nvar );
        vector<bool> used( nvar, false);

        {
            ios::ocstream fp("ind.dat");
            for( Iterator i=db.begin();i!=db.end();++i)
            {
                Record &r = **i;
                Sample &sample = *r.sample;
                r.vars = gvars;
                std::cerr << "-- Individual for " << r.filename << std::endl;
                std::cerr << r.name << ".vars=" << r.vars << std::endl;

                gvars(aorg,"p")  = 1.0;
                gvars(aorg,"q")  = r.thalf-r.t_min;
                gvars(aorg,"vi") = r.pH_min;
                gvars(aorg,"ve") = r.pH_end;

                gvars(aorg,"ve") = 7.1;


                tao::ld(used,false);


                int level = 0;
                gvars.diplay(std::cerr,aorg);
                {
                    //gvars.on(used,"p");
                    gvars.on(used,"q");
                    ++level;
                    std::cerr << "starting level " << level << " with: " << std::endl;
                    gvars.diplay(std::cerr,aorg);
                    if( !ls.fit( sample, F, aorg, aerr, used) )
                    {
                        throw exception("couldn't fit %s @level-%d", *(r.name), level );
                    }
                    //gvars.diplay(std::cerr, aorg, aerr, "\t");
                    //std::cerr << "\tR2=" << sample.computeR2() << std::endl;
                }

                if(false)
                {
                    gvars.on(used,"p");
                    ++level;
                    std::cerr << "starting level " << level << " with: " << std::endl;
                    gvars.diplay(std::cerr,aorg);
                    if( !ls.fit( sample, F, aorg, aerr, used) )
                    {
                        throw exception("couldn't fit %s @level-%d", *(r.name), level );
                    }
                }

                {
                    gvars.on(used,"t0");
                    gvars.on(used,"vi");
                    gvars.on(used,"ve");
                    ++level;
                    std::cerr << "starting level " << level << " with: " << std::endl;
                    gvars.diplay(std::cerr,aorg);
                    if( !ls.fit( sample, F, aorg, aerr, used) )
                    {
                        throw exception("couldn't fit %s @level-%d", *(r.name), level );
                    }
                }

                if(true)
                {
                    tao::ld(used,false);
                    gvars.on(used,"ve");
                    ++level;
                    std::cerr << "starting level " << level << " with: " << std::endl;
                    gvars.diplay(std::cerr,aorg);
                    if( !ls.fit( sample, F, aorg, aerr, used) )
                    {
                        throw exception("couldn't fit %s @level-%d", *(r.name), level );
                    }
                }

                if(true)
                {
                    tao::ld(used,false);
                    gvars.on(used,"p:q");
                    ++level;
                    std::cerr << "starting level " << level << " with: " << std::endl;
                    gvars.diplay(std::cerr,aorg);
                    if( !ls.fit( sample, F, aorg, aerr, used) )
                    {
                        throw exception("couldn't fit %s @level-%d", *(r.name), level );
                    }
                    //ls.verbose = false;
                }


                gvars.diplay(std::cerr, aorg, aerr, "\t");
                std::cerr << "\tR2=" << sample.computeR2() << std::endl;
                for(size_t j=1;j<=sample.count();++j)
                {
                    fp("%g %g %gn\n", sample.X[j], -log10(sample.Y[j]), -log10(sample.Yf[j]));
                }
                fp << '\n';

                C.push_back(   r.conc );
                P.push_back(   gvars(aorg, "p" ) );
                Q.push_back(   gvars(aorg, "q" ) );
                pHi.push_back( gvars(aorg, "vi" ) );
                pHe.push_back( gvars(aorg, "ve" ) );

            }
        }
    }

    std::cerr << "Global R2=" << samples.computeR2() << std::endl;
    for(size_t i=1;i<=nrec;++i)
    {
        std::cerr << C[i] << " " << P[i] << " " << Q[i] << " " << pHi[i] << " " << pHe[i] << std::endl;
    }


    ////////////////////////////////////////////////////////////////////////////
    //
    // Second pass, linked field
    //
    ////////////////////////////////////////////////////////////////////////////
    gvars.free();
    const double p_ini = average_of(P);
    for( Iterator i=db.begin();i!=db.end();++i)
    {
        Record    &r = **i;
        Variables &v = r.sample->variables;
        v.free();

    }


}
Y_PROGRAM_END()


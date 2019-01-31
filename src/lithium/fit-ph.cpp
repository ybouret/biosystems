#include "y/program.hpp"

#include "y/associative/set.hpp"
#include "y/math/fit/ls.hpp"
#include "y/ptr/intr.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/math/io/data-set.hpp"
#include "y/string/tokenizer.hpp"
#include "y/fs/vfs.hpp"

using namespace upsylon;
using namespace math;
using namespace Lang;

typedef vector<double>                Vector;
typedef vector<string,memory::pooled> Strings;


class Proton
{
public:
    explicit Proton() throw()
    {
    }

    virtual ~Proton() throw()
    {
    }

    double Compute( double t, const array<double> &a, const Fit::Variables &vars )
    {
        const double vi = vars(a,"Hini");
        const double ve = vars(a,"Hend");
        const double t0 = vars(a,"t0");
        if( t<= t0 )
        {
            return vi;
        }
        else
        {
            const double p = vars(a,"p");
            const double q = vars(a,"q");
            const double tp = pow(t-t0,p);
            const double qp = pow(q,p);
            return vi + (ve-vi) * tp/(tp+qp);
        }
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Proton);
};

static bool isWS( const char C) { return ' ' == C || '\t' == C; }

class Record : public counted_object
{
public:
    typedef intr_ptr<string,Record> Pointer;
    typedef set<string,Pointer>     DataBase;

    const string file_name;
    const string Li;           //!< concentration string, used as key
    Vector       t;
    Vector       pH;
    const size_t N;
    double       t0;
    double       pH_asymp;
    Vector       H;
    double       Hmin;
    double       Hmax;
    double       Hend;
    Vector       Y;
    string       save_name;
    Vector       Yf;
    Fit::Sample<double> &sample;

    inline const string & key() const throw() { return Li; }

    inline explicit Record( const string &fn, Fit::Samples<double> &samples ) :
    file_name( fn ),
    Li( GetLiFromString(file_name) ),
    t(),
    pH(),
    N(0),
    t0(0),
    pH_asymp(0),
    H(),
    Hmin(0),
    Hmax(0),
    Hend(0),
    save_name( vfs::base_name_from(fn) ),
    Yf(),
    sample( samples.add(t,Y,Yf) )
    {
        std::cerr << "|_Record: [" << Li << "] from '" << file_name << "'" << std::endl;
        vfs::change_extension(save_name,"dat");
        // Load data
        {
            data_set<double> ds;
            ds.use(1,t);
            ds.use(2,pH);
            ios::icstream fp(file_name);
            {
                // reading first line
                string line;
                if(!fp.gets(line)) throw exception("missing first line in '%s'", *file_name);
                Strings words;
                tokenizer<char>::split(words,line,isWS);
                if(words.size()<2) throw exception("invalid #tokens on first line");
                pH_asymp = string_convert::to<double>( words[2], "pH_asymp" );
            }
            ds.load(fp);
        }
        (size_t&)N = t.size();

        // computing H and metrics
        H.make(N,0);
        t0 = t[1];
        Hmin = Hmax = H[1] = pow(10,-pH[1]);
        Hend = pow(10,-pH_asymp);
        for(size_t i=2;i<=N;++i)
        {
            const double h = (H[i]=pow(10.0,-pH[i]));
            if(h>Hmax)
            {
                t0=t[i];
                Hmax=h;
            }
            else if(h<Hmin)
            {
                Hmin = h;
            }
            else
            {
                // do nothing
            }
        }

        // Building fitting curve
        Y.make(N,0);
        Yf.make(N,0);
        for(size_t i=1;i<=N;++i)
        {
            if(t[i]<=t0||i<=1)
            {
                Y[i] = Hmax;
            }
            else
            {
                assert(i>1);
                Y[i] = min_of( H[i], Y[i-1] );assert(Y[i]<=Y[i-1]);
                //Y[i] = H[i];
            }
        }

        std::cerr << "| \\_Loaded #" << N << std::endl;
        std::cerr << "| \\_pH_asymp = " << pH_asymp  << std::endl;
        std::cerr << "| \\_Hmin     = " << Hmin      << std::endl;
        std::cerr << "| \\_Hmax     = " << Hmax      << std::endl;
        std::cerr << "| \\_Hend     = " << Hend      << std::endl;
        std::cerr << "| \\_t0       = " << t0        << std::endl;
        std::cerr << "| \\_save     = " << save_name << std::endl;
        std::cerr << "|" << std::endl;

        save();
    }

    inline virtual ~Record() throw()
    {
    }

    inline static string GetLiFromString(const string &s )
    {
        MatchString                  match = "[:digit:]+([.][:digit:]+)?mM";
        vector<Token,memory::pooled> tokens;
        if( 1 != match(tokens,s) )
        {
            throw exception("need one concentration in file name '%s'", *s);
        }
        return tokens.front().to_string(0,2);
    }

    void save() const
    {
        ios::ocstream fp(save_name);
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g %g %g %g\n", t[i], H[i], Y[i], Yf[i], pH[i], -log10(Yf[i]));
        }
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Record);
};

static inline int compare_natural( const string &lhs, const string &rhs )
{
    const double l = atof(*lhs);
    const double r = atof(*rhs);
    return comparison::increasing(l,r);
}

class Records : public Record::DataBase
{
public:
    inline explicit Records() throw() {}
    inline virtual ~Records() throw() {}

    inline void add( const char *fn, Fit::Samples<double> &samples )
    {
        const string          _(fn);
        const Record::Pointer p = new Record(_,samples);
        if( !insert(p) )
        {
            throw exception("Multiple concentration %s", *(p->Li));
        }
        sort_keys(compare_natural);
    }

    void processParameters(const array<double> &aorg)
    {
        std::cerr << "Saving Parameters..." << std::endl;
        Vector Li(size(),as_capacity);
        {
            ios::ocstream fp("pH_params.dat");
            //fp("0 0\n");
            for(iterator i=begin();i!=end();++i)
            {
                const Record         &r    = **i;
                const Fit::Variables &vars = r.sample.variables;
                Li.push_back( string_convert::to<double>(r.Li,"Li") );
                const double Hini = vars(aorg,"Hini");
                const double Hend = vars(aorg,"Hend");
                const double dpH  = -log10(Hend) + log10(Hini);
                const double q    = vars(aorg,"q");
                fp("%.15g %.15g %.15g\n", Li.back(),dpH,q);
            }

        }
    }


private:
    Y_DISABLE_COPY_AND_ASSIGN(Records);
};

typedef Records::iterator Iterator;

Y_PROGRAM_START()
{
    Fit::Samples<double> samples;
    Records              db;
    for(int i=1;i<argc;++i)
    {
        db.add(argv[i],samples);
    }

    // building variables
    Fit::Variables &gvars = samples.variables;
    gvars << "p";// << "q";
    for(Iterator i=db.begin();i!=db.end();++i)
    {
        Record         &r = **i;
        Fit::Variables &lvars = r.sample.variables;
        std::cerr << "using " << r.file_name << std::endl;

        lvars("p",gvars);
       // lvars("q",gvars);

        static const char *vnames[] =
        {
            "q",
            "Hini",
            "Hend",
            "t0"
        };

        for(size_t j=0;j<sizeof(vnames)/sizeof(vnames[0]);++j)
        {
            const string name = vnames[j] + r.Li;
            gvars << name;
            lvars(vnames[j],gvars[name]);
        }

        std::cerr << "lvars@" << r.Li << "=" << lvars << std::endl;
    }

    std::cerr << "gvars=" << gvars << std::endl;

    // initializing variables
    const size_t nv = gvars.size();

    Vector       aorg( nv );
    Vector       aerr( nv );
    vector<bool> used(nv,false);

    gvars.on(used,"p");
    gvars(aorg,"p") = 1;

    for(Iterator i=db.begin();i!=db.end();++i)
    {
        Record         &r      = **i;
        Fit::Variables &lvars = r.sample.variables;

        lvars(aorg,"Hini") = r.Hmax;
        lvars(aorg,"Hend") = r.Hend;
        lvars(aorg,"t0"  ) = r.t0;
        lvars(aorg,"q")    = 30;
        lvars.on(used,"q");
    }

    gvars.display(std::cerr,aorg);

    Proton                      proton;
    Fit::Type<double>::Function F( &proton, & Proton::Compute );
    Fit::LeastSquares<double>   LS;

    int level=0;
    if( !LS.fit(samples, F, aorg, aerr, used) )
    {
        throw exception("couldn't fit level-%d",++level);
    }

    gvars.display(std::cerr,aorg,aerr);
    for(Iterator i=db.begin();i!=db.end();++i)
    {
        Record         &r      = **i;
        r.save();
    }


    if( true )
    {
        std::cerr << std::endl;
        for(Iterator i=db.begin();i!=db.end();++i)
        {
            Record         &r      = **i;
            Fit::Variables &lvars = r.sample.variables;
            lvars.on(used,"Hini:Hend");
        }
        if( !LS.fit(samples, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit level-%d",++level);
        }
        gvars.display(std::cerr,aorg,aerr);
        for(Iterator i=db.begin();i!=db.end();++i)
        {
            Record         &r      = **i;
            r.save();
        }
    }

    if( true )
    {
        std::cerr << std::endl;
        for(Iterator i=db.begin();i!=db.end();++i)
        {
            Record         &r      = **i;
            Fit::Variables &lvars = r.sample.variables;
            lvars.on(used,"t0");
        }
        if( !LS.fit(samples, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit level-%d",++level);
        }
        gvars.display(std::cerr,aorg,aerr);
        for(Iterator i=db.begin();i!=db.end();++i)
        {
            Record         &r      = **i;
            r.save();
        }
    }

    gvars.off(used,"p");
    gvars(aorg,"p") = 1;

    if( true )
    {
        if( !LS.fit(samples, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit level-%d",++level);
        }
        gvars.display(std::cerr,aorg,aerr);
        for(Iterator i=db.begin();i!=db.end();++i)
        {
            Record         &r      = **i;
            r.save();
        }
    }

    correlation<double> cctx;
    const double R2 = samples.computeR2();
    const double cr = samples.compute_correlation(cctx);
    std::cerr << std::endl;
    std::cerr << "-- R2=" << R2 << ", corr=" << cr << std::endl;
    db.processParameters(aorg);





}
Y_PROGRAM_END()


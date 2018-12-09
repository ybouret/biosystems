
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

static inline bool is_sep(const char C)
{
    return ' '==C || '\t' == C;
}

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

    double         Hini;
    double         Hend;

    Sample::Pointer sample;
    Variables      &vars;

    inline const string & key() const throw() { return name; }

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
        const double thalf = tz[1];
        size_t       ihalf = 0;
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
    db.sort_keys(comparison::increasing<string>);

    Samples samples;
    for( Iterator i=db.begin();i!=db.end();++i)
    {
        Record &r = **i;
        std::cerr << r.filename << std::endl;
        samples.push_back(r.sample);
    }
}
Y_PROGRAM_END()


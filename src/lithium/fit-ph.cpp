#include "y/program.hpp"

#include "y/associative/set.hpp"
#include "y/math/fit/ls.hpp"
#include "y/ptr/intr.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/math/io/data-set.hpp"
#include "y/string/tokenizer.hpp"

using namespace upsylon;
using namespace math;
using namespace Lang;

typedef vector<double>                Vector;
typedef vector<string,memory::pooled> Strings;

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

    inline const string & key() const throw() { return Li; }

    inline explicit Record( const string &fn ) :
    file_name( fn ),
    Li( GetLiFromString(file_name) ),
    t(),
    pH(),
    N(0),
    t0(0),
    pH_asymp(0)
    {
        std::cerr << "|_Record: [" << Li << "] from '" << file_name << "'" << std::endl;

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
        std::cerr << "| \\_Loaded #" << N << std::endl;
        std::cerr << "| \\_pH_asymp=" << pH_asymp << std::endl;
        std::cerr << "|" << std::endl;

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

    inline void add( const char *fn )
    {
        const string _(fn);
        const Record::Pointer p = new Record(_);
        if( !insert(p) )
        {
            throw exception("Multiple concentration %s", *(p->Li));
        }
        sort_keys(compare_natural);
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Records);
};

typedef Records::iterator Iterator;

Y_PROGRAM_START()
{
    Records db;
    for(int i=1;i<argc;++i)
    {
        db.add(argv[i]);
    }

    for(Iterator i=db.begin();i!=db.end();++i)
    {
        std::cerr << "using " << (**i).file_name << std::endl;
    }

}
Y_PROGRAM_END()


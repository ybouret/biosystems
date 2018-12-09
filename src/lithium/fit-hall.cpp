
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

using namespace upsylon;
using namespace math;

typedef vector<double> Vector;
typedef array<double> Array;


class Record : public counted_object
{
public:
    typedef intr_ptr<string,Record>     Pointer;
    typedef set<string,Pointer>         DataBase;

    const string filename;
    const string name; //!< from concentration
    const double conc; //!< value

    inline const string & key() const throw() { return name; }

    inline Record(const char *id ) :
    filename(id),
    name( GetNameFromFilename() ),
    conc( string_convert::to<double>(name,"conc") )
    {

    }

    string GetNameFromFilename() const
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
    for( Iterator i=db.begin();i!=db.end();++i)
    {
        std::cerr << (**i).filename << std::endl;
    }
}
Y_PROGRAM_END()


#include "y/program.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/sequence/vector.hpp"
#include "y/math/io/data-set.hpp"

using namespace upsylon;



Y_PROGRAM_START()
{

    vector<string> files;
    vector<double> concs;

    // finding files parameters
    {
        Lang::MatchString     match( "[:digit:]+" );
        vector<Lang::Token>   tokens(4,as_capacity);

        for(int i=1;i<argc;++i)
        {
            const string fn = argv[i];
            if( match(tokens,fn) <= 0 )
            {
                throw exception("no concentration found in '%s'", *fn );
            }
            const string first_C = tokens.front().to_string();
            std::cerr << "found '" << first_C << "'" << std::endl;

            files.push_back(fn);
            concs.push_back( string_convert::to<double>(first_C,"concentration") );
        }
    }

    std::cerr << "Processing " << files << ", with C=" << concs << std::endl;

}
Y_PROGRAM_END()


#include "y/program.hpp"
#include "y/math/io/data-set.hpp"
#include "y/sequence/vector.hpp"
#include "y/string/convert.hpp"
#include "y/string/tokenizer.hpp"
#include "y/ios/icstream.hpp"

using namespace upsylon;
using namespace math;


Y_PROGRAM_START()
{
    for(int iarg=1;iarg<argc;++iarg)
    {
        const string  filename = argv[iarg];
        ios::icstream fp(filename);

    }
}
Y_PROGRAM_END()


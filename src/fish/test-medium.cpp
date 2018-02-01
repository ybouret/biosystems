#define MEDIUM_TEST 1
#include "medium.h"
#include "yocto/program.hpp"
#include "yocto/code/alea.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;

typedef Medium::NodeOf<int>  Node;
typedef Medium::ListOf<Node> List;
typedef Medium::PoolOf<Node> Pool;

#define NN 1000
YOCTO_PROGRAM_START()
{
    alea.initialize();

    List l;
    Pool p;
    for(size_t i=0;i<NN;++i)
    {
        p.store( new Node() );
    }
    if(p.size!=NN) throw exception("Pool error");

    while(p.size>0)
    {
        if( alea.nextBool() )
        {
            l.push_back( p.query() );
        }
        else
        {
            l.push_front( p.query() );
        }
    }

    while(l.size>0)
    {
        if( alea.nextBool() ) delete l.pop_front(); else delete l.pop_back();
    }

    Medium medium;

    ios::wcstream fp("waves.dat");
    const double period=10;
    for(double x=-2*period;x<=2*period;x+=0.01)
    {
        fp("%g %g %g\n",x,Medium::TriangleWave(x,period),Medium::CosWave(x,period));
    }
}
YOCTO_PROGRAM_END()



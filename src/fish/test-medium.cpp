#define MEDIUM_TEST 1
#include "medium.h"
#include "y/program.hpp"
#include "y/alea.hpp"
#include "y/ios/ocstream.hpp"

using namespace upsylon;

typedef Medium::NodeOf<int>  Node;
typedef Medium::ListOf<Node> List;
typedef Medium::PoolOf<Node> Pool;

#define NN 1000
Y_PROGRAM_START()
{
    alea_init();

    List l;
    Pool p;
    for(size_t i=0;i<NN;++i)
    {
        p.store( new Node() );
    }
    if(p.size!=NN) throw exception("Pool error");

    while(p.size>0)
    {
        if( alea.choice() )
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
        if( alea.choice() ) delete l.pop_front(); else delete l.pop_back();
    }

    Medium medium;

    ios::ocstream fp("waves.dat");
    const double period=10;
    for(double x=-2*period;x<=2*period;x+=0.01)
    {
        fp("%g %g %g\n",x,Medium::TriangleWave(x,period),Medium::CosWave(x,period));
    }
}
Y_PROGRAM_END()



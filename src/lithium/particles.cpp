#include "yocto/math/point3d.hpp"
#include "yocto/program.hpp"
#include "yocto/core/list.hpp"
#include "yocto/threading/scheme/simd.hpp"
#include "yocto/counted-object.hpp"
#include "yocto/ptr/arc.hpp"
#include "yocto/sequence/vector.hpp"

using namespace yocto;

class Particle : public object
{
public:
    Particle *next;
    Particle *prev;
    inline explicit Particle() throw() : next(0), prev(0)
    {
    }

    inline virtual ~Particle() throw() {}

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Particle);
};

typedef core::list_of_cpp<Particle> Particles;

class Worker : public counted_object
{
public:
    typedef arc_ptr<Worker> Pointer;

    inline Worker() throw()
    {
    }

    inline ~Worker() throw()
    {
    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Worker);
};

typedef vector<Worker::Pointer> __Workers;
class Workers : public __Workers
{
public:

    inline explicit Workers(const size_t numWorkers): __Workers()
    {
        for(size_t i=0;i<numWorkers;++i)
        {
            const Worker::Pointer pW( new Worker() );
            push_back(pW);
        }
    }

    inline virtual ~Workers() throw()
    {

    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Workers);
};

class Simulation
{
public:
    Particles particles;  //!< original particles
    Particles inactive;   //!< if any
    Workers   workers;    //!< workers, one per thread

    inline explicit Simulation(const size_t numParticles,const size_t numThreads) :
    particles(),
    inactive(),
    workers(numThreads)
    {
        for(size_t i=0;i<numParticles;++i)
        {
            particles.push_back( new Particle() );
        }
        std::cerr << "[SIMULATION] #particles = " << particles.size << std::endl;
        std::cerr << "[SIMULATION] #workers   = " << workers.size() << std::endl;
    }

    inline virtual ~Simulation() throw()
    {
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};


YOCTO_PROGRAM_START()
{
    threading::SIMD simd(true);

    Simulation sim( 100000, simd.size );


}
YOCTO_PROGRAM_END()



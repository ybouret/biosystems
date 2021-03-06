#include "yocto/math/point3d.hpp"
#include "yocto/program.hpp"
#include "yocto/core/list.hpp"
#include "yocto/threading/scheme/simd.hpp"
#include "yocto/counted-object.hpp"
#include "yocto/ptr/arc.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/parallel/basic.hpp"
#include "yocto/random/default.hpp"
#include "yocto/code/rand32.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"

using namespace         yocto;
typedef point3d<double> Vertex;
typedef point2d<double> V2D;
typedef Random::Default RandomGenerator;

static Vertex Box(10,10,4);// initial: [-Box.x/2,Box.x/2]x[-Box.y/2,Box.y/2]x[-Box.z,0]

static inline double __ANINT( double x ) throw() { return floor(x+0.5); }


class Particle : public object
{
public:
    enum Status
    {
        Reservoir,
        Channel
    };
    Particle *next;
    Particle *prev;
    Vertex    r;
    Status    status;

    inline explicit Particle() throw() :
    next(0),
    prev(0),
    r(),
    status(Reservoir)
    {
    }

    inline virtual ~Particle() throw() {}

    inline void setInReservoir() throw()
    {
        r.x -= Box.x * __ANINT(r.x/Box.x);
        r.y -= Box.y * __ANINT(r.y/Box.y);

    CHECK_Z:
        {
            bool modified = false;
            if(r.z<-Box.z)
            {
                r.z = -(Box.z+Box.z+r.z);
                modified = true;
            }

            if(r.z>0)
            {
                r.z = -r.z;
                modified = true;
            }
            if(modified) goto CHECK_Z;
        }
    }

    bool checkReservoirToChannel(const Vertex &r_old,
                                 const Vertex &delta)  throw()
    {
        assert(r_old.z<=0);
        if(r.z>0)
        {
            // find intersection
            const double alpha = clamp<double>(0,-r_old.z/delta.z,1);
            const Vertex I(r_old.x+alpha*delta.x,r_old.y+alpha*delta.y,0);
            if(I.x*I.x+I.y*I.y<=1.0)
            {
                r      = I;
                status = Channel;
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }

    inline int move(const double     h,
                    Random::Uniform &ran,
                    unit_t          &count) throw()
    {


        switch(status)
        {
            case Reservoir: {
                const Vertex delta = h* ran.getUnit3D<Vertex>();
                const Vertex r_old = r;
                r += delta;
                if( checkReservoirToChannel(r_old,delta) )
                {
                    ++count;
                    return 1;
                }
                setInReservoir();
            } break;

            case Channel:
                break;
        }
        return 0;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Particle);
};

typedef core::list_of_cpp<Particle> Particles;

class Worker : public counted_object
{
public:
    typedef arc_ptr<Worker> Pointer;

    Particles       particles;
    Particles       inactive;
    unit_t          count; //!< last
    RandomGenerator ran;

    inline Worker() throw() : particles(), inactive(), count(0), ran( Random::SimpleBits() )
    {
    }

    inline ~Worker() throw()
    {
    }

    inline void step() throw()
    {
        count = 0;
        for(Particle *p = particles.head; p; p=p->next )
        {
            p->move(0.2, ran,count);
        }
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
    Particles       particles;  //!< original particles
    Particles       inactive;   //!< if any
    Workers         workers;    //!< workers, one per thread
    RandomGenerator ran;
    inline explicit Simulation(const size_t numParticles,const size_t numThreads) :
    particles(),
    inactive(),
    workers(numThreads),
    ran( Random::SimpleBits() )
    {
        for(size_t i=0;i<numParticles;++i)
        {
            particles.push_back( new Particle() );
        }
        std::cerr << "[SIMULATION] #particles = " << particles.size << std::endl;
        std::cerr << "[SIMULATION] #workers   = " << workers.size() << std::endl;
    }

    void setup() throw()
    {
        // gather all particles
        const size_t nw = workers.size();
        for(size_t i=nw;i>0;--i)
        {
            particles.merge_back(workers[i]->particles);
            particles.merge_back(workers[i]->inactive );
            workers[i]->count = 0; //TODO: no necessary?
        }
        particles.merge_back(inactive);
        const size_t np = particles.size;

        // put them in place
        for(Particle *p = particles.head;p;p=p->next)
        {
            p->status = Particle::Reservoir;
            p->r.x    = -Box.x/2 + Box.x * ran();
            p->r.y    = -Box.y/2 + Box.y * ran();
            p->r.z    = -Box.z * ran();

            p->setInReservoir();
        }

        // dispatch
        for(size_t i=0;i<nw;++i)
        {
            size_t offset = 0;
            size_t length = np;
            parallel::basic_split(i,nw,offset,length);
            //std::cerr << "[SIMULATION]\t-> #particles=" << length << std::endl;
            Worker &worker = *workers[i+1];
            while(length-- > 0)
            {
                worker.particles.push_back( particles.pop_back() );
            }
        }

        for(size_t i=1;i<=nw;++i)
        {
            //std::cerr << "[SIMULATION] worker[" << i << "] has " << workers[i]->particles.size << std::endl;
        }

    }

    inline void step( threading::context &ctx ) throw()
    {
        workers[ctx.indx]->step();
    }

    inline unit_t collect() const throw()
    {
        unit_t count = 0;
        for(size_t i=workers.size();i>0;--i)
        {
            count += workers[i]->count;
        }
        return count;
    }

    inline void saveXYZ( ios::ostream &fp ) const
    {
        unsigned np = 0;
        for(size_t i=workers.size();i>0;--i)
        {
            np += workers[i]->particles.size;
        }
        fp("%u\n",np);
        fp("\n");
        for(size_t i=workers.size();i>0;--i)
        {
            for(const Particle *p=workers[i]->particles.head;p;p=p->next)
            {
                fp("H %g %g %g\n",p->r.x,p->r.y,p->r.z);
            }
        }
    }

    inline virtual ~Simulation() throw()
    {
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};


YOCTO_PROGRAM_START()
{
    threading::SIMD  crew(true);
    const size_t     NP = 10000;
    Simulation        sim( NP, crew.size );
    threading::kernel k( &sim, & Simulation::step );

    sim.setup();
    {
        ios::wcstream fp("output.xyz");
        sim.saveXYZ(fp);
    }

    unit_t count = 0;
    {
        ios::wcstream fp("count.dat");
        fp("0 0\n");
    }

    for(size_t iter=1;iter<=100000;++iter)
    {
        crew.run(k);
        const unit_t local_count = sim.collect();
        if(local_count!=0)
        {
            count += local_count;
            ios::acstream fp("count.dat");
            fp("%u %g\n", unsigned(iter), double(count)/NP);
        }
        if(iter%1000==0)
        {
            ios::acstream fp("output.xyz");
            sim.saveXYZ(fp);
        }
    }



}
YOCTO_PROGRAM_END()



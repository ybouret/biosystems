#include "yocto/core/list.hpp"
#include "yocto/math/pbc.hpp"
#include "yocto/program.hpp"
#include "yocto/threading/scheme/simd.hpp"
#include "yocto/random/default.hpp"
#include "yocto/counted-object.hpp"
#include "yocto/ptr/arc.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/sort/merge.hpp"
#include "yocto/comparator.hpp"
#include "yocto/parallel/basic.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/round.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/eta.hpp"

using namespace yocto;

// types definition
typedef point3d<double> Vertex;
typedef Random::Default RandomGenerator;

// global parameters
Vertex Box(2,2,2);



class Particle : public object
{
public:
    enum Status
    {
        InReservoir
    };
    Vertex    r;
    Particle *next;
    Particle *prev;
    Status    status;
    int       kind;

    inline explicit Particle() throw() :
    r(),
    next(0), prev(0),
    status(InReservoir),
    kind(0)
    {
    }

    inline virtual ~Particle() throw()
    {
    }

    inline void inReservoir() throw()
    {
        PBC::OnXY(r,Box);

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

    static inline int compareByAddress( const Particle *lhs, const Particle *rhs, void *) throw()
    {
        const void *L = lhs;
        const void *R = rhs;
        return __compare(L,R);
    }



    inline int moveAndDetectImpact(const double     delta_lam,
                                   Random::Uniform &ran ) throw()
    {
        assert(InReservoir==status);
        int count = 0;
        Vertex delta = ran.getUnit3D<Vertex>(); delta *= delta_lam;
        Vertex r_old = r;
        assert(r_old.z<=0);
        Vertex r_new = r_old + delta;
        if(r_new.z>0)
        {
            assert(delta.z>0);
            const double    alpha = -r_old.z/delta.z;
            point2d<double> Ixy(r_old.x+alpha*delta.x,
                                r_old.y+alpha*delta.y);
            if(Ixy.norm2()< 1)
            {
                ++count;
                kind = 1;
            }
        }

        r = r_new;
        inReservoir();

        return count;
    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Particle);
};

typedef core::list_of_cpp<Particle> Particles;

class Worker : public counted_object
{
public:
    typedef arc_ptr<Worker> Pointer;
    unit_t          count;
    Particles       running;
    Particles       waiting;
    RandomGenerator ran;

    inline explicit Worker() :
    count(0),
    running(),
    waiting(),
    ran( Random::SimpleBits() )
    {
    }

    inline virtual ~Worker() throw()
    {
    }

    void detectImpact(const double     delta_lam,
                      Random::Uniform &ran) throw()
    {
        count = 0;
        Particles tmp;
        while(running.size>0)
        {
            Particle *p = running.pop_back();
            const int ans = p->moveAndDetectImpact(delta_lam,ran);
            if(ans!=0)
            {
                count += ans;
                waiting.push_back(p);
            }
            else
            {
                tmp.push_front(p);
            }
        }
        tmp.swap_with(running);
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Worker);
};

typedef vector<Worker::Pointer> Workers;

class Simulation
{
public:
    const double    delta_lam;
    const double    delta_tau;
    Workers         workers;
    Particles       running;
    Particles       waiting;
    RandomGenerator ran;

    inline explicit Simulation(const size_t np,
                               const size_t nw,
                               const double delta_lambda) :
    delta_lam(delta_lambda),
    delta_tau(delta_lam*delta_lam),
    workers(nw,as_capacity),
    running(),
    waiting(),
    ran( Random::SimpleBits() )
    {
        for(size_t i=np;i>0;--i) { running.push_back( new Particle() ); }
        for(size_t i=nw;i>0;--i) { const Worker::Pointer pW( new Worker() ); workers.push_back(pW); }
    }

    inline void setup() throw()
    {
        // gather particles
        const size_t nw = workers.size();
        running.merge_back(waiting);
        for(size_t i=nw;i>0;--i)
        {
            running.merge_back(workers[i]->running);
            running.merge_back(workers[i]->waiting);
        }
        core::merging<Particle>::sort(running, Particle::compareByAddress, NULL);

        // initialize them
        for(Particle *p=running.head;p;p=p->next)
        {
            p->r.x = -Box.x/2 + Box.x * ran();
            p->r.y = -Box.y/2 + Box.y * ran();
            p->r.z = -Box.z * ran();
            p->status = Particle::InReservoir;
            p->kind   = 0;
            p->inReservoir();
        }

        // dispatch particles
        const size_t np = running.size;
        for(size_t i=0;i<nw;++i)
        {
            size_t offset = 0;
            size_t length = np;
            parallel::basic_split(i,nw, offset, length);
            Worker &w = *workers[i+1];
            while(length-->0)
            {
                w.running.push_back(running.pop_front());
            }
            //std::cerr << "#running[" << i+1 << "]=" << w.running.size << std::endl;
        }
    }

    inline virtual ~Simulation() throw()
    {
    }

    void saveXYZ( ios::ostream &fp ) const
    {
        unsigned np = 0;
        const size_t nw = workers.size();
        for(size_t i=nw;i>0;--i)
        {
            np += workers[i]->running.size;
        }
        fp("%u\n", np);
        fp("\n");
        for(size_t i=nw;i>0;--i)
        {
            for(const Particle *p=workers[i]->running.head;p;p=p->next)
            {
                const char *id = "H";
                if(p->kind) id = "Li";
                fp("%s %g %g %g\n",id,p->r.x,p->r.y,p->r.z);
            }
        }
    }

    // parallel kernel
    inline void run( threading::context &ctx ) throw()
    {
        workers[ctx.indx]->detectImpact(delta_lam,ran);
    }

    inline unit_t getCount() const throw()
    {
        unit_t count = 0;
        for(size_t i=workers.size();i>0;--i)
        {
            count += workers[i]->count;
        }
        return count;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};
#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    double lambda = 0.1;
    if(argc>1)
    {
        lambda = strconv::to_double(argv[1],"lambda");
    }
    if(lambda<=0) throw exception("lambda<=0");

    size_t NP = 10000;
    threading::SIMD   engine(true);
    Simulation        sim(NP,engine.size,lambda);
    threading::kernel kernelRun( &sim, &Simulation::run);


    const double dtau     = sim.delta_tau;
    double       tauMax   = 1;
    double       dtauSave = 0.01;

    const size_t every = math::simulation_save_every_(dtau, dtauSave);
    const size_t iter  = math::simulation_iter(tauMax,dtau,every);

    std::cerr << "tauMax = " << tauMax << std::endl;
    std::cerr << "dtau   = " << dtau   << std::endl;
    std::cerr << "iter   = " << iter   << std::endl;
    std::cerr << "every  = " << every  << std::endl;

    eta ETA;
    vector<unit_t> F(iter,0);
    vector<double> Q(iter,0);
    const string   dataName = vformat("q%g.dat",dtau);

    for(size_t cycle=1;cycle<=10;++cycle)
    {
        ETA.reset();
        sim.setup();



        double all      = 0;
        double tau      = 0;

        for(size_t i=1;i<=iter;++i)
        {
            tau += dtau;
            engine.run(kernelRun);
            const unit_t count = sim.getCount();
            F[i] += count;
            all  += double(F[i])/cycle;
            Q[i]  = all/NP;
            if( (i%every) == 0)
            {

                {
                    ios::wcstream fp(dataName);
                    for(size_t j=1;j<=iter;++j)
                    {
                        fp("%g %g\n",j*dtau,Q[j]);
                    }
                }
                ETA.progress(double(i)/iter);
            }
        }
        ETA.progress_flush();
        std::cerr << std::endl;
    }
}
YOCTO_PROGRAM_END()



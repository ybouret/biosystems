#include "yocto/program.hpp"
#include "yocto/threading/crew.hpp"
#include "yocto/math/point3d.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;


typedef float         Real;
typedef point3d<Real> Point;

static Real L=0, xmin=0, xmax=0, ymin=0, ymax=0, zmin=0;
static Real H = 0, Radius=0;

static void SetupGeometry()
{
    xmin = -L/2;
    xmax =  L/2;
    ymin = -L/2;
    ymax =  L/2;
    zmin = -L;
}


#define IN_BOX 0x0001
#define IN_CYL 0x0002

#define LI7 7
#define LI6 6

class Particle
{
public:
    Point     r;
    Real      step_length;
    int       flags;
    int       type;

    inline Particle() throw() : r(), step_length(0), flags(0), type(0)
    {
    }

    inline ~Particle() throw() {}

    void move()
    {
        
    }
    
private:
    YOCTO_DISABLE_ASSIGN(Particle);
};


class Simulation
{
public:
    typedef threading::context Context;
    typedef Context::range     Range;
    virtual ~Simulation() throw() {}

    vector<Particle>                      particles;
    threading::crew                       engine;
    uniform_generator<Real,rand32_kiss>   ran;
    threading::kernel                     kStep;

    explicit Simulation(const size_t n) :
    particles(n),
    engine(true),
    ran(),
    kStep(this, & Simulation::StepCall )
    {

        for(size_t i=0;i<engine.size;++i)
        {
            threading::context &ctx = engine[i];
            ctx.create_range(n);
        }
    }


    void initialize()
    {
        for(size_t i=particles.size();i>0;--i)
        {
            Particle &p = particles[i];
            p.r.x   = clamp<Real>(xmin,xmin + L * ran(),xmax);
            p.r.y   = clamp<Real>(ymin,ymin + L * ran(),ymax);
            p.r.z   = clamp<Real>(zmin,zmin + L * ran(),0);
            p.flags = IN_BOX;
            p.type  = LI7;
        }
    }

    void append_to( const string &filename, const Real tau = 0) const
    {
        ios::acstream fp(filename);
        const unsigned n = particles.size();
        fp("%u\n",n);
        fp("t=%.15g\n", tau);
        for(size_t i=1;i<=n;++i)
        {
            const Particle &p = particles[i];
            switch(p.type)
            {
                case LI7: fp("LI7"); break;
                case LI6: fp("LI6"); break;
                default:  fp("H");   break;
            }
            fp(" %.7g %.7g %.7g\n", p.r.x, p.r.y, p.r.z);
        }
    }

    void StepCall( threading::context &context  ) throw()
    {
        const Range &range = context.data.as<Range>();
        if(false)
        {
            //YOCTO_LOCK(context.access);
            //std::cerr << range << std::endl;
        }
        for(size_t i=range.offset,count=range.length;count>0;--count,++i)
        {
            Particle &p = particles[i];
            p.move();
        }
    }


    void step()
    {
        engine(kStep);
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};

YOCTO_PROGRAM_START()
{
    L = 1;
    H = 1;
    Radius = 0;
    SetupGeometry();

    Simulation sim(100);
    sim.initialize();
    ios::ocstream::overwrite("sim.xyz");
    sim.append_to("sim.xyz");
    sim.step();
    
}
YOCTO_PROGRAM_END()

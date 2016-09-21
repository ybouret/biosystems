#include "yocto/program.hpp"
#include "yocto/threading/vpu.hpp"
#include "yocto/math/point3d.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/code/utils.hpp"

using namespace yocto;
using namespace math;


typedef float         Real;
typedef point3d<Real> Point;

static Real L=0, xmin=0, xmax=0, ymin=0, ymax=0, zmin=0;
static Real H=0, Radius=0;

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

class Particle
{
public:
    Point     r;
    Real      step_length;
    int       flags;
    inline Particle() throw() : r(), step_length(0), flags(0)
    {
    }

    inline ~Particle() throw() {}

private:
    YOCTO_DISABLE_ASSIGN(Particle);
};

class Processor
{
public:
    inline Processor() throw()
    {
    }

    inline ~Processor() throw()
    {
    }

    inline
    void run( threading::context &ctx, array<Particle> &particles, void *) throw()
    {
        size_t offset = 1;
        size_t length = particles.size();
        ctx.split(offset, length);
        
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Processor);
};

class Simulation
{
public:
    virtual ~Simulation() throw() {}

    vector<Particle>                      particles;
    uniform_generator<Real,rand32_kiss>   ran;
    threading::processing_unit<Processor> cpu;

    explicit Simulation(const size_t n) :
    particles(n),
    ran(),
    cpu( new threading::crew(true) )
    {
        // prepare processors
        for(size_t i=0;i<cpu.cores;++i)
        {
            cpu.push_back();// create Processor()
        }
        cpu.compile<Particle>();


        for(size_t i=particles.size();i>0;--i)
        {
            Particle &p = particles[i];
            p.r.x   = clamp<Real>(xmin,xmin + L * ran(),xmax);
            p.r.y   = clamp<Real>(ymin,ymin + L * ran(),ymax);
            p.r.z   = clamp<Real>(zmin,zmin + L * ran(),0);
            p.flags = IN_BOX;
        }
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

    Simulation sim(10);


}
YOCTO_PROGRAM_END()

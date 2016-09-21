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
static Real H = 0;

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

#define Li7 7
#define Li6 6

typedef uniform_generator<Real,rand32_kiss> Rand32;

class RunTime
{
public:
    size_t          offset;
    size_t          length;
    mutable Rand32  ran;

    RunTime() throw() : offset(0), length(0), ran()
    {
        ran.wseed();
    }

    ~RunTime() throw()
    {
    }

    void prepare(const threading::context &ctx, const size_t n)
    {
        offset = 1;
        length = n;
        ctx.split<size_t>(offset,length);
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(RunTime);
};

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

    void move(Rand32 &ran)
    {
        // create displacement
        Point dr = Point::on_unit_sphere(ran);
        dr *= step_length;
        Point src = r;
        Point tgt = r+dr;

        if(flags==IN_BOX)
        {
        CHECK_BOX:
            if(tgt.z<zmin)
            {
                const Real dz = zmin - tgt.z;
                tgt.z = zmin + dz;
                goto CHECK_BOX;
            }

            if(tgt.x>xmax)
            {
                const Real dx = tgt.x - xmax;
                tgt.x = xmax - dx;
                goto CHECK_BOX;
            }

            if(tgt.x<xmin)
            {
                const Real dx = xmin - tgt.x;
                tgt.x = xmin + dx;
                goto CHECK_BOX;
            }

            if(tgt.y>ymax)
            {
                const Real dy = tgt.y - ymax;
                tgt.y = ymax - dy;
                goto CHECK_BOX;
            }

            if(tgt.y<ymin)
            {
                const Real dy = ymin - tgt.y;
                tgt.y = ymin + dy;
                goto CHECK_BOX;
            }


        }

        // assign new point
        r = tgt;

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
    threading::kernel                     kStep;

    explicit Simulation(const size_t n) :
    particles(n),
    engine(true),
    kStep(this, & Simulation::StepCall )
    {

        for(size_t i=0;i<engine.size;++i)
        {
            threading::context &ctx = engine[i];
            ctx.data.make<RunTime>().prepare(ctx,n);
        }
    }


    void initialize()
    {
        Rand32 &ran = engine[0].data.as<RunTime>().ran;
        for(size_t i=particles.size();i>0;--i)
        {
            Particle &p = particles[i];
            p.r.x   = clamp<Real>(xmin,xmin + L * ran(),xmax);
            p.r.y   = clamp<Real>(ymin,ymin + L * ran(),ymax);
            p.r.z   = clamp<Real>(zmin,zmin + L * ran(),0);
            p.flags = IN_BOX;
            p.type  = Li7;
        }

        for(size_t i=particles.size();i>0;--i)
        {
            Particle &p = particles[i];
            p.step_length = 0;
            switch(p.type)
            {
                case Li6: break;
                case Li7: p.step_length = 0.5; break;
                default:  break;
            }
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
                case Li7: fp("Li7"); break;
                case Li6: fp("Li6"); break;
                default:  fp("H");   break;
            }
            fp(" %.7g %.7g %.7g\n", p.r.x, p.r.y, p.r.z);
        }
    }

    void StepCall( threading::context &context  ) throw()
    {
        const RunTime &range = context.data.as<RunTime>();
        if(false)
        {
            //YOCTO_LOCK(context.access);
            //std::cerr << range << std::endl;
        }
        Rand32 &ran = range.ran;
        for(size_t i=range.offset,count=range.length;count>0;--count,++i)
        {
            Particle &p = particles[i];
            p.move(ran);
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
    SetupGeometry();

    Simulation sim(1000);
    sim.initialize();
    ios::ocstream::overwrite("sim.xyz");
    sim.append_to("sim.xyz");

    for(size_t i=0;i<200;++i)
    {
        sim.step();
        sim.append_to("sim.xyz");
    }
    
}
YOCTO_PROGRAM_END()

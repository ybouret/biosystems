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
    mutable Real    tau;
    mutable vector<Real>    tau6;
    mutable vector<Real>    tau7;

    RunTime() throw() : offset(0), length(0), ran(), tau(0)
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
        tau6.free();
        tau7.free();
        tau6.ensure(length);
        tau7.ensure(length);
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

    inline Particle() throw() : r(), step_length(0),  flags(0), type(0)
    {
    }

    inline ~Particle() throw() {}

    void move(const Real      tau,
              Rand32         &ran,
              sequence<Real> &tau6,
              sequence<Real> &tau7)
    {


        if(flags<=0) return;

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
                src=tgt; goto CHECK_BOX;
            }

            if(tgt.x>xmax)
            {
                const Real dx = tgt.x - xmax;
                tgt.x = xmax - dx;
                src=tgt; goto CHECK_BOX;
            }

            if(tgt.x<xmin)
            {
                const Real dx = xmin - tgt.x;
                tgt.x = xmin + dx;
                src=tgt; goto CHECK_BOX;
            }

            if(tgt.y>ymax)
            {
                const Real dy = tgt.y - ymax;
                tgt.y = ymax - dy;
                src=tgt; goto CHECK_BOX;
            }

            if(tgt.y<ymin)
            {
                const Real dy = ymin - tgt.y;
                tgt.y = ymin + dy;
                src=tgt; goto CHECK_BOX;
            }

            if(tgt.z>0)
            {
                // do we enter the cylinder
                dr = tgt-src;
                const Real  lam = clamp<Real>(0,-src.z/dr.z,1);
                const Point I   = src + lam * dr;
                const Real  r   = Hypotenuse(I.x,I.y);
                if(r<1)
                {
                    flags = 0;
                    switch(type)
                    {
                        case Li6: tau6.push_back(tau); break;
                        case Li7: tau7.push_back(tau); break;
                        default:
                            break;
                    }
                }
                tgt.z=-tgt.z;
                src=tgt; goto CHECK_BOX;
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
    size_t                                active;

    explicit Simulation(const size_t n) :
    particles(n),
    engine(true),
    kStep(this, & Simulation::StepCall ),
    active(0)
    {

        for(size_t i=0;i<engine.size;++i)
        {
            threading::context &ctx = engine[i];
            ctx.make<RunTime>().prepare(ctx,n);
        }
    }


    void initialize()
    {
        Rand32 &ran = engine[0].as<RunTime>().ran;
        active = particles.size();
        for(size_t i=active;i>0;--i)
        {
            Particle &p = particles[i];
            p.r.x   = clamp<Real>(xmin,xmin + L * ran(),xmax);
            p.r.y   = clamp<Real>(ymin,ymin + L * ran(),ymax);
            p.r.z   = clamp<Real>(zmin,zmin + L * ran(),0);
            p.flags = IN_BOX;
            p.type  = Li7;
            if(ran.get<double>()<0.5) p.type = Li6;
        }

        const Real step_length7 = 0.1;
        const Real step_length6 = step_length7 * sqrt(7.0/6.0);
        for(size_t i=particles.size();i>0;--i)
        {
            Particle &p = particles[i];
            p.step_length = 0;
            switch(p.type)
            {
                case Li6: p.step_length = step_length6; break;
                case Li7: p.step_length = step_length7; break;
                default:  break;
            }
        }

    }

    void append_to( const string &filename, const Real tau = 0) const
    {
        ios::acstream fp(filename);
        fp("%u\n",unsigned(active));
        fp("t=%.15g\n", tau);
        for(size_t i=1;i<=particles.size();++i)
        {
            const Particle &p = particles[i];
            if(p.flags<=0) continue;

            switch(p.type)
            {
                case Li7: fp("Li"); break;
                case Li6: fp("He"); break;
                default:  fp("H");   break;
            }
            fp(" %.7g %.7g %.7g\n", p.r.x, p.r.y, p.r.z);
        }
    }

    void StepCall( threading::context &context) throw()
    {
        const RunTime &rt = context.as<RunTime>();
        if(false)
        {
            //YOCTO_LOCK(context.access);
            //std::cerr << range << std::endl;
        }
        Rand32 &ran = rt.ran;
        rt.tau6.free();
        rt.tau7.free();
        for(size_t i=rt.offset,count=rt.length;count>0;--count,++i)
        {
            Particle &p = particles[i];
            p.move(rt.tau,ran,rt.tau6,rt.tau7);
        }
    }


    void step(const Real tau)
    {
        for(size_t i=0;i<engine.size;++i)
        {
            const RunTime &rt = engine[i].as<RunTime>();
            rt.tau = tau;
        }

        engine(kStep);

        for(size_t i=0;i<engine.size;++i)
        {
            const RunTime &rt = engine[i].as<RunTime>();
            const array<Real> &tau6 = rt.tau6;
            const array<Real> &tau7 = rt.tau7;
            assert(tau6.size()+tau7.size()<=active);
            active -= tau6.size();
            active -= tau7.size();
        }
        std::cerr << "#active=" << active << std::endl;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};

YOCTO_PROGRAM_START()
{
    L = 4;
    H = 1;
    SetupGeometry();

    Simulation sim(1000);
    sim.initialize();
    ios::ocstream::overwrite("sim.xyz");
    sim.append_to("sim.xyz");

    const Real dtau = 1;
    for(size_t i=1;i<=20000;++i)
    {
        sim.step(dtau*i);
        if(0==(i%100))
        {
            sim.append_to("sim.xyz");
        }
    }
    
}
YOCTO_PROGRAM_END()

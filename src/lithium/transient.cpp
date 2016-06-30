#include "yocto/program.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/math/ode/explicit/driver-ck.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;

#define NVAR      9

#define _EE        1

#define _Li6_out   2
#define _Li6EE_out 3
#define _Li6EE_in  4
#define _Li6_in    5

#define _Li7_out   6
#define _Li7EE_out 7
#define _Li7EE_in  8
#define _Li7_in    9

typedef numeric<double>::function    Function;
typedef ode::Field<double>::Equation Equation;
typedef ode::Field<double>::Callback Callback;

class LiSystem
{
public:
    lua_State      *L;
    Equation       eqs;
    vector<double> C;
    ode::driverCK<double>::type odeint;
    const double   Li6_ini;
    const double   Li7_ini;
    const double   E_ini;
    const double   ka6;
    const double   kd6;
    const double   kf6;
    const double   kr6;
    const double   kt6;

    const double   ka7;
    const double   kd7;
    const double   kf7;
    const double   kr7;
    const double   kt7;

    double         ratio0;

#define LiLoad(NAME) NAME( Lua::Config::Get<lua_Number>(L, #NAME ) )

    explicit LiSystem( Lua::State &VM) :
    L( VM() ),
    eqs( this, & LiSystem::equations ),
    C(NVAR),
    odeint(1e-5),
    Li6_ini(Lua::Config::Get<lua_Number>(L, "Li6" )),
    Li7_ini(Lua::Config::Get<lua_Number>(L, "Li7" )),
    E_ini  (Lua::Config::Get<lua_Number>(L, "E"   )),
    LiLoad(ka6),
    LiLoad(kd6),
    LiLoad(kf6),
    LiLoad(kr6),
    LiLoad(kt6),
    LiLoad(ka7),
    LiLoad(kd7),
    LiLoad(kf7),
    LiLoad(kr7),
    LiLoad(kt7),
    ratio0(0)
    {
        odeint.start(C.size());
    }

    virtual ~LiSystem() throw()
    {
    }

    void equations( array<double> &dYdt, const double, const array<double> &Y )
    {
        tao::ld(dYdt,0);
        const double  EE      = Y[_EE];
        const double  Li6_out = Y[_Li6_out];
        const double  Li7_out = Y[_Li7_out];

        const double  va6  = ka6 * EE * Li6_out;
        const double  va7  = ka7 * EE * Li7_out;

        const double  Li6EE_out = Y[_Li6EE_out];
        const double  Li7EE_out = Y[_Li7EE_out];

        const double  vd6       = kd6 * Li6EE_out;
        const double  vd7       = kd7 * Li7EE_out;

        const double  vf6       = kf6 * Li6EE_out;
        const double  vf7       = kf7 * Li7EE_out;

        const double  Li6EE_in = Y[_Li6EE_in];
        const double  Li7EE_in = Y[_Li7EE_in];

        const double  vr6       = kr6 * Li6EE_in;
        const double  vr7       = kr7 * Li7EE_in;

        const double  vt6       = kt6 * Li6EE_in;
        const double  vt7       = kt7 * Li7EE_in;


        dYdt[_EE]        = (vt6+vt7) + (vd6+vd7) - (va6+va7);

        dYdt[_Li6EE_out] = (va6+vr6) - (vd6+vf6);
        dYdt[_Li7EE_out] = (va7+vr7) - (vd7+vf7);

        dYdt[_Li6EE_in] = vf6 - (vr6+vt6);
        dYdt[_Li7EE_in] = vf7 - (vr7+vt7);

        dYdt[_Li6_in]   = vt6;
        dYdt[_Li7_in]   = vt7;


    }


    inline void initialize()
    {
        tao::ld(C,0);
        C[_Li6_out] = Li6_ini;
        C[_Li7_out] = Li7_ini;
        C[_EE]      = E_ini;
        ratio0      = 0;

        if(Li7_ini>0)
        {
            ratio0 = Li6_ini/Li7_ini;
        }

    }
    
    inline void step(const double t1,const double t2)
    {
        double ctrl = (t2-t1)/100.0;
        odeint(eqs,C,t1,t2,ctrl,NULL);
    }

    inline void prolog( ios::ostream &fp ) const
    {
        fp("#t E(2) Li6_out(3) Li6EE_out(4) Li6EE_in(5) Li6_in(6) Li7_out(7) Li7EE_out(8) Li7EE_in(9) Li7_in(10) delta\n");
    }

    inline void save(const double t, ios::ostream &fp) const
    {
        fp("%.15g", t);
        for(size_t i=1;i<=C.size();++i)
        {
            fp(" %.15g", C[i]);
        }
        const double Li7 = C[_Li7_in];
        if(Li7>0)
        {
            const double ratio1 = C[_Li6_in]/Li7;
            double       delta7 = 0;
            if(ratio0>0)
            {
                delta7 = (ratio1-ratio0)/ratio0;
            }
            fp(" %g", delta7);
        }
        fp("\n");
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(LiSystem);
};


YOCTO_PROGRAM_START()
{
    Lua::State VM;
    lua_State *L = VM();

    for(int i=1;i<argc;++i)
    {
        Lua::Config::DoFile(L,argv[i]);
    }

    const double Tmax = Lua::Config::Get<lua_Number>(L,"Tmax");

    LiSystem sys(VM);

    ios::wcstream fp("sim.dat");
    sys.prolog(fp);
    sys.initialize();
    double       t1 = 0;
    const double dt = 0.01;
    size_t it = 0;
    sys.save(t1,fp);
    while(t1<Tmax)
    {
        const double t2 = ++it * dt;
        sys.step(t1, t2);
        t1 = t2;
        sys.save(t1,fp);
    }



}
YOCTO_PROGRAM_END()

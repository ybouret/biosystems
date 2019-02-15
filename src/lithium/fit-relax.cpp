
#include "y/program.hpp"
#include "y/lang/pattern/matching.hpp"
#include "y/math/io/data-set.hpp"
#include "y/math/fit/ls.hpp"
#include "y/ptr/shared.hpp"
#include "y/sort/unique.hpp"
#include "y/math/fit/vectors.hpp"
#include "y/math/stat/metrics.hpp"
#include "y/math/fcn/zfind.hpp"
#include "y/sort/sorted-sum.hpp"
#include "y/sequence/maintain.hpp"

using namespace upsylon;
using namespace math;

typedef Fit::Sample<double>     Sample;
typedef Fit::Samples<double>    Samples;
typedef vector<double>          Vector;
typedef array<double>           Array;
typedef Fit::Variables          Variables;

static const double d7Out = 14.57;
static const double sigma = 1.0029;

class Leak
{
public:
    inline Leak()
    {
    }

    inline ~Leak() throw()
    {
    }

    static inline
    double delta_of( const double r )
    {
        return 1000.0 * ( (1.0+d7Out/1000.0) * r - 1.0 );
    }

    static inline
    double ratio_of( const double d )
    {
        return (1.0+d/1000.0)/(1.0+d7Out/1000.0);
    }

    double Compute(double           t,
                   const Array     &aorg,
                   const Variables &vars)
    {
        const double B_star = vars( aorg, "B" );
        const double r_star = vars( aorg, "r" );
        const double k7     = vars( aorg, "k7" );
        const double t0     = vars( aorg, "t0" );

        if(t<=t0)
        {
            return delta_of( r_star );
        }
        else
        {
            const double rB  = r_star * B_star;
            const double rp1 = r_star+1.0;
            const double tau = k7*(t-t0);
            const double b7  = rB     + (rp1 - rB    )*(1.0-exp(      -tau) );
            const double b6  = B_star + (rp1 - B_star)*(1.0-exp(-sigma*tau) );
            const double rr  = b7/b6;
            return delta_of( rr );
        }

    }

    static inline bool Keep( const double t )
    {
        return t>=60;
    }

    static inline void Save( const Sample &s )
    {
        ios::ocstream fp("fit-relax.dat");
        for(size_t i=1;i<=s.X.size();++i)
        {

            fp("%.15g %.15g %.15g\n", s.X[i], s.Y[i], s.Yf[i] );
        }
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Leak);
};


Y_PROGRAM_START()
{

    if(argc<=1)
    {
        return 0;
    }

    const string file_name = argv[1];
    Vector t;
    Vector d7;
    {
        ios::icstream fp(file_name);
        data_set<double> ds;
        ds.use(1,t);
        ds.use(2,d7);
        ds.load(fp);
    }

    // cut

    std::cerr << "#loaded=" << t.size() << std::endl;
    {
        vector<size_t> good;
        maintain::build_indices(good, t, Leak::Keep,true);
        std::cerr << "#good=" << good.size() << std::endl;

        maintain::strip(t,good);
        maintain::strip(d7,good);

    }

    const size_t np = t.size();
    Vector d7fit(np,0);

    Sample sample(t,d7,d7fit);


    Variables &vars = sample.variables;
    vars << "B" << "r" << "k7" << "t0";

    const size_t nv = vars.size();

    Vector aorg(nv,0);
    Vector aerr(nv,0);
    vector<bool> used(nv,false);

    vars(aorg,"k7") = 7e-5;
    vars(aorg,"t0") = t[1];
    vars(aorg,"r")  = Leak::ratio_of(d7[1]);
    vars(aorg,"B")  = 0.5;

    Leak leak;
    Fit::LeastSquares<double>::Function F( &leak, & Leak::Compute );
    Fit::LeastSquares<double> LS;

    sample.computeD2(F,aorg);

    Leak::Save(sample);

    vars.on(used,"k7");

    int level=0;

    if(!LS.fit(sample, F, aorg, aerr, used) )
    {
        ++level;
        std::cerr << "couldn't fit level-" << level << std::endl;
    }

    Leak::Save(sample);
    vars.display(std::cerr, aorg,aerr);
    std::cerr << std::endl;

    vars.on(used,"B");
    if(!LS.fit(sample, F, aorg, aerr, used) )
    {
        ++level;
        std::cerr << "couldn't fit level-" << level << std::endl;
    }
    Leak::Save(sample);
    vars.display(std::cerr, aorg,aerr);
    std::cerr << std::endl;

    vars.on(used,"r");
    if(!LS.fit(sample, F, aorg, aerr, used) )
    {
        ++level;
        std::cerr << "couldn't fit level-" << level << std::endl;
    }
    Leak::Save(sample);
    vars.display(std::cerr, aorg,aerr);
    std::cerr << std::endl;

}
Y_PROGRAM_END()


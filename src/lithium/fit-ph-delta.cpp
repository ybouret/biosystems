
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


class delta_pH
{
public:
    inline  delta_pH() throw() {}
    inline ~delta_pH() throw() {}

    inline double compute( double Li, const Array &aorg, const Variables &vars )
    {
        const double L = vars(aorg,"L");
        const double p = vars(aorg,"p");
        const double A = vars(aorg,"A");
        const double xx = pow((Li/L),p);
        return A * xx / (1.0+xx);
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(delta_pH);
};

static inline
void save( const Sample &sample )
{
    ios::ocstream fp("pH_delta.fit");
    for(size_t i=1;i<=sample.count();++i)
    {
        fp("%.15g %.15g %.15g\n", sample.X[i], sample.Y[i], sample.Yf[i]);
    }

}

Y_PROGRAM_START()
{
    if(argc<=1)
    {
        throw exception("need a pH_delta.dat");
    }

    // load data
    Vector       Li;
    Vector       dpH;
    const string file_name = argv[1];
    const size_t n = data_set<double>::loadXY(file_name, 1, Li, 2, dpH);
    Vector       Yf(n,0);

    // create sample
    Sample   sample(Li,dpH,Yf);

    // create variables
    Variables &vars = sample.variables;

    vars << "L" << "A" << "p";

    const size_t nv = vars.size();
    Vector       aorg(nv,0);
    Vector       aerr(nv,0);
    vector<bool> used(nv,false);

    vars(aorg,"A") = dpH[n];
    vars(aorg,"L") = Li[1+n/2];
    vars(aorg,"p") = 1;

    vars.on(used,"A:L");

    // create fit function
    delta_pH dd;
    Fit::Type<double>::Function F( &dd, & delta_pH::compute );

    // create leasy square
    int level = 0;
    Fit::LeastSquares<double>   ls;

    {
        std::cerr << std::endl;
        ++level;
        if( !ls.fit(sample, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit level-%d", level);
        }

        vars.display(std::cerr, aorg, aerr);
        save(sample);
    }


    if(false)
    {
        vars.on(used,"p");
        vars(aorg,"p") =  0.5;
        std::cerr << std::endl;
        ++level;
        if( !ls.fit(sample, F, aorg, aerr, used) )
        {
            throw exception("couldn't fit level-%d", level);
        }

        vars.display(std::cerr, aorg, aerr);
        save(sample);
    }



}
Y_PROGRAM_END()


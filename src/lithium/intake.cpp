#include "y/program.hpp"
#include "y/math/ode/explicit/driver-ck.hpp"

using namespace upsylon;
using namespace math;

typedef array<double> Array;

class Intake
{
public:
    static const size_t NVAR = 3;
    static const size_t I_ALPHA = 1;
    static const size_t I_LI6   = 2;
    static const size_t I_LI7   = 3;

    double k6;
    double k7;

    void Compute( Array &dYdt, double t, Array &Y)
    {
        const double alpha = Y[I_ALPHA];
        const double Li6   = Y[I_LI6];
        const double Li7   = Y[I_LI7];
        const double h     = get_h(t);

        dYdt[I_ALPHA] = 0;
        dYdt[I_LI6]   = 0;
        dYdt[I_LI7]   = 0;
    }

    double get_h(double)
    {
        return pow(1.0,-5.9);
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Intake);
};

Y_PROGRAM_START()
{
    ODE::DriverCK<double>::Type       odeint;

    odeint.start( Intake::NVAR );

}
Y_PROGRAM_END()



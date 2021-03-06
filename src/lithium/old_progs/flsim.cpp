#include "y/program.hpp"
#include "sim.h"
#include "y/ptr/auto.hpp"
#include "y/alea.hpp"
#include "fcnli.hpp"

static const size_t NP = 1000;



void update()
{
    Ca->xaxis.set_min( logMin->value() );
    Ca->xaxis.set_range( logMin->value(), logMax->value() );

    FLTK::Curve &beta7 = Ca->curves["beta7"];
    FLTK::Curve &beta6 = Ca->curves["beta6"];
    FLTK::Curve &rho   = Ca->curves2["rho"];

    beta7.free();
    beta6.free();
    rho.free();
    beta7.ensure(NP+1);
    beta6.ensure(NP+1);
    rho.ensure(NP+1);

    beta6.color = FL_RED;
    beta7.color = FL_GREEN;

    const double sigma = inputSigma->value();
    const double A     = pow(10.0,inputLog10A->value());
    const double c2    = inputC2->value();
    const double s2    = 1-c2;
    const double Omega = inputOmega->value();
    const double kappa = inputKappa->value();

    for(size_t i=0;i<=NP;++i)
    {
        const double lt  = Ca->xaxis.vmin + i*(Ca->xaxis.vmax-Ca->xaxis.vmin)/NP;
        const double tau = exp(lt);
        const double b7 = Lithium::Grow(tau) * (1+A*c2) + A*s2 * Lithium::Bump(Omega,tau);
        const double b6 = Lithium::Grow(tau*sigma) * (1+kappa*A*c2) + kappa *A *s2 * Lithium::Bump(Omega/sigma,tau*sigma);

        beta6.push_back( FLTK::Point(lt,b6) );
        beta7.push_back( FLTK::Point(lt,b7) );
        rho.push_back( FLTK::Point(lt,b7/b6) );
    }

    double rmin = rho.front().y;
    double rmax = rmin;

    double bmin = beta6.front().y;
    double bmax = bmin;

    for(size_t i=2;i<=rho.size();++i)
    {
        const double r = rho[i].y;
        rmin = min_of(r,rmin);
        rmax = max_of(r,rmax);

        bmin = min_of(bmin,min_of( beta6[i].y, beta7[i].y ));
        bmax = max_of(bmax,max_of( beta6[i].y, beta7[i].y ));

    }

    Ca->y2axis.set_range(rmin,rmax);
    Ca->yaxis.set_range(bmin, bmax);

    Ca->redraw();

}


Y_PROGRAM_START()
{
    alea_init();

    auto_ptr<Fl_Window> win( MakeSimWindow() );

    // initialize
    update();
    win->show(argc,argv);



    (void) Fl::run();
}
Y_PROGRAM_END()


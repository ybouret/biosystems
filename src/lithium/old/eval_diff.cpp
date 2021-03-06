#include "yocto/code/rand32.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/program.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/point3d.hpp"

using namespace yocto;
using namespace math;

namespace
{
    static inline
    double compute_variance(const array<double> &u, double &average)
    {
        const size_t n = u.size();
        average = 0;
        for(size_t i=n;i>0;--i)
        {
            average += u[i];
        }
        average/=n;
        double var = 0;
        for(size_t i=n;i>0;--i)
        {
            var += math::Square( u[i] - average );
        }
        return var/n;
    }
}

YOCTO_PROGRAM_START()
{
    size_t       N = 10000;
    rand32_kiss ran;
    ran.wseed();

    const double Ddt   = 1.0;    // scaling
    const double gvar  = 2*Ddt;  // gaussian variance
    const double fac   = sqrt(gvar);

    // 1D
    {
        std::cerr << std::endl;
        vector<double> X(N);
        for(size_t i=1;i<=N;++i)
        {
            X[i] = Ddt * ran.sym1<double>();
        }
        double ave=0;
        double var=compute_variance(X,ave);

        for(size_t i=1;;)
        {
            double a,b;
            ran.gaussian(a,b);
            a *= fac;
            b *= fac;
            if(i<=N)
            {
                X[i++] = a;
            }
            else
            {
                break;
            }

            if(i<=N)
            {
                X[i++] = b;
            }
            else
            {
                break;
            }
        }
        double aveG = 0;
        double varG = compute_variance(X,aveG);

        std::cerr << "1D uniform : average=" << ave << ", variance=" << var << std::endl;
        std::cerr << "1D gaussian: average=" << aveG << ", variance=" << varG << std::endl;
        std::cerr << "1D decrease=" << floor(varG/var+0.5) << std::endl;

    }

    // 2D
    typedef point2d<double> P2D;
    {
        std::cerr << std::endl;
        vector<double> X(N),Y(N);
        for(size_t i=1;i<=N;++i)
        {
            const P2D p = Ddt * P2D::on_unit_circle(ran);
            X[i]  = p.x;
            Y[i]  = p.y;
        }
        P2D ave,var;
        var.x = compute_variance(X,ave.x);
        var.y = compute_variance(Y,ave.y);
        double r_var = 0;
        for(size_t i=1;i<=N;++i)
        {
            r_var += Square(X[i]-ave.x) + Square(Y[i]-ave.y);
        }
        r_var/=N;


        for(size_t i=1;i<=N;++i)
        {
            double a,b;
            ran.gaussian(a,b);
            a *= fac;
            b *= fac;
            X[i] = a;
            Y[i] = b;
        }

        P2D aveG,varG;
        varG.x = compute_variance(X,aveG.x);
        varG.y = compute_variance(Y,aveG.y);
        double r_varG = 0;
        for(size_t i=1;i<=N;++i)
        {
            r_varG += Square(X[i]-ave.x) + Square(Y[i]-ave.y);
        }
        r_varG/=N;

        std::cerr << "2D uniform:  average=" << ave << ", variance=" << var << std::endl;
        std::cerr << "2D uniform   r_var  =" << r_var << std::endl;
        std::cerr << "2D gaussian: average=" << aveG << ", variance=" << varG << std::endl;
        std::cerr << "2D uniform   r_var  =" << r_varG << std::endl;


        P2D decr(floor(varG.x/var.x+0.5),floor(varG.y/var.y+0.5));
        std::cerr << "2D decrease=" << decr << ", r_decrease=" << floor(r_varG/r_var+0.5) << std::endl;
    }

    // 2D
    typedef point3d<double> P3D;
    {
        std::cerr << std::endl;
        vector<double> X(N),Y(N),Z(N);
        for(size_t i=1;i<=N;++i)
        {
            const P3D p = Ddt * P3D::on_unit_sphere(ran);
            X[i]  = p.x;
            Y[i]  = p.y;
            Z[i]  = p.z;
        }
        P3D ave,var;
        var.x = compute_variance(X,ave.x);
        var.y = compute_variance(Y,ave.y);
        var.z = compute_variance(Z,ave.z);
        double r_var = 0;
        for(size_t i=1;i<=N;++i)
        {
            r_var += Square(X[i]-ave.x) + Square(Y[i]-ave.y) + Square(Z[i]-ave.z);
        }
        r_var/=N;

        const size_t ng = (2*(3*N)+1)/2;
        vector<double> g(ng,as_capacity);
        for(size_t i=0;i<ng/2;++i)
        {
            double a,b;
            ran.gaussian(a,b);
            a *= fac;
            b *= fac;
            g.push_back(a);
            g.push_back(b);
        }
        assert(g.size()==ng);

        for(size_t i=1,j=1;i<=N;++i)
        {
            X[i] = g[j++];
            Y[i] = g[j++];
            Z[i] = g[j++];
        }
        P3D aveG,varG;
        varG.x = compute_variance(X,aveG.x);
        varG.y = compute_variance(Y,aveG.y);
        varG.z = compute_variance(Z,aveG.z);
        double r_varG = 0;
        for(size_t i=1;i<=N;++i)
        {
            r_varG += Square(X[i]-ave.x) + Square(Y[i]-ave.y) + Square(Z[i]-ave.z);
        }
        r_varG/=N;



        std::cerr << "3D uniform:  average=" << ave << ", variance=" << var << std::endl;
        std::cerr << "3D uniform   r_var="   << r_var << std::endl;
        std::cerr << "3D gaussian: average=" << aveG << ", variance=" << varG << std::endl;
        std::cerr << "3D gaussian  r_var="   << r_varG << std::endl;

        P3D decr(floor(varG.x/var.x+0.5),floor(varG.y/var.y+0.5),floor(varG.z/var.z+0.5));
        std::cerr << "3D decrease=" << decr << std::endl;
    }



}
YOCTO_PROGRAM_END()


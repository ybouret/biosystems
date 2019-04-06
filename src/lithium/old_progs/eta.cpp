#include "li-common.hpp"
#include "y/program.hpp"

class Data
{
public:
    Variables vars;
    vector<double> aorg;

    Data() : vars(), aorg()
    {
        vars << "pH_ini" << "pH_end";
        vars << "pH_eta" << "pw_eta";


        aorg.make( vars.size(), 0);

        vars(aorg,"pH_eta") = 6.39;
        vars(aorg,"pw_eta") = 1.7;

        vars(aorg,"pH_ini") = 5.8;
        vars(aorg,"pH_end") = 6.8;
    }

    double get_eta( double h ) const
    {
        return __lithium::eta(h, aorg, vars);
    }

    double h_ini() const
    {
        return __lithium::h_ini(aorg,vars);
    }

    double h_end() const
    {
        return __lithium::h_end(aorg,vars);
    }

    double get_eta_ini() const
    {
        return get_eta( h_ini() );
    }

    double get_eta_end() const
    {
        return get_eta( h_end() );
    }

    void save_info(double pH_min, double pH_max)
    {
        if(pH_min>=pH_max) cswap(pH_max,pH_min);

        {
            vars(aorg,"pH_ini") = pH_min;
            vars(aorg,"pH_end") = pH_max;
            const string file_name = vformat("eta%1.1fto%1.1f.dat",pH_min,pH_max);
            ios::ocstream fp(file_name);
            const size_t  N  = 200;

            for(size_t i=0;i<=N;++i)
            {
                const double u   = double(i)/N;
                const double h   = h_ini() + u * (h_end()-h_ini());
                const double pH  = -log10(h);
                const double eta = get_eta(h );
                fp("%g %g %g %g\n",eta);
            }
        }

    }

    ~Data() throw()
    {
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Data);
};

Y_PROGRAM_START()
{
    Data data;

    const double pH_end[] = { 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2 };

    for(unsigned i=1;i<sizeof(pH_end)/sizeof(pH_end[0]);++i)
    {
        data.save_info(5.8,pH_end[i]);
    }

    std::cerr << "plot 'eta5.8to6.0.dat' w l, 'eta5.8to6.2.dat' w l, 'eta5.8to6.4.dat' w l, 'eta5.8to6.6.dat' w l, 'eta5.8to6.8.dat' w l, 'eta5.8to7.0.dat' w l, 'eta5.8to7.2.dat' w l" << std::endl;


}
Y_PROGRAM_END()


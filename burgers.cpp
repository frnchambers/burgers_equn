#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include <Derivatives/Types.hpp>
#include <Derivatives/Grid.hpp>
#include <Derivatives/Lagrange.hpp>

#include <boost/numeric/odeint.hpp>




template <size_t N>  // ------------------------------------------------------------------------- //
class burgers_equn {  // ------------------------------------------------------------------------ //
public:

// numerical accuracies
  const double
    epsrel=1.0e-6,
    epsabs=0.0;

// grid quantities
  Grid_base<N> &grid;
  lagrange<N> diff;
  algebra::vector dudx, d2udx2;

// parameters in solution
  const double nu=1.0;

  burgers_equn ( Grid_base<N> &grid_in, double nu_in )
    : grid(grid_in),
      diff(grid), dudx(N), d2udx2(N),
      nu(nu_in)
  {}

  size_t size() const {
    return grid.size();
  }

  algebra::vector init_u () {
    algebra::vector u_pts(N);
    for ( size_t i=0; i<N; ++i )
      u_pts[i] = std::sin( M_PI * grid[i] );
    return u_pts;
  }

  void operator() ( const algebra::vector &u , algebra::vector &dudt , double t ) {

    dudx   = diff.derivs(u);
    d2udx2 = diff.derivs2(u);

  // boundary points constant
    dudt[0] = 0.0;
    dudt[N-1] = 0.0;
    
  // inner points
    for ( size_t i=1; i<N-1; ++i )
      dudt[i] = -u[i]*dudx[i] + nu*d2udx2[i];
  }

};  // ------------------------------------------------------------------------------------------ //


template <size_t N>
class observer {  // ---------------------------------------------------------------------------- //
public:

  const burgers_equn<N> &ode;
  std::ofstream &out;

  observer ( const burgers_equn<N> &ode_in, std::ofstream &out_in )
    : ode(ode_in), out(out_in)
  {}

  void operator() ( const algebra::vector &u, double t ) {
    static int count=0, every=100;

    std::cout << "t = " << t << "\r";

    if ( count % every == 0 ) {
      out << "\"t = " << t << "\"\n";
      for ( size_t i=0; i<ode.size(); ++i )
        out << ode.grid[i] << ' ' << u[i] << ' ' << ode.dudx[i] << ' ' << ode.d2udx2[i] << '\n';
      out << "\n" << std::endl;
    }

    count++;
  }

};  // ------------------------------------------------------------------------------------------ //


int main( int argc, char * argv[] ) {
  std::cout << std::scientific << std::setprecision(5);

#if defined NGRID
  const size_t N_pts=NGRID;
#else
#warning no grid size specified, using default value of 21
  const size_t N_pts=21;
#endif

  const double a=-1.0, b=1.0, nu=1.0;

  std::cout << "# Initialising grid...\n";
  Cheb_2<N_pts> grid(a, b);

  std::cout << "# Initialising equation set-up...\n";
  burgers_equn<N_pts> burgers(grid, nu);
  
  std::cout << "# Initialising initial u...\n";
  algebra::vector u_pts = burgers.init_u();

  std::ofstream output("data/burgers.dat");
  output << std::scientific;
  observer<N_pts> obs(burgers,output);


  std::cout << std::scientific << std::setprecision(6) <<
    "# solution to buger's equation\n" <<
    "# N      = " << burgers.size() << '\n' <<
    "# nu     = " << burgers.nu << '\n' <<
    "# (a, b) = (" << a << ", " << b << ")\n";

  double t0=0.0, t1=1.0, dt_init=1.0e-9;
  {
    namespace ode = boost::numeric::odeint;
    size_t n_steps = integrate_adaptive(
                                        ode::make_controlled<
                                        ode::runge_kutta_dopri5<
                                        algebra::vector> >(
                                                           burgers.epsabs,
                                                           burgers.epsrel ),
                                        burgers,
                                        u_pts, t0, t1, dt_init,
                                        obs );
    obs(u_pts, t1);
    std::cout << "# N steps = " << n_steps << std::endl;
  }

}

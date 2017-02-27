#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include <utility>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint.hpp>

#include <Derivatives/Grid.hpp>
#include <Derivatives/Lagrange.hpp>



template <size_t N>  // ------------------------------------------------------------------------- //
class burgers_equn {  // ------------------------------------------------------------------------ //
public:

// numerical accuracies
  const double
    epsrel=1.0e-3,
    epsabs=0.0;
  
// parameters in solution
  double nu=1.0;

// grid quantities
  Grid_base<N> &grid;
  lagrange<N> diff;
  algebra::vector dudx, d2udx2;

  burgers_equn ( Grid_base<N> &grid_in, double nu_in )
    : nu(nu_in),
      grid(grid_in), diff(grid), dudx(N), d2udx2(N)
  {}

  size_t size() const {
    return grid.size();
  }

  algebra::vector init_u () {
    algebra::vector u_pts(N);
    for ( size_t i=0; i<N; ++i )
      u_pts[i] = - sin( grid[i] );
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



class observer {  // ---------------------------------------------------------------------------- //
public:

  std::vector<std::pair<algebra::vector, double> > &ut_save;

  observer ( std::vector<std::pair<algebra::vector, double> > &ut_in )
    : ut_save(ut_in)
  {}

  void operator() ( const algebra::vector &u, double t ) {
    std::cout << "\"t = " << t << "\r";
    ut_save.push_back( std::make_pair(u, t) );
  }

  template <class ode_type>
  void output ( const ode_type &ode, const std::string &filename ) {
    std::ofstream out(filename);
    out << std::scientific << std::setprecision(5) <<
      "# solution to buger's equation\n" <<
      "# N  = " << ode.size() << '\n' <<
      "# nu = " << ode.nu << '\n';
    for ( const auto pair : ut_save ) {
      const algebra::vector u(pair.first);

      out << "\"t = " << pair.second << "\"\n";
      for ( size_t i=0; i<ode.size(); ++i )
        out << ode.grid[i] << ' ' << u[i] << '\n'; // << ode.dudx[i] << ' ' << ode.d2udx2[i] << '\n';
      out << "e" << std::endl;

    }

  }

};  // ------------------------------------------------------------------------------------------ //


int main( int argc, char * argv[] ) {
  std::cout << std::scientific << std::setprecision(5);

  const size_t N_pts=21;
  const double a=-1.0, b=1.0, nu=1.0;

  std::cout << "# Initialising grid...\n";
  Cheb_2<N_pts> grid(a, b);

  std::cout << "# Initialising equation set-up...\n";
  burgers_equn<N_pts> burgers(grid, nu);

  std::cout << "# Initialising initial u...\n";
  algebra::vector u_pts = burgers.init_u();

  // std::cout << "# Printing:\n";
  // for ( size_t i=0; i<N_pts; ++i )
  //   std::cout << grid[i] << ' ' << u_pts[i] << '\n';

  std::vector<std::pair<algebra::vector, double> > ut_save;
  observer obs(ut_save);


  double t0=0.0, t1=1.0, dt_init=1.0e-3;
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

  obs.output(burgers, "data/burgers.dat");

}

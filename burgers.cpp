#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

#include <boost/numeric/odeint.hpp>
#include <boost/program_options.hpp>

#include <Derivatives/Types.hpp>
#include <Derivatives/Grid.hpp>
#include <Derivatives/Lagrange.hpp>



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

  double a, b, nu;
  std::string outfile;

  {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("LHS-boundary,a", po::value<double>(&a)->default_value(-1.0), "Left hand side boundary")
      ("RHS-boundary,b", po::value<double>(&b)->default_value(1.0), "Right hand side boundary")
      ("diffusion-coefficient,nu", po::value<double>(&nu)->default_value(1.0),
       "Diffusion coefficient in equation")
      ("output-file,O", po::value<std::string>(&outfile)->default_value(std::string("data/burgers.dat")),
       "File to output solution")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
    }
  }

  std::cout << "# Initialising grid...\n";
  Cheb_2<N_pts> grid(a, b);

  std::cout << "# Initialising equation and initial conditions...\n";
  burgers_equn<N_pts> burgers(grid, nu);
  algebra::vector u_pts = burgers.init_u();

  std::ofstream output(outfile.c_str());
  output << std::scientific << std::setprecision(6);
  observer<N_pts> obs(burgers,output);

  std::cout <<
    "# solution to burger's equation\n" <<
    "# N      = " << burgers.size() << '\n' <<
    "# nu     = " << burgers.nu << '\n' <<
    "# (a, b) = (" << a << ", " << b << ")\n" <<
    "# outputting to file: " << outfile << '\n';

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

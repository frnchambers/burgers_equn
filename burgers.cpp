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



////////////////////////////////////////////////////////////////////////////////////////////////////
template <size_t N>  // ------------------------------------------------------------------------- //
struct burgers_equn {  // ----------------------------------------------------------------------- //

// numerical accuracies for time integration
  const double
    epsrel=1.0e-6,
    epsabs=0.0;

// grid quantities
  const Grid_base<N> &grid;
  const lagrange<N> diff;

// space to save derivative vectors
  mutable algebra::vector dudx, d2udx2, ududx;

// parameters in solution
  const double nu=1.0;

  burgers_equn ( const Grid_base<N> &grid_in, const double &nu_in )
    : grid(grid_in),
      diff(grid), dudx(N), d2udx2(N), ududx(N),
      nu(nu_in)
  {}

  size_t size() const {
    return grid.size();
  }

  algebra::vector init_u () const {
    algebra::vector u_pts(N);
    for ( size_t i=0; i<N; ++i )
      // u_pts[i] = - std::cos( M_PI * grid[i] );
      u_pts[i] = std::sin( M_PI * grid[i] );
    return u_pts;
  }

  void operator() ( const algebra::vector &u , algebra::vector &dudt , double t ) const {

    diff.derivs(u, dudx);
    diff.derivs2(u, d2udx2);

    for ( size_t i=0; i<N; ++i ) {
      ududx[i] = 0.0;
      for ( size_t n=0; n<i; ++n ) {
        ududx[i] += u[n] * dudx[i-n];
      }
    }

  // boundary points constant
    dudt[0] = 0.0;
    dudt[N-1] = 0.0;
    
  // inner points
    for ( size_t i=1; i<N-1; ++i )
      dudt[i] = -ududx[i] + nu*d2udx2[i];
  }

};  // ------------------------------------------------------------------------------------------ //
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
template <size_t N>  // ------------------------------------------------------------------------- //
struct observer {  // --------------------------------------------------------------------------- //

  const burgers_equn<N> &ode;
  std::ofstream &out;

  observer ( const burgers_equn<N> &ode_in, std::ofstream &out_in )
    : ode(ode_in), out(out_in)
  {}

  void operator() ( const algebra::vector &u, const double &t ) {
    static int count=0, every=10;

    std::cout << "t = " << t << "\r";

    if ( count % every == 0 ) {
      out << "t: " << t << ' ' << ode.nu << '\n';

      ode.diff.derivs(u, ode.dudx);
      ode.diff.derivs2(u, ode.d2udx2);

      for ( size_t i=0; i<ode.size(); ++i )
        out << "d: " <<
          u[i] << ' ' <<
          ode.dudx[i] << ' ' << ode.d2udx2[i] <<
          '\n';
      out << "-----" << std::endl;

    }
    count++;
  }

};  // ------------------------------------------------------------------------------------------ //
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char * argv[] ) {  // ------------------------------------------------------- //
  std::cout << std::scientific << std::setprecision(5);

#if defined NGRID
  const size_t N_pts=NGRID;
#else
#warning no grid size specified, using default value of 21
  const size_t N_pts=21;
#endif

// -- Parameters used in solution, location to save solution -- //
  double a, b, nu, t1;
  std::string outfile, griddir;
// ------------------------------------------------------------ // 

  {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("LHS-boundary,a", po::value<double>(&a)->default_value(-1.0), "Left hand side boundary")
      ("RHS-boundary,b", po::value<double>(&b)->default_value(1.0), "Right hand side boundary")
      ("diffusion-coefficient,nu", po::value<double>(&nu)->default_value(1.0), "Diffusion coefficient in equation")
      ("t-end,t", po::value<double>(&t1)->default_value(1.0),"Time to integrate until")
      ("output-file,O", po::value<std::string>(&outfile)->default_value(std::string("data/solution.dat")), "File to output solution")
      ("grid-directory,G", po::value<std::string>(&griddir)->default_value(std::string("data/")), "Directory where to save grid data")
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
  const Cheb_2<N_pts> grid(a, b);

  std::cout << "# Initialising equation and initial conditions...\n";
  const burgers_equn<N_pts> burgers(grid, nu);
  algebra::vector u_pts = burgers.init_u();

  std::cout << "# Initialising output utilities...\n";
  burgers.diff.save_data(griddir);
  std::cout << "# outputting solution to file: " << outfile << "...\n";
  std::ofstream output(outfile);
  output << std::scientific << std::setprecision(6) <<
    "# solution to burger's equation\n" <<
    "# N_grid = "   << burgers.size() << '\n' <<
    "# nu     = "   << burgers.nu << '\n' <<
    "# (a,b)  = ( " << a << " , " << b << " )\n" <<
    "# t_end  = "   << t1 << '\n';
  observer<N_pts> obs(burgers,output);

  double t0=0.0, dt_init=1.0e-9;
  {
    namespace ode = boost::numeric::odeint;
    auto stepper = ode::make_controlled<ode::runge_kutta_dopri5<algebra::vector> >
                                                                   (burgers.epsabs, burgers.epsrel);
    size_t n_steps = integrate_adaptive( stepper,
                                         burgers,
                                         u_pts, t0, t1, dt_init,
                                         obs );
    obs(u_pts, t1);
    std::cout << "# N steps = " << n_steps << std::endl;
  }

}  // ------------------------------------------------------------------------------------------- //
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <functional>

#include <boost/program_options.hpp>

#include <Derivatives/Types.hpp>
#include <Derivatives/Grid.hpp>
#include <Derivatives/Lagrange.hpp>

#include <boost/numeric/odeint.hpp>

// #include <Derivatives/Blaze_range_algebra.hpp>
// #include <boost/numeric/odeint/external/blaze/blaze_resize.hpp>
// #include <boost/numeric/odeint/external/blaze/blaze_algebra_dispatcher.hpp>


////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
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

  void operator() ( const algebra::vector &u , algebra::vector &dudlogt , double t ) const {

    diff.derivs(u, dudx);
    diff.derivs2(u, d2udx2);

    for ( size_t i=0; i<N; ++i ) {
      ududx[i] = 0.0;
      for ( size_t n=0; n<i; ++n ) {
        ududx[i] += u[n] * dudx[i-n];
      }
    }

  // boundary points constant
    dudlogt[0] = 0.0;
    dudlogt[N-1] = 0.0;
    
  // inner points
    for ( size_t i=1; i<N-1; ++i )
      dudlogt[i] = -ududx[i] + nu*d2udx2[i];
  }

};  // ------------------------------------------------------------------------------------------ //
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
template <size_t N>  // ------------------------------------------------------------------------- //
struct observer {  // --------------------------------------------------------------------------- //

  const burgers_equn<N> &ode;
  std::ofstream &out;
  const int every;

  observer ( const burgers_equn<N> &ode_in, const std::string &griddir,
             const double t_end, const double dt, const int n_prints,
             std::ofstream &out_in )
    : ode(ode_in), out(out_in), every( static_cast<int>(t_end/(dt*static_cast<double>(n_prints)) ) )
  {
    out << std::scientific << std::setprecision(6) <<
      "# solution to burger's equation\n" <<
      "# N_grid = "   << N << '\n' <<
      "# nu     = "   << ode.nu << '\n' <<
      "# (a,b)  = ( " << ode.grid[0] << " , " << ode.grid[N-1] << " )\n" <<
      "# t_end  = "   << t_end << '\n' <<
      "# dt     = "   << dt << '\n' <<
      "# every  = "   << every << '\n';
    ode.diff.save_data(griddir);
  }

  void operator() ( const algebra::vector &u, const double &t ) {
    static int count=0;

    std::cout << "t = " << t << ", count = " << count << ", saves = " << count/every << "\r";

    if ( count % every == 0 ) {
      out << "t: " << t << ' ' << count << '\n';

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
// ---------------------------------------------------------------------------------------------- //
int main( int argc, char * argv[] ) {  // ------------------------------------------------------- //
  std::cout << std::scientific << std::setprecision(5);

#if defined NGRID
  const size_t N_pts=NGRID;
#else
#warning no grid size specified, using default value of 31
  const size_t N_pts=31;
#endif

// -- Parameters used in solution, location to save solution -- //
  double a, b, nu, t1, dt;
  int n_prints;
  std::string outdir;
// ------------------------------------------------------------ // 

  {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("LHS-boundary,a", po::value<double>(&a)->default_value(-1.0), "Left hand side boundary")
      ("RHS-boundary,b", po::value<double>(&b)->default_value(1.0), "Right hand side boundary")
      ("ode-coeff,nu", po::value<double>(&nu)->default_value(1.0), "Diffusion coefficient in equation")
      ("dt,s", po::value<double>(&dt)->default_value(1.0e-5),"Initial time step")
      ("t-end,t", po::value<double>(&t1)->default_value(1.0),"Time to integrate until")
      ("n-prints,n", po::value<int>(&n_prints)->default_value(100),"Number of steps to print out")
      ("solun-directory,D", po::value<std::string>(&outdir)->default_value(std::string("data/")), "Directory where to save solution and grid data")
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
  std::ofstream output(outdir+"solution.dat");
  observer<N_pts> obs(burgers, outdir, t1, dt, n_prints, output);
  std::cout << "# -> outputting solution to file: " << outdir << "solution.dat...\n";

  {
    namespace ode = boost::numeric::odeint;

    typedef ode::adams_bashforth_moulton<5, algebra::vector> stepper_type;
    //auto stepper = ode::make_controlled( burgers.epsabs, burgers.epsrel, stepper_type() );
    size_t n_steps = ode::integrate_const( stepper_type(), burgers, u_pts, 0.0, t1, dt, obs );

    obs(u_pts, t1);
    std::cout << "# N steps = " << n_steps << std::endl;
  }

}  // ------------------------------------------------------------------------------------------- //
////////////////////////////////////////////////////////////////////////////////////////////////////

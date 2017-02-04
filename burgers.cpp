#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint.hpp>

#include <Grid/Grid.hpp>
#include <Interpolation/Interp_lagrange.hpp>


typedef lagrange_types::vec_type vec_type;


template <class grid_type>  // ----------------------------------------------------------------------- //
class burgers_equn {  // ----------------------------------------------------------------------------- //
public:

// numerical accuracies
  const double
    epsrel=1.0e-6,
    epsabs=0.0;
  
// parameters in solution
  double nu=1.0;

// grid quantities
  size_t N_pts;
  lagrange<grid_type> diff;
  vec_type dudx, d2udx2;

  burgers_equn ( size_t N_in, const double nu_in )
    : nu(nu_in), N_pts(N_in), diff(N_in), dudx(N_in), d2udx2(N_in)
  {}

  vec_type init_u () {
    vec_type u_pts(N_pts);
    for ( size_t i=0; i<N_pts; ++i )
      u_pts[i] = - sin( diff[i] );
    return u_pts;
  }


  void operator() ( const vec_type &u , vec_type &dudt , double t ) {

    dudx   = diff.derivs(u);
    d2udx2 = diff.derivs2(u);

  // boundary points constant
    dudt[0] = 0.0;
    dudt[N_pts-1] = 0.0;
    
  // inner points
    for ( size_t i=1; i<N_pts-1; ++i )
      dudt[i] = -u[i]*dudx[i] + nu*d2udx2[i];
  }

};  // ------------------------------------------------------------------------------------------ //


// class observer {  // ---------------------------------------------------------------------------- //
// public:

//   void operator() ( const vec_type &u, double t ) {
//     std::cout << "t = " << t << "\n";
//     // for ( size_t i=0; i<u.size(); ++i )
//     //   std::cout << u[i] << '\n';
//   }

// };  // ------------------------------------------------------------------------------------------ //

template <class ode_type>  // ------------------------------------------------------------------- //
class observer {  // ---------------------------------------------------------------------------- //
public:

  const ode_type & ode;
  std::ofstream out;

  observer ( const ode_type & ode_in, std::string filename )
    : ode(ode_in), out(filename)
  {
    out << std::scientific << std::setprecision(5) <<
      "# solution to buger's equation\n" <<
      "# N  = " << ode.N_pts << '\n' <<
      "# nu = " << ode.nu << '\n';
  }

  void operator() ( const vec_type &u, double t ) {
    out << "\"t = " << t << "\"\n";
    for ( size_t i=0; i<ode.N_pts; ++i )
      out << ode.diff[i] << ' ' << u[i] << ' ' << ode.dudx[i] << ' ' << ode.d2udx2[i] << '\n';
    out << "e" << std::endl;
  }

};  // ------------------------------------------------------------------------------------------ //


int main( int argc, char * argv[] ) {

  typedef burgers_equn<Cheb_1> ode_type;

  std::cout << std::scientific << std::setprecision(5);

  size_t N_pts=21;
  if ( argc>1 )
    N_pts = atoi(argv[1]);

  std::cout << "# initialised burgers" << std::endl;
  ode_type burgers( N_pts, 1.0 );
  vec_type u_pts = burgers.init_u();

  observer<ode_type> obs(burgers, "data/solun.dat");
  // observer obs(burgers,"data/solun.dat");

  double t0=0.0, t1=1.0, dt_init=1.0e-3;
  {
    namespace ode = boost::numeric::odeint;
    size_t n_steps = integrate_adaptive(
                                        ode::make_controlled<
                                        ode::runge_kutta_dopri5<vec_type> >(
                                                                            burgers.epsabs,
                                                                            burgers.epsrel ),
                                        burgers,
                                        u_pts, t0, t1, dt_init,
                                        obs );

    std::cout << n_steps << std::endl;
  }


}

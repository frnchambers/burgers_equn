#ifndef __DERIVATIVES_LAGRANGE_HPP__
#define __DERIVATIVES_LAGRANGE_HPP__

#include <Derivatives/Types.hpp>
#include <Derivatives/Grid.hpp>

////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
template <size_t N>  // -------------------------------------------------------- //
class lagrange { // ----------------------------------------------------------------------------- //
// Spectral Methods routines given in:
//          Berrut & Trefethen 2004 "Barycentric Lagrange Interpolation"
//                             - SIAM Review Vol. 46, No. 3, pp. 501-517
public:
//private:

////// grid structure
  const Grid_base<N> &grid;
////// matrices used to calculate derivatives at gridpoints
  algebra::matrix
    d1lpts,  // first derivative
    d2lpts;  // second derivative

public:

  lagrange ( const Grid_base<N> &grid_in )
    : grid(grid_in), d1lpts(N,N), d2lpts(N,N)
  {
    std::cout << "# Initialising Derivative space\n";
    update_derivs();
  }

  void update_derivs () {
    std::cout << "# First derivative\n";
    update_first_deriv();
    std::cout << "# Second derivative\n";
    update_second_deriv();
  }

  size_t size () const {
    return N;
  }

  // algebra::vector derivs  ( const algebra::vector &fpts ) const {
  //   // return algebra::prod( d1lpts, fpts );
  //   return prod( d1lpts, fpts );
  // }
  // algebra::vector derivs2 ( const algebra::vector &fpts ) const {
  //   // return algebra::prod( d2lpts, fpts );
  //   return prod( d2lpts, fpts );
  // }

  void derivs  ( const algebra::vector &fpts, algebra::vector &dfptsdx ) const {
    // return algebra::prod( d1lpts, fpts );
    dfptsdx = prod( d1lpts, fpts );
  }
  void derivs2 ( const algebra::vector &fpts, algebra::vector &d2fptsdx2 ) const {
    // return algebra::prod( d2lpts, fpts );
    d2fptsdx2 = prod( d2lpts, fpts );
  }

  double df_FD_d ( const size_t &i, const size_t &j ) const {
    return d1lpts(i,j);
  }
  double df_FD_dd ( const size_t &i, const size_t &j ) const {
    return d2lpts(i,j);
  }

private:

////// interpolate, at points x, arbitrary function on grid of interpolation points
// using formula (4.2) - the big one in the box
  double interpolate ( const algebra::vector &fpts, const double x ) const {
    double numer=0, denom=0;
    for( size_t j=0; j<N; j++ ) {
      numer += grid.w(j) / ( x - grid.x(j) ) * fpts[j];
      denom += grid.w(j) / ( x - grid.x(j) );
    }
    return numer / denom;
  }

// --------------------------------------------------------------------------------------- //
// --------- update the matrix d1lpts and d2lpts for derivatives at points grid ---------- //
// -------- equations (9.4) and (9.5), TYPO IN PAPER!!! use formulae given below --------- //
// --------------------------------------------------------------------------------------- //

////// first derivative
// D^{(1)}_{ij} =   l_j'(x_i)
//                            =   \frac{ w_j/w_i }{ x_i-x_j }  ... i \neq j
//                            = - \sum_{k \neq j} l_j'(x_k)    ... i  = j
  void update_first_deriv () {

  // off diagonal components, formula (9.4)
    for ( size_t i=0; i<N; i++ )
      for ( size_t j=0; j<N; j++ )
        if ( j!=i )
          d1lpts(i,j) = ( grid.w(j)/grid.w(i) ) / ( grid.x(i)-grid.x(j) );

  // diagonal components, formula (9.5)
    for ( size_t j=0; j<N; j++ )
      for ( size_t i=0; i<N; i++ )
        if ( i!=j )
          d1lpts(i,i) -= d1lpts(i,j);

  } // end update_first_deriv ()


////// second derivative
// D^{(2)}_{ij} =    l_j''(x_i)
//                            = -2 l_j'(x_i) \left( \frac{1}{x_i-x_j}
//                              + \sum_{k \neq i} l_k'(x_i) \right)  ... i \neq j
//                            = -  \sum_{k \neq j} l_j''(x_k)  ... i = j
  void update_second_deriv () {

  // off diagonal components, formula (9.4)
    for ( size_t i=0; i<N; i++ )
      for ( size_t j=0; j<N; j++ )
        if ( i!=j ) {

          double sum=0.0;
          for ( size_t k=0; k<N; k++ )
            if ( k!=i )
              sum += d1lpts(i,k);

          d2lpts(i,j) = -2.0 * d1lpts(i,j) * ( sum + ( 1.0/( grid.x(i)-grid.x(j) ) ) );
        }

  // diagonal components, formula (9.5)
    for ( size_t j=0; j<N; j++ )
      for ( size_t i=0; i<N; i++ )
        if ( i!=j )
          d2lpts(j,j) -= d2lpts(j,i);
    
  }  // end update_deriv ()

public:

////// print grid, wpts and derivative matrices
  void save_data ( const std::string & dirname ) const {
    std::cout <<
      "# Saving grid data to files:\n" <<
      "# -> grid.dat, weights.dat, d1_lpts.dat, and d2_lpts.dat\n" <<
      "# in directory: " << dirname << std::endl;

    std::ofstream
      xpts(dirname+"grid.dat"),
      wpts(dirname+"weights.dat"),
      d1_matrix(dirname+"d1_lpts.dat"),
      d2_matrix(dirname+"d2_lpts.dat");

    xpts << std::scientific;
    wpts << std::scientific;
    d1_matrix << std::scientific;
    d2_matrix << std::scientific;
    
    for ( size_t i=0; i<N; i++ ) {
      xpts << grid.x(i) << '\n';
      wpts << grid.w(i) << '\n';
      for ( size_t j=0; j<N; j++ ) {
        d1_matrix << d1lpts(i,j) <<  ' ';
        d2_matrix << d2lpts(i,j) << ' ';
      }
      d1_matrix << '\n';
      d2_matrix << '\n';
    }

    xpts << std::flush;
    wpts << std::flush;
    d1_matrix << std::flush;
    d2_matrix << std::flush;
  } // end print_data

}; // --------------------------------------------------------------------------------------------------------- //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////






// double operator() ( const algebra::vector & fpts, const double & x ) const {
//   return interpolate( fpts, x );
// }



#endif

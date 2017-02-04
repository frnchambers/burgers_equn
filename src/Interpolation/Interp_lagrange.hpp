#ifndef __INTERP_LAGRANGE_H__
#define __INTERP_LAGRANGE_H__

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <Grid/Grid.hpp>

////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
namespace lagrange_types {
  typedef boost::numeric::ublas::vector< double > vec_type;
  typedef boost::numeric::ublas::matrix< double > mat_type;
}

template <class grid_type>
class lagrange { // ----------------------------------------------------------------------------- //
  // Spectral Methods routines given in:
  //          Berrut & Trefethen 2004 "Barycentric Lagrange Interpolation"
  //                             - SIAM Review Vol. 46, No. 3, pp. 501-517
public:

  typedef lagrange_types::vec_type vec_type;
  typedef lagrange_types::mat_type mat_type;

  //private:

////// datapoints for interpolation
  size_t N_points;  
////// grid structure
  grid_type xpts;
////// matrices used to calculate derivatives at gridpoints
  mat_type
    d1lpts,  // first derivative
    d2lpts;  // second derivative

  const std::string name="lagrange";

public:

  lagrange ( size_t N_in )
    : N_points( N_in ), xpts(N_in),
      d1lpts( N_in,N_in,0.0 ), d2lpts( N_in,N_in,0.0 )
  {
    update_derivs();
  }

  lagrange ( grid_type & xpts_in )
    : N_points( xpts_in.size() ), xpts(xpts_in),
      d1lpts( xpts_in.size(),xpts_in.size(),0.0 ), d2lpts( xpts_in.size(),xpts_in.size(),0.0 )
  {
    update_derivs();
  }

  size_t size () {
    return N_points;
  }
  double operator[] ( const size_t & i ) const {
    return xpts[i];
  }
  double x ( const size_t & i ) const {
    return xpts[i];
  }
  std::string diff_name () const {
    return name;
  }
  std::string grid_name () const {
    return xpts.name;
  }

  ////// interpolate, at points x, arbitrary function on grid of interpolation points
  double operator() ( const vec_type & fpts, const double & x ) const {
    return interpolate( fpts, x );
  }

  vec_type derivs ( const vec_type & fpts ) const {
    return prod( d1lpts, fpts );
  }
  vec_type derivs2 ( const vec_type & fpts ) const {
    return prod( d2lpts, fpts );
  }

  double df_FD_d ( const size_t &i, const size_t &j ) const {
    return d1lpts(i,j);
  }
  double df_FD_dd ( const size_t &i, const size_t &j ) const {
    return d2lpts(i,j);
  }

  void update_derivs () {
    update_first_deriv();
    update_second_deriv();
  }


  ////// print grid, wpts and derivative matrices
  void print_data () const {

    std::cout << "points:\n";
    for ( size_t i=0; i<N_points; i++ )
      std::cout <<
        xpts[i] << "\t" <<
        std::endl;

    std::cout << "weights:\n";
    for ( size_t i=0; i<N_points; i++ )
      std::cout <<
        xpts.w(i) << "\t" <<
        std::endl;

    std::cout << "d1lpts:\n";
    for ( size_t i=0; i<N_points; i++ ) {
      for ( size_t j=0; j<N_points; j++ )
        std::cout << d1lpts(i,j) << "\t";
      std::cout << std::endl;
    }

    std::cout << "d2lpts:\n";
    for ( size_t i=0; i<N_points; i++ ) {
      for ( size_t j=0; j<N_points; j++ )
        std::cout << d2lpts(i,j) << "\t";
      std::cout << std::endl;
    }


  } // end print_data




private:

  ////// interpolate, at points x, arbitrary function on grid of interpolation points
  // using formula (4.2) - the big one in the box
  double interpolate ( const vec_type & fpts, const double x ) const {
    double numer=0, denom=0;

    for( size_t j=0; j<N_points; j++ ) {
      numer += xpts.w(j) / ( x - xpts[j] ) * fpts[j];
      denom += xpts.w(j) / ( x - xpts[j] );
    }

    return numer / denom;
  }


  // --------------------------------------------------------------------------------------- //
  // --------- update the matrix d1lpts and d2lpts for derivatives at points xpts ---------- //
  // -------- equations (9.4) and (9.5), TYPO IN PAPER!!! use formulae given below --------- //
  // --------------------------------------------------------------------------------------- //

  ////// first derivative
  // D^{(1)}_{ij} =   l_j'(x_i)
  //                            =   \frac{ w_j/w_i }{ x_i-x_j }  ... i \neq j
  //                            = - \sum_{k \neq j} l_j'(x_k)    ... i  = j
  void update_first_deriv () {

    //// off diagonal components, formula (9.4)
    for ( size_t i=0; i<N_points; i++ )
      for ( size_t j=0; j<N_points; j++ )
        if ( j!=i )
          d1lpts(i,j) = ( xpts.w(j)/xpts.w(i) ) / ( xpts[i]-xpts[j] );

    //// diagonal components, formula (9.5)
    for ( size_t j=0; j<N_points; j++ ) {
      for ( size_t i=0; i<N_points; i++ ) {
        if ( i!=j ) {
          d1lpts(i,i) -= d1lpts(i,j);
        }
      }
    }
  } // end update_first_deriv ()


  ////// second derivative
  // D^{(2)}_{ij} =    l_j''(x_i)
  //                            = -2 l_j'(x_i) \left( \frac{1}{x_i-x_j} + \sum_{k \neq i} l_k'(x_i) \right)  ... i \neq j
  //                            = -  \sum_{k \neq j} l_j''(x_k)  ... i = j
  void update_second_deriv () {

    //// off diagonal components, formula (9.4)
    for ( size_t i=0; i<N_points; i++ ) {
      for ( size_t j=0; j<N_points; j++ ) {
        if ( i!=j ) {

          double sum=0.0;
          for ( size_t k=0; k<N_points; k++ ) {
            if ( k!=i ) {
              sum += d1lpts(i,k);
            }
          }

          d2lpts(i,j) = -2.0 * d1lpts(i,j) * ( sum + ( 1.0/( xpts[i]-xpts[j] ) ) );
        }
      }
    }

    //// diagonal components, formula (9.5)
    for ( size_t j=0; j<N_points; j++ ) {
      for ( size_t i=0; i<N_points; i++ ) {
        if ( i!=j ) {
          d2lpts(j,j) -= d2lpts(j,i);
        }
      }
    }

  } // end update_deriv ()


}; // --------------------------------------------------------------------------------------------------------- //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////





#endif

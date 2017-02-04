#ifndef __INTERP_SPLINE_HPP__
#define __INTERP_SPLINE_HPP__

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

template <class x_type=boost::numeric::ublas::vector<double> >
class spline_interpolate {
public:

  typedef boost::numeric::ublas::vector< double > vec_type;
  typedef boost::numeric::ublas::matrix< double > mat_type;

  size_t N_points;
  x_type xpts;
  vec_type ypts;
  vec_type y2ppts;

  std::string name = "spline";
  std::string interp_name () {
    return name;
  }

  ////// value of derivative at boundary, set when performing update
  double yp0=0.0,ypNm1=0.0;

  spline_interpolate ( const x_type &xpts_new, const vec_type &ypts_new )
    : N_points(xpts_new.size()), xpts(xpts_new), ypts(ypts_new), y2ppts(xpts_new.size()) {
  }

  ////// often require this expression, make a function
  double rat ( size_t j ) {
    // if ( j>N_points-2 ) {
    // 	std::cerr << " -> bad intex in spline_interpolation: j=" << j << std::endl;;
    // 	throw(0);
    // }
    return ( ypts[j+1]-ypts[j] )/( xpts[j+1]-xpts[j] );
  }


  ////// find j for which: x_j < x < x_{j+1}
  size_t locate( double x ) {

    if ( x<xpts[0] || x>xpts[N_points-1] ) {
      std::cerr << " -> x=" << x << " out of range: (" << xpts[0] << "," << xpts[N_points-1] << ")" << std::endl;
      throw(0);
    }

    size_t jup=N_points-1,jlo=0,jmid=0;;
    while ( jup > 1+jlo ) {

      // std::cout <<
      // 	jlo << ' ' << jup << ' ' << jup-jlo << ' ' << jmid << std::endl;

      jmid = 0.5*(jup+jlo);
      if ( x >= xpts[jmid] )
        jlo=jmid;
      else
        jup=jmid;
    }
    return jmid;
  }

  ////// perform interpolation
  double operator() ( double x ) {

    size_t j = locate(x);

    double A = ( xpts[j+1]-x )/( xpts[j+1]-xpts[j] );
    double B = 1.0-A;

    double fact = ( xpts[j]-xpts[j] ) * ( xpts[j]-xpts[j] ) / 6.0;
    double
      C = A*( A*A - 1.0 ) * fact,
      D = B*( B*B - 1.0 ) * fact;

    return A*ypts[j] + B*ypts[j+1] + C*y2ppts[j] + D*y2ppts[j+1];
  }

  ////// perform interpolation for first derivative
  double deriv ( double x ) {

    size_t j = locate(x);

    double A = ( xpts[j+1]-x )/( xpts[j+1]-xpts[j] );
    double B = 1.0-A;

    return
      rat(j) +
      ( xpts[j+1]-xpts[j] ) / 6.0 * (
                                     ( (3.0*B*B)-1.0 )*y2ppts[j+1] - ( (3.0*A*A)-1.0 )*y2ppts[j]
                                     );
  }



  ////// update y'' using "natural boundary conditions"
  void update () {
    update_y2p(0,0);
  }

  ////// update y'' using knowledge of derivatives at boundary
  void update ( double yp0_in, double ypNm1_in ) {
    yp0   = yp0_in;
    ypNm1 = ypNm1_in;

    double A0   = 3.0 / (xpts[1]-xpts[0])           * ( rat(0) - yp0 );
    double ANm1 = 3.0 / (xpts[N_points-1]-xpts[N_points-2]) * ( ypNm1 - rat(N_points-2) );

    update_y2p(A0,ANm1);
  }

  ////// update y'', must give boundary points A(0,0),A(N-1,N-1)
  void update_y2p ( double A0, double ANm1 ) {
    mat_type A(N_points,N_points);

    // boundary points
    A(0,0)           = A0;
    A(N_points-1,N_points-1) = ANm1;

    // inside matrix points
    for ( size_t j=1; j<N_points-1; ++j ) {
      A(j,j-1) = ( xpts[j]   - xpts[j-1] ) / 6.0;
      A(j,j)   = ( xpts[j+1] - xpts[j-1] ) / 3.0;
      A(j,j+1) = ( xpts[j+1] - xpts[j]   ) / 6.0;

      y2ppts[j] = rat(j) - rat(j-1);
    }

    // Perform lu decomposition and invert, in place
    {
      using namespace boost::numeric::ublas;
      permutation_matrix<size_t> pm(A.size1());
      lu_factorize(A, pm);
      lu_substitute(A, pm, y2ppts);
    }

  }


  vec_type derivs () {
    vec_type yppts(N_points);

    yppts[0]=yp0;
    yppts[N_points-1]=ypNm1;
    for ( size_t i=1; i<N_points-1; ++i ) {
      yppts[i] = rat(i) - ( xpts[i+1]-xpts[i] ) * ( 2.0*y2ppts[i] + y2ppts[i+1] ) / 6.0;
    }

    return yppts;
  }

  vec_type derivs2 () {
    return y2ppts;
  }


};

#endif

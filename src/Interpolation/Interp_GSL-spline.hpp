#ifndef __INTERP_GSL_SPLINE_HPP__
#define __INTERP_GSL_SPLINE_HPP__

#include <iostream>
#include <string>
// #include <boost/numeric/ublas/vector.hpp>
// typedef boost::numeric::ublas::vector< double > vec_type;

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


template<class vec_type>
class GSL_spline {

private:
  gsl_interp_accel *acc;
  gsl_spline *spline;

public:

  size_t N_points=0;
  std::string name = "GSL-spline";
  std::string interp_name() const {
    return name;
  }

  GSL_spline ()
  {
    acc = gsl_interp_accel_alloc();
  }

  GSL_spline ( const vec_type & xpts, const vec_type & fpts ) :
    N_points(xpts.size())
  {
    acc = gsl_interp_accel_alloc();
    set_points(xpts,fpts);
  }

  // ~GSL_spline () {
  //   gsl_spline_free (spline);
  //   gsl_interp_accel_free (acc);
  // }

  void set_points ( const vec_type & xpts, const vec_type & fpts ) {
    if ( xpts.size() != fpts.size() ) {
      std::cerr << "ERROR IN GSL SPLINE, mismatch in size" << std::endl;
      throw(0);
    }
    N_points=xpts.size();

    spline = gsl_spline_alloc( gsl_interp_cspline, N_points );
    gsl_spline_init( spline, &xpts[0], &fpts[0], N_points );
  }

  double operator() ( const double & x ) const {
    return gsl_spline_eval( spline, x, acc );
  }
  double deriv  ( const double & x ) const {
    return gsl_spline_eval_deriv( spline, x, acc );
  }
  double deriv2 ( const double & x ) const {
    return gsl_spline_eval_deriv2( spline, x, acc );
  }

  void operator() ( const double & x,
                    double &f, double &dfdx ) const {
    f    = gsl_spline_eval( spline, x, acc );
    dfdx = gsl_spline_eval_deriv( spline, x, acc );
  }

};

#endif

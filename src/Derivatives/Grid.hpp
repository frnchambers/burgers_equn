// Copyright 2016 frnchambers

#ifndef __DERIVATIVES_GRID_HPP__
#define __DERIVATIVES_GRID_HPP__

#include <cmath>
#include <Derivatives/Types.hpp>



template <typename vec_type>
void norm_vec (vec_type &v, double norm=1.0) {  // normalise vector so that largest value is norm
  double v_biggest=0.0;
  for ( const auto & vi : v )
    if ( std::abs(vi) > std::abs(v_biggest) )
      v_biggest = vi;
  for ( double &vi : v )
    vi *= norm / v_biggest;
}

template <typename vec_type>
void w_from_x (vec_type &w, const vec_type &x) {  // for grids that don't have a simple analytic w
  size_t N=w.size();
  for ( size_t i=0; i < N; ++i ) {
    w[i] = 1.0;
    for ( size_t k=0; k < N; ++k ) {
      if (k != i)
        w[i] *= x[i]-x[k];
    }
    w[i] = 1.0/w[i];
  }
  norm_vec(w);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
template <size_t N>  // ------------------------------------------------------------------------- //
class Grid_base {  // -------------------------------------------- Base Class to define a grid -- //
protected:
  grid_types::vector<N> xpts, wpts;
  // const double a = -1.0, x0 = 0.0, b = 1.0;  // end points and 'center' of the grid 
  virtual double x_func(const size_t & i) const = 0;  // depends on the grid
  virtual double w_func(const size_t & i) const = 0;  // depends on the grid
  // Grid_base ( double a_in, double x0_in, double b_in )
  //   : a(a_in), x0(x0_in), b(b_in)
  // {}
  void auto_grid() {
    for (size_t i = 0; i < N; ++i)
      xpts[i] = x_func(i);
    w_from_x(wpts, xpts);
  }
  void set_grid() {
    for (size_t i = 0; i < N; ++i) {
      xpts[i] = x_func(i);
      wpts[i] = w_func(i);
    }
  }
public:  // ------------------------------------------------------------------------------------- //
////// return size of grid
  size_t size() const {
    return N;
  }
////// call element i of x or w grid
  double operator[] (const size_t & i) const {
    return xpts[i];
  }
  double &w (const size_t & i) {
    return wpts[i];
  }
  double &x (const size_t & i) {
    return xpts[i];
  }
}; // ------------------------------------------------------------------------------------------- //
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
template <size_t N>  // ------------------------------------------------------------------------- //
class Even : public Grid_base<N> {  // ------------------------------------- Evenly space grid -- //
  using Base = Grid_base<N>;
private:
  const double a = -1.0, x0 = 0.0, b = 1.0;  // end points and 'center' of the grid 
////// default parameters for grid centered at 0, with boundaries: -1,1
  double x_func (const size_t &i) const {
    return a + ( b-a ) * static_cast<double>(i) / static_cast<double>(N-1);
  }
  double w_func (const size_t &i) const {
    return 0;
  }
 public:
  Even ( double a_in = -1.0, double x0_in = 0.0, double b_in = 1.0 )
    : a(a_in), x0(x0_in), b(b_in)
  {
    Base::auto_grid();
  }
  // void reset(double a_in, double x0_in, double b_in) {
  //   a = a_in; x0 = x0_in; b = b_in;
  //   Grid_base<N>::auto_grid();
  // }
};  // ------------------------------------------------------------------------------------------ //
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
template <size_t N>  // ------------------------------------------------------------------------- //
class Sinh : public Grid_base<N> {  // -------------------------- Shiny log -> A.L. Watts 2016 -- //
  using Base = Grid_base<N>;
private:
////// default parameters for grid centered at 0, with boundaries: -1,1
  const double
    a = -1.0, x0 = 0.0, b = 1.0,
    d = -1.0, m = 2.0;
  //  formulae for x0, m, and d: x = x0 + sinh( d + m * i/(N-1) )
  //  for x 'element of' [a,b]:
  //                -> d = asinh( a-x0 )
  //                -> m = asinh( b-x0 ) - asinh( a-x0 )
  double x_func(const size_t & i) const {
    double z = d + m * static_cast<double>(i) / static_cast<double>(N-1);
    return x0 + std::sinh( z );
  }
  double w_func (const size_t & i) const {
    return 0;
  }
////// set up parameters that characterise the grid
  double d_param (const double & a_in, const double & x0_in, const double & b_in) {
    return std::asinh(a_in - x0_in);
  }
  double m_param (const double & a_in, const double & x0_in, const double & b_in) {
    return std::asinh(b_in - x0_in) - std::asinh(a_in - x0_in);
  }
public:
  Sinh ( double a_in = -1.0, double x0_in = 0.0, double b_in = 1.0)
    : a(a_in), x0(x0_in), b(b_in),
      d(d_param(a,x0,b)), m(m_param(a,x0,b))
  {
    Grid_base<N>::auto_grid();
  }
  // void reset(double a_in, double x0_in, double b_in) {
  //   a = a_in; x0 = x0_in; b = b_in;
  //   m=m_param(a, x0, b);
  //   d=d_param(a, x0, b);
  //   auto_grid();
  // }
};  // ------------------------------------------------------------------------------------------ //
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
template <size_t N>  // ------------------------------------------------------------------------- //
class Cheb_1 : public Grid_base<N> {  // ----------------------------------- gauss-lobato mesh -- //
  using Base = Grid_base<N>;
private:
  const double a = -1.0, b = 1.0;
  double x_func(const size_t & i) const {
    double t = std::cos(M_PI * static_cast<double>(2*i + 1)
                        / static_cast<double>(2*N) );
    return 0.5*(a+b) - 0.5*(a-b)*t;
  }
  double w_func(const size_t & i) const {
    return ( i%2 == 1 ? -1.0 : 1.0 ) *
        std::sin(M_PI * static_cast<double>(2*i + 1)
                 / static_cast<double>(2*N) );
  }
 public:
  Cheb_1(double a_in = -1.0, double b_in = 1.0)
    : a(a_in), b(b_in)
  {
    Base::set_grid();
  }
  // void reset(double a_in, double b_in) {
  //   a = a_in; b = b_in;
  //   set_grid();
  // }
};  // ------------------------------------------------------------------------------------------ //
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
template <size_t N>  // ------------------------------------------------------------------------- //
class Cheb_2 : public Grid_base<N> {  // -------------------------------------------------------- //
  using Base = Grid_base<N>;
private:
  const double a = -1.0, b = 1.0;
  double x_func(const size_t & i) const {
    double t = std::cos(M_PI * static_cast<double>(i)
                        / (static_cast<double>(N-1) ) );
    return 0.5*(a+b) + 0.5*(a-b)*t;
  }
  double w_func(const size_t & i) const {
    return ( i == 0 || i == N-1 ? 0.5 : 1.0 )
        * (i%2 == 1 ? -1.0 : 1.0);
  }
 public:
  Cheb_2(double a_in = -1.0, double b_in = 1.0)
    : a(a_in), b(b_in)
  {
    Base::set_grid();
  }
  // void reset(double a_in, double x0_in, double b_in) {
  //   x0_in=x0_in;
  //   a = a_in; b = b_in;
  //   set_grid();
  // }
};  // ------------------------------------------------------------------------------------------ //
////////////////////////////////////////////////////////////////////////////////////////////////////










// // ---------------------------------------------------------------------------------------------- //
// class Arb : public Grid {  // -- For arbitrary grid --------------------------------------------- //
// public:

//  Arb ( const vec_type & xpts_in, const std::string name_in="" ) : xpts(xpts_in), name(name_in) {
//    a=xpts[0];
//    b=xpts[N_points-1];   
// ////// initialise weights
//    w_from_x();
//  }
// private:
//  double x_func ( size_t i ) const {
//    return xpts[i];
//  }
//  double w_func ( size_t i ) const {
//    return wpts[i];
//  }
// };  // ------------------------------------------------------------------------------------------ //




// // ---------------------------------------------------------------------------------------------- //
// class P7z : public Grid {  // ------------------------------------------------------------------- //
// public:
//  T7z ( size_t N_in, double a_in=-1.0, double x0_in=0.0, double b_in=1.0 ) : Grid(N_in,"P7z") {
//    set_params();
//    auto_grid();
//  }
//  void reset ( double a_in, double x0_in, double b_in ) {
//    a=a_in; x0=x0_in; b=b_in;
//    set_params(a,x0,b);
//    auto_grid();
//  }

// private:
// ////// parameters of the mesh:
//  double
//    zm=0,   // offset of the center points
//    m1=0.0, // slope at end points, 
//    m2=0.7, // slope in center, 
//    s=0.0;  // second derivative at end
// ////// coefficients in the polynomial:
//  double c7,c6,c5,c4,c3,c2,c1,c0;

//  double x_func ( size_t i ) const {
//    double z = (1-zm) + double(i)/double(N_points-1);
//    double xt =
//        c7*std::pow(z,7)
//      + c6*std::pow(z,6)
//      + c5*std::pow(z,5)
//      + c4*std::pow(z,4)
//      + c3*std::pow(z,3)
//      + c2*std::pow(z,2)
//      + c1*z
//      + c0;
//    return 0.5*(a+b) - 0.5*(a-b)*xt;
//  }
// ////// coefficients in polynomial
//  void set_params ( double a_in, double x0_in, double b_in ) {
//    zm = ( x0-a ) / ( b-a );

//  }
// };  // ------------------------------------------------------------------------------------------ //






// // ------------------------------------------------------------------------------------------------------------ //
// class T7z : public Grid { // ---------------------------------------------------------------------------------- //
// public:
//  T7z ( size_t N_in ) : Grid(N_in,"T7z") {}
//  T7z ( size_t N_in, double a_in, double b_in ) : Grid(N_in,"T7z") {
//    a=a_in; b=b_in;
//    set_params();
//    auto_grid();
//  }
//  void reset ( double a_in, double x0_in, double b_in ) {
//    a=a_in; x0=x0_in; b=b_in;
//    set_params();
//    auto_grid();
//  }
// private:
// ////// parameters of the mesh:
//  double
//    m1=0.0, // slope at ends
//    m2=0.7, // slope in center 
//    s=-1.0; // second derivative at ends
// ////// coefficients in the polynomial:
//  double c7,c5,c3,c1;
//  double x_func ( size_t i ) const {
//    double z = double(2*i)/double(N_points-1) - 1.0;
//    double xt =
//        c7*std::pow(z,7)
//      + c5*std::pow(z,5)
//      + c3*std::pow(z,3)
//      + c1*z;
//    return 0.5*(a+b) - 0.5*(a-b)*xt;
//  }
// ////// coefficients in polynomial
//  void set_params () {
//    c7 = 15.0/8.0  - m1*7.0/8.0  - m2     + s/8.0;
//    c5 = -21.0/4.0 + m1*9.0/4.0  + m2*3.0 - s/4.0;
//    c3 = 35.0/8.0  - m1*11.0/8.0 - m2*3.0 + s/8.0;
//    c1 = m2;
//  }
// }; // --------------------------------------------------------------------------------------------------------- //



#endif

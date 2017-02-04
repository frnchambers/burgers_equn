// Copyright 2016 frnchambers

#ifndef __GRID_GRID_HPP__
#define __GRID_GRID_HPP__

#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>


// ---------------------------------------------------------------------------------------------- //
//template <class vec_type>  // ------------------------------------------------------------------- //
class Grid {  // -- Base Class to define a grid ------------------------------------------------- //
 public:
  typedef std::vector<double> vec_type;

 protected:
  size_t N_points;
  vec_type xpts, wpts;
////// initialise the grid memory, doesn't initialise points xpts or wpts
  // Grid () {}
  Grid(size_t N_in, std::string name_in = "",
       double a_in = -1.0, double x0_in = 0.0, double b_in = 1.0 )
    : N_points(N_in), xpts(N_in), wpts(N_in), name(name_in),
      a(a_in), x0(x0_in), b(b_in)
  {}
////// to be defined and depend on the grid
  virtual double x_func(const size_t & i) const = 0;
////// initialise x_grid and weights from this grid
  void auto_grid() {
    for (size_t i = 0; i < N_points; ++i)
      xpts[i] = x_func(i);
    w_from_x();
  }
////// for grids that don't have a simple analytic w
  void w_from_x() {
    for (size_t i = 0; i < N_points; i++) {
      wpts[i] = 1.0;
      for (size_t k = 0; k < N_points; k++)
        if (k != i)
          wpts[i] *= xpts[i]-xpts[k];
      wpts[i] = 1.0/wpts[i];
    }
    scale_w();
  }
////// scale wpts so that largest value is 1
  void scale_w() {
    double w_big = 0.0;
    for ( const double &w : wpts )
      if ( std::abs(w) > w_big )
        w_big = std::abs(w);
    for ( double &w : wpts )
      w /= w_big;
  }

 public:  // ------------------------------------------------------------------------------------ //
    std::string name = "";  // name of the grid....
  double a = -1.0, x0 = 0.0, b = 1.0;  // end points and 'center' of the grid
////// return size of grid
  size_t size() const {
    return N_points;
  }
////// call element i of x or w grid
  double operator[] (const size_t & i) const {
    return xpts[i];
  }
  double w(const size_t & i) const {
    return wpts[i];
  }

}; // ------------------------------------------------------------------------------------------- //



// ---------------------------------------------------------------------------------------------- //
class Even : public Grid {  // -- Evenly space grid --------------------------------------------- //
 public:
  Even(size_t N_in,
        double a_in = -1.0, double x0_in = 0.0, double b_in = 1.0)
    : Grid(N_in, "Even", a_in, x0_in, b_in) {
    auto_grid();
  }

  void reset(double a_in, double x0_in, double b_in) {
    a = a_in; x0 = x0_in; b = b_in;
    auto_grid();
  }

 private:
////// default parameters for grid centered at 0, with boundaries: -1,1
  double x_func(const size_t & i) const {
    return a + ( b-a ) * static_cast<double>(i) / static_cast<double>(N_points-1);
  }
};  // ------------------------------------------------------------------------------------------ //



// ---------------------------------------------------------------------------------------------- //
class Sinh : public Grid {  // -- Shiny log - A.L. Watts 2016 ----------------------------------- //
 public:
  Sinh(size_t N_in,
        double a_in = -1.0, double x0_in = 0.0, double b_in = 1.0)
    : Grid(N_in, "Sinh", a_in, x0_in, b_in) {
    set_params(a, x0, b);
    auto_grid();
  }

  void reset(double a_in, double x0_in, double b_in) {
    a = a_in; x0 = x0_in; b = b_in;
    set_params(a, x0, b);
    auto_grid();
  }

 private:
////// default parameters for grid centered at 0, with boundaries: -1,1
  double d = -1.0, m = 2.0;
  //  formulae for x0, m, and d: x = x0 + sinh( d + m * i/(N-1) )
  //  for x 'element of' [a,b]:
  //                -> d = asinh( a-x0 )
  //                -> m = asinh( b-x0 ) - asinh( a-x0 )
  double x_func(const size_t & i) const {
    double z = d + m * static_cast<double>(i) / static_cast<double>(N_points-1);
    return x0 + std::sinh( z );
  }
////// set up parameters that characterise the grid
  void set_params(const double & a_in, const double & x0_in, const double & b_in) {
    d = std::asinh(a_in - x0_in);
    m = std::asinh(b_in - x0_in) - d;
  }
};  // ------------------------------------------------------------------------------------------ //




// ---------------------------------------------------------------------------------------------- //
class Cheb_1 : public Grid {  // -- gauss-lobato mesh ------------------------------------------- //
 public:
  Cheb_1(size_t N_in,
         double a_in = -1.0, double x0_in = 0.0, double b_in = 1.0)
      : Grid(N_in, "Cheb_1", a_in, 0.5*(a_in+b_in), b_in) {
    set_grid();
  }
  void reset(double a_in, double b_in) {
    a = a_in; b = b_in;
    set_grid();
  }

 private:
  double x_func(const size_t & i) const {
    double t = std::cos(M_PI * static_cast<double>(2*i + 1)
                        / static_cast<double>(2*N_points) );
    return 0.5*(a+b) - 0.5*(a-b)*t;
  }
  double w_func(const size_t & i) const {
    return ( i%2 == 1 ? -1.0 : 1.0 ) *
        std::sin(M_PI * static_cast<double>(2*i + 1)
                 / static_cast<double>(2*N_points) );
  }
  void set_grid() {
    for ( size_t i = 0; i < N_points; ++i ) {
      xpts[i] = x_func(i);
      wpts[i] = w_func(i);
    }
  }
};  // ------------------------------------------------------------------------------------------ //



// ---------------------------------------------------------------------------------------------- //
class Cheb_2 : public Grid {  // ---------------------------------------------------------------- //
 public:
  Cheb_2(size_t N_in,
         double a_in = -1.0, double x0_in = 0.0, double b_in = 1.0)
    : Grid(N_in, "Cheb_2", a_in, 0.5*(a_in+b_in), b_in) {
    a = a_in; b = b_in;
    set_grid();
  }
  void reset(double a_in, double x0_in, double b_in) {
    x0_in=x0_in;
    a = a_in; b = b_in;
    set_grid();
  }

 private:
  double x_func(const size_t & i) const {
    double t = std::cos(M_PI * static_cast<double>(i)
                        / (static_cast<double>(N_points-1) ) );
    return 0.5*(a+b) + 0.5*(a-b)*t;
  }
  double w_func(const size_t & i) const {
    return ( i == 0 || i == N_points-1 ? 0.5 : 1.0 )
        * (i%2 == 1 ? -1.0 : 1.0);
  }
  void set_grid() {
    for (size_t i = 0; i < N_points; ++i) {
      xpts[i] = x_func(i);
      wpts[i] = w_func(i);
    }
  }
};  // ------------------------------------------------------------------------------------------ //



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

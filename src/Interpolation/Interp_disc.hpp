#ifndef __INTERPOLATION_INTERP_DISC_HPP__
#define __INTERPOLATION_INTERP_DISC_HPP__

#include <iostream>
#include <vector>

template<typename vec_type, class base_interp>
class Discontinuity {

public:

  size_t
  // peak/trough about discontinuity
    ipeak=0, itrough=0,
    N_lhs=0, N_mid=0, N_rhs=0;
  base_interp
  // splines for each of the three
    f_lhs_interp, f_mid_interp, f_rhs_interp;
  vec_type
  // vectors for either side and center of discontinuity
    x_lhs, x_mid, x_rhs,
    f_lhs, f_mid, f_rhs;

  Discontinuity () {}

  Discontinuity ( const vec_type & xpts, const vec_type & fpts ) {
    set_points(xpts,fpts);
  }

private:
  size_t i_peak ( const vec_type & fpts ) const {
  // value of i for peak temperature
    size_t ip=0;
    while ( fpts[ip] < fpts[ip+1] )
      ip++;
    return ip;
  }
  size_t i_trough ( const vec_type & fpts ) const {
  // value of i to return trough
    size_t it=fpts.size()-1;
    while ( fpts[it] > fpts[it-1] )
      it--;
    return it;
  }

public:
  void set_points ( const vec_type & xpts, const vec_type & fpts ) {  // -------------------- //
    if ( xpts.size() != fpts.size() ) {
      std::cerr << "ERROR IN DISC SPLINE, mismatch in size" << std::endl;
      throw(0);
    }
    size_t N_points = xpts.size();

  ////// find locations of peaks/troughs and sizes of domains
    ipeak = i_peak(fpts);
    itrough = i_trough(fpts);
    N_lhs = ipeak + 1;
    N_rhs = N_points - itrough;
    N_mid = N_points + 2 - (N_lhs + N_rhs); // +2 from shared peak and trough

  ////// reset the size of different objects
    x_lhs.resize(N_lhs);
    x_mid.resize(N_mid);
    x_rhs.resize(N_rhs);

    f_lhs.resize(N_lhs);
    f_mid.resize(N_mid);
    f_rhs.resize(N_rhs);

  ////// put data in the correct domain
    for ( size_t i=0; i<N_lhs; ++i ) {
      x_lhs[i] = xpts[i];
      f_lhs[i] = fpts[i];
    }
    for ( size_t i=0; i<N_mid; ++i ) {
      x_mid[i] = xpts[ipeak+i];
      f_mid[i] = fpts[ipeak+i];
    }
    for ( size_t i=0; i<N_rhs; ++i ) {
      x_rhs[i] = xpts[itrough+i];
      f_rhs[i] = fpts[itrough+i];
    }

  ////// reset the interpolatino routines
    f_lhs_interp.set_points(x_lhs,f_lhs);
    f_mid_interp.set_points(x_mid,f_mid);
    f_rhs_interp.set_points(x_rhs,f_rhs);
  }

  double operator() ( const double & x ) const {
  ////// calculate f using interpolation, decide which domain
    double f;
    if ( x > x_mid[N_mid-1] ) {
      f    = f_rhs_interp(x);
    } else if ( x < x_mid[0] ) {
      f    = f_lhs_interp(x);
    } else if ( x <= x_mid[N_mid-1] && x >= x_mid[0] ) {
      f    = f_mid_interp(x);
    } else {
      std::cerr << "-> ERROR with interpolation in Back_mode, given x not in any of the possible ranges\n";
      throw(0);
    }
    return f;
  }

  void operator() ( const double & x,
                    double & f, double & dfdx ) const {
  ////// calculate f using interpolation, decide which domain
    if ( x > x_mid[N_mid-1] ) {
      f    = f_rhs_interp(x);
      dfdx = f_rhs_interp.deriv(x);
    } else if ( x < x_mid[0] ) {
      f    = f_lhs_interp(x);
      dfdx = f_lhs_interp.deriv(x);
    } else if ( x <= x_mid[N_mid-1] && x >= x_mid[0] ) {
      f    = f_mid_interp(x);
      dfdx = f_mid_interp.deriv(x);
    } else {
      std::cerr << "-> ERROR with interpolation in Back_mode, given x not in any of the possible ranges\n";
      throw(0);
    }
  }

  void plot_partition () {
    std::ofstream
      plt_out("data/plot_partition.gp"),
      lhs_out("data/lhs.dat"), dis_out("data/dis.dat"), rhs_out("data/rhs.dat");

    plt_out <<
      "set terminal qt size 500,450 persist\n" <<
      "set log x\n" << "set log y\n" <<
      "set xlabel 'x'\n" <<
      "set ylabel 'f'\n" <<
      "plot 'lhs.dat' w linespoints title 'lhs', 'dis.dat' w linespoints title 'dis', 'rhs.dat' w linespoints title 'rhs'\n";

    for ( size_t i=0; i<x_lhs.size(); ++i )
      lhs_out << x_lhs[i] << ' ' << f_lhs[i] << '\n';
    for ( size_t i=0; i<x_mid.size(); ++i )
      dis_out << x_mid[i] << ' ' << f_mid[i] << '\n';
    for ( size_t i=0; i<x_rhs.size(); ++i )
      rhs_out << x_rhs[i] << ' ' << f_rhs[i] << '\n';
  }

};

#endif

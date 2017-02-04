#include <iostream>
#include <iomanip>
#include <fstream>

#include <Grid/Grid.hpp>
#include "Function.h"

#include <Interpolation/Interp_GSL-spline.hpp>
#include <Interpolation/Interp_disc.hpp>

//#include <boost/numeric/ublas/vector.hpp>
//typedef boost::numeric::ublas::vector< double > vec;

#include <vector>
typedef std::vector<double> vec;

int main( int argc, char * argv[] )
{

  vec xpts = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
  vec fpts = {2.0,2.1,2.2,2.3,2.4,2.5, 2.2,1.9,1.6,1.3, 1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9};

  std::string pt_file = "data/points_disc.dat";
  std::ofstream points(pt_file);  
  points << std::scientific << std::setprecision(6);

  std::cout <<
    "# starting to populate x and f vectors, and printing points to file: " <<
    pt_file << '\n';

  for ( size_t i=0; i<xpts.size(); ++i )
    points << xpts[i] << ' ' << fpts[i] << '\n';


  //GSL_spline<vec,vec> interp(xpts,fpts);
  Discontinuity<vec,GSL_spline<vec> > interp(xpts,fpts);


  std::string ln_file = "data/line_disc.dat";
  std::ofstream lines(ln_file);
  points << std::scientific << std::setprecision(6);

  std::cout <<
    "# initialising interpolation object and printing interpolating line to file: " <<
    ln_file << '\n';

  double a=xpts[0], b=xpts[xpts.size()-1];

  size_t N_dense=101;
  for ( size_t i=0; i<N_dense; ++i ) {
    double t = double(2*i)/double(N_dense-1) - 1.0;
    double x = 0.5*(a+b) - 0.5*(a-b)*t;

    lines <<
      x << ' ' <<
      interp(x) << ' ' <<
  
      '\n';
  }


}

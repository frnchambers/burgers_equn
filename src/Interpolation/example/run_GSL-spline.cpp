#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>

#include <Grid/Grid.hpp>
#include <Interpolation/Interp_GSL-spline.hpp>
#include "Function.h"


typedef boost::numeric::ublas::vector< double > vec;


int main( int argc, char * argv[] )
{
  size_t N_data=11;
  if (argc==2)
    N_data=atoi(argv[1]);

  Cheb_2 xpts(N_data);

  std::string pt_file = "data/points_"+xpts.name+".dat";
  std::ofstream points(pt_file);  
  points << std::scientific << std::setprecision(6);

  std::cout <<
    "# starting to populate x and f vectors, and printing points to file: " <<
    pt_file << '\n';

  vec xpts_vec(N_data), fpts(N_data);
  for ( size_t i=0; i<N_data; ++i ) {
    xpts_vec[i] = xpts[i];
    fpts[i] = f(xpts[i]);
    points << xpts[i] << ' ' << fpts[i] << '\n';
  }

  // ------------------------------ begin interpolation utilities ------------------------------ //
  std::cout << "# initialising interpolation object and printing interpolating line to file: line_GSL-spline.dat\n";

  GSL_spline<vec,vec> interp(xpts_vec,fpts);

  std::string ln_file = "data/line_"+xpts.name+"_"+interp.interp_name()+".dat";
  std::ofstream lines(ln_file);
  points << std::scientific << std::setprecision(6);

  std::cout <<
    "# initialising interpolation object and printing interpolating line to file: " <<
    ln_file << '\n';

  size_t N_dense=(10*N_data)+1;
  for ( size_t i=0; i<N_dense; ++i ) {
    double x = double(2*i)/double(N_dense-1) - 1.0;

    lines <<
      x << ' ' <<
      interp(x) << ' ' <<
      f(x) << ' ' <<
      interp(x)-f(x) << ' ' <<

      '\n';
  }


}

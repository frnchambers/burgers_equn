// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <cmath>

// #include <map>
// #include <string>
// #include <vector>

// #include <boost/numeric/ublas/vector.hpp>
// #include <boost/numeric/ublas/matrix.hpp>

#include <Grid/Grid.hpp>


typedef std::vector<double> vec_type;


// int main( int argc, char * argv[] ) {

//   size_t N_data=11;
//   if ( argc==2 )
//     N_data=atoi(argv[1]);

//   ////// make the different grid
//   std::map<std::string, Grid< vec_type > *> grid_funcs;
//   grid_funcs["Cheb_1"]  = new Cheb_1<vec_type>(N_data);
//   grid_funcs["Cheb_2"]  = new Cheb_2<vec_type>(N_data);
//   grid_funcs["Shiny"]   = new Shiny<vec_type>(N_data);
//   grid_funcs["T7z"]     = new T7z<vec_type>(N_data);

//   for ( const auto &grid : grid_funcs ) {

//     ////// print out points
//     std::ofstream g_points("data/gridpoints_"+grid.first+".dat");
//     g_points <<
//       std::scientific <<
//       std::setprecision(6);

//     for( size_t i=0; i<grid.second->size(); i++ )
//       g_points <<
//         i << ' ' <<
//         grid.second->x(i) << ' ' <<
//         grid.second->w(i) << '\n';

//   }

// }



template<class grid_type>
void print_points ( size_t N_data, double a, double x0, double b ) {

  grid_type grid(N_data,a,x0,b);

////// print out points
  std::cout << "# Printing data to data/gridpoints_" << grid.name << ".dat\n";
  std::ofstream g_points("data/gridpoints_"+grid.name+".dat");
  g_points <<
    std::scientific <<
    std::setprecision(6);

  for( size_t i=0; i<grid.size(); i++ )
    g_points <<
      i << ' ' <<
      grid.x(i) << ' ' <<
      grid.w(i) << '\n';
}


int main ( int argc, char * argv[] ) {

  size_t N_data=11;
  if ( argc==2 )
    N_data=atoi(argv[1]);
  double a=-10.0, b=10.0, x0=0.0;
	

  print_points<Sinh>(N_data,a,x0,b);
  print_points<Cheb_1>(N_data,a,x0,b);
  print_points<Cheb_2>(N_data,a,x0,b);

}

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <Grid/Grid.hpp>
#include <Interpolation/Interp_lagrange.hpp>

#include "Function.h"


int main( int argc, char * argv[] ) {

	std::cout << std::scientific << std::setprecision(5);

	size_t N_data=11;
	if ( argc>1 )
		N_data=atoi(argv[1]);

////// make the grid
	//T7z<lagrange_interpolate::vec_type> grid(N_data);
	//Cheb_2<lagrange_interpolate::vec_type> grid(N_data);
	//Shiney<lagrange_interpolate::vec_type> grid(N_data);

	//lagrange_interpolate::vec_type x_data(N_data);
	//for( size_t i=0; i<N_data; i++ )
	//	x_data[i] = -std::cos( M_PI * double(i) / ( double(N_data-1) ) );
	//Arb<lagrange_interpolate::vec_type> grid(x_data);

	std::cout << "# starting to populate f vector from grid, and printing points to file: points_lagrange.dat\n";

////// interpolation object
	// lagrange_interpolate<Shiny> interp(N_data);
	// interp.grid.reset(-1,0,1);
	lagrange_interpolate<Cheb_1> interp(N_data);
	interp.grid.reset(-1,1);

	interp.update_derivs();



////// make the functions and print out points
	std::ofstream points("data/points_lagrange.dat");
	points << std::scientific << std::setprecision(6);
	lagrange_interpolate<Shiny>::vec_type fpts( interp.size() );
	for ( size_t i=0; i<interp.size(); i++ ) {
		fpts[i] = f(interp.grid[i]);

		points <<
			interp.grid[i] << ' ' <<
			fpts[i] << ' ' <<
			std::endl;
	}


// ------------------------------ begin interpolation utilities ------------------------------ //
	std::cout << "# initialising interpolation object and printing interpolating line to file: line_lagrange.dat\n";

////// set up grid for printing out interpolation
	size_t N_dense=(10*N_data)+1;
////// print out interpolation
	std::ofstream line("data/line_lagrange.dat");
	line << std::scientific << std::setprecision(6);

	for( size_t i=0; i<N_dense; i++ ) {
		double x = interp.grid.a + ( ( interp.grid.b-interp.grid.a ) * double(i) / double(N_dense-1) );

		line <<
			x << ' ' <<
			interp(fpts,x) << ' ' <<
			f(x) << ' ' <<
			interp(fpts,x)-f(x) << ' ' <<

			'\n';
	}


	std::cout << "# calculating derivatives and printing files: derivs{1,2}_lagrange.dat\n";

////// derivative points
	lagrange_interpolate<Shiny>::vec_type
		dfpts  = interp.derivs(fpts),
		d2fpts = interp.derivs2(fpts);

////// print out derivs
	std::ofstream
		derivs1("data/derivs1_lagrange.dat"),
		derivs2("data/derivs2_lagrange.dat");

	derivs1 << std::scientific << std::setprecision(6);
	derivs2 << std::scientific << std::setprecision(6);

	for( size_t i=0; i<interp.size(); i++ ) {

		derivs1 <<
			interp.grid[i] << ' ' <<

			dfpts[i] << ' ' <<
			df(interp.grid[i]) << ' ' <<
			dfpts[i]-df(interp.grid[i]) << ' ' <<

			std::endl;

		derivs2 <<
			interp.grid[i] << ' ' <<

			d2fpts[i] << ' ' <<
			d2f(interp.grid[i]) << ' ' <<
			d2fpts[i]-d2f(interp.grid[i]) << ' ' <<

			std::endl;
	}


}

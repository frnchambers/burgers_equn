#ifndef __DERIVATIVES_TYPES_HPP__
#define __DERIVATIVES_TYPES_HPP__


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <blaze/Math.h>

#include <array>

////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
namespace grid_types {  // ---------------------------------------------------------------------- //
  template <std::size_t N>
  using vector = std::array<double, N>;
  //template <std::size_t N>
  //using matrix = boost::numeric::ublas::matrix< double > mat_type;
}  // ------------------------------------------------------------------------------------------- //
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //
namespace algebra {  // ------------------------------------------------------------------- //
  // using vector = boost::numeric::ublas::vector< double >;
  // using matrix = boost::numeric::ublas::matrix< double >;

  using vector = blaze::DynamicVector<double>;
  using matrix = blaze::DynamicMatrix<double>;
}  // ------------------------------------------------------------------------------------------- //
////////////////////////////////////////////////////////////////////////////////////////////////////

// namespace boost { namespace numeric { namespace odeint {

// struct is_resizeable< algebra::vector >
// {
//     typedef boost::true_type type;
//     static const bool value = type::value;
// };

// } } }

#endif

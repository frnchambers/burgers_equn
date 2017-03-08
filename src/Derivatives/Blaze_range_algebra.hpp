#ifndef __DERIVATIVES_BLAZE_RANGE_ALGEBRA_HPP__
#define __DERIVATIVES_BLAZE_RANGE_ALGEBRA_HPP__

#include <blaze/Math.h>
#include <Derivatives/Types.hpp>


////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------------------------- //

#include <boost/range.hpp>
#include <boost/mpl/size_t.hpp>

#include <boost/numeric/odeint/algebra/detail/macros.hpp>
#include <boost/numeric/odeint/algebra/detail/for_each.hpp>
#include <boost/numeric/odeint/algebra/detail/norm_inf.hpp>
#include <boost/numeric/odeint/algebra/norm_result_type.hpp>

namespace ode = boost::numeric::odeint;


template< class State > struct vector_space_norm_inf;


struct blaze_range_algebra
{
    template< class S1 , class Op >
    static void for_each1( S1 &s1 , Op op )
    {
      ode::detail::for_each1( begin( s1 ) , end( s1 ) ,
                op );
    }

    template< class S1 , class S2 , class Op >
    static void for_each2( S1 &s1 , S2 &s2 , Op op )
    {
        ode::detail::for_each2( begin( s1 ) , end( s1 ) ,
                begin( s2 ) , op );
    }

    template< class S1 , class S2 , class S3 , class Op >
    static void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
    {
        ode::detail::for_each3( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class Op >
    static void for_each4( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , Op op )
    {
        ode::detail::for_each4( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class Op >
    static void for_each5( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , Op op )
    {
        ode::detail::for_each5( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 , class Op >
    static void for_each6( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , Op op )
    {
        ode::detail::for_each6( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class Op >
    static void for_each7( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , Op op )
    {
        ode::detail::for_each7( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , begin( s7 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class Op >
    static void for_each8( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , Op op )
    {
        ode::detail::for_each8( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , begin( s7 ) , begin( s8 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class S9 , class Op >
    static void for_each9( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , S9 &s9 , Op op )
    {
        ode::detail::for_each9( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , begin( s7 ) , begin( s8 ) , begin( s9 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class S9 , class S10 , class Op >
    static void for_each10( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , S9 &s9 , S10 &s10 , Op op )
    {
        ode::detail::for_each10( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , begin( s7 ) , begin( s8 ) , begin( s9 ) , begin( s10 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class S9 , class S10 , class S11 , class Op >
    static void for_each11( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , S9 &s9 , S10 &s10 , S11 &s11 , Op op )
    {
        ode::detail::for_each11( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , begin( s7 ) , begin( s8 ) , begin( s9 ) , begin( s10 ) , begin( s11 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class S9 , class S10 , class S11 , class S12 , class Op >
    static void for_each12( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , S9 &s9 , S10 &s10 , S11 &s11 , S12 &s12 , Op op )
    {
        ode::detail::for_each12( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , begin( s7 ) , begin( s8 ) , begin( s9 ) , begin( s10 ) , begin( s11 ) , begin( s12 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class S9 , class S10 , class S11 , class S12 , class S13 , class Op >
    static void for_each13( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , S9 &s9 , S10 &s10 , S11 &s11 , S12 &s12 , S13 &s13 , Op op )
    {
        ode::detail::for_each13( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , begin( s7 ) , begin( s8 ) , begin( s9 ) , begin( s10 ) , begin( s11 ) , begin( s12 ) , begin( s13 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class S9 , class S10 , class S11 , class S12 , class S13 , class S14 , class Op >
    static void for_each14( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , S9 &s9 , S10 &s10 , S11 &s11 , S12 &s12 , S13 &s13 , S14 &s14 , Op op )
    {
        ode::detail::for_each14( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , begin( s7 ) , begin( s8 ) , begin( s9 ) , begin( s10 ) , begin( s11 ) , begin( s12 ) , begin( s13 ) , begin( s14 ) , op );
    }

    template< class S1 , class S2 , class S3 , class S4 , class S5 , class S6 ,class S7 , class S8 , class S9 , class S10 , class S11 , class S12 , class S13 , class S14 , class S15 , class Op >
    static void for_each15( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , S5 &s5 , S6 &s6 , S7 &s7 , S8 &s8 , S9 &s9 , S10 &s10 , S11 &s11 , S12 &s12 , S13 &s13 , S14 &s14 , S15 &s15 , Op op )
    {
        ode::detail::for_each15( begin( s1 ) , end( s1 ) , begin( s2 ) , begin( s3 ) , begin( s4 ) , begin( s5 ) , begin( s6 ) , begin( s7 ) , begin( s8 ) , begin( s9 ) , begin( s10 ) , begin( s11 ) , begin( s12 ) , begin( s13 ) , begin( s14 ) , begin( s15 ) , op );
    }

  
    static double norm_inf( const blaze::DynamicVector<double> &s )
    {
      return ode::detail::norm_inf( std::begin( s ) , std::end( s ) ,
                                    static_cast< double >( 0 ) );
    }

};





// ---------------------------------------------------------------------------------------------- //
////////////////////////////////////////////////////////////////////////////////////////////////////


#endif

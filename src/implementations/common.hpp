#ifndef GOLDENROCEKEFELLER_FAST_ADDITIVE_COMPARISON_CONSTANTS_HPP
#define GOLDENROCEKEFELLER_FAST_ADDITIVE_COMPARISON_CONSTANTS_HPP

#include <cmath>

namespace goldenrockefeller{ namespace fast_additive_comparison{
    // template <typename T>
    // struct typed_constants{
    //     static const T pi;
    //     static const T tau;
    //     static const T inv_tau;
    // };

    template <typename T>
    struct typed_constants{
        static constexpr T pi = T(3.141592653589793238462643383279502884197L);
        static constexpr T tau = T(6.283185307179586476925286766559005768394L);
        static constexpr T inv_tau = T(0.1591549430918953357688837633725143620345L);
    };

    // template <typename T>
    // struct typed_constants{
    //     static T pi; 
    //     static T tau;
    //     static T inv_tau;
    // };

    // template <>
    // struct typed_constants<float> : typed_constants_numeric<float> {};

    // template <>
    // struct typed_constants<double> : typed_constants_numeric<double> {};

    // template <>
    // struct typed_constants<long double> : typed_constants_numeric<long double> {};

    template <typename T>
    inline T wrap_phase(const T& phase) { 
        using k = typed_constants<T>;
        /* Wrap phase between -pi, and pi */
        return phase - round(phase * k::inv_tau) * k::tau;
    }

    template <typename T>
    inline T wrap_phase_bounded(const T& phase) {
        
        using k = typed_constants<T>;
        /* Wrap phase between -pi, and pi, but raw phase must be between -pi and 3 * pi */
        return phase - T(phase > k::pi) * k::tau;
    }   

    template <typename T>
    inline T wrap_phase_offset(const T& phase) {
        
        using k = typed_constants<T>;
        /* Wrap phase between 0, and 2 * pi, but raw phase must be between -pi and 3 * pi */
        return phase - floor(phase * k::inv_tau) * k::tau;
    } 

    template <typename T>
    inline T cos(const T& x);

    using std::cos;

    template <typename sample_type, typename operand_type>
    inline void load(const sample_type* ptr, operand_type& operand);

    template <typename sample_type>
    inline void load(const sample_type* ptr, sample_type& operand) {
        operand = *ptr;
    }

    // template <typename sample_type, typename operand_type>
    // inline void store(sample_type* ptr, operand_type& operand);

    template <typename sample_type>
    inline void store(sample_type* ptr, const sample_type& operand) {
        *ptr = operand;
    }

}}

#endif
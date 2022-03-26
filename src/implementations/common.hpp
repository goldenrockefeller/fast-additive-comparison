#ifndef GOLDENROCEKEFELLER_FAST_ADDITIVE_COMPARISON_CONSTANTS_HPP
#define GOLDENROCEKEFELLER_FAST_ADDITIVE_COMPARISON_CONSTANTS_HPP

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
        return phase - k::tau * (phase > k::pi);
    }   

    template <typename T>
    inline T wrap_phase_offset(const T& phase) {
        
        using k = typed_constants<T>;
        /* Wrap phase between 0, and 2 * pi, but raw phase must be between -pi and 3 * pi */
        return phase - floor(phase * k::inv_tau) * k::tau;
    } 
}}

#endif
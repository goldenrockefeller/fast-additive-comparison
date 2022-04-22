#include <iostream>
#include <cstddef>
#include <vector>
#include <sstream>
#include <algorithm>
#include <numeric>
#include "nanobench.h"
#include "../implementations/phase-to-amplitude.hpp"
#include "../implementations/oscillator-bank.hpp"
#include "xsimd/xsimd.hpp"

namespace xs = xsimd;

namespace gfac = goldenrockefeller::fast_additive_comparison;

using std::cout;
using std::vector;
using std::size_t;
using std::ostringstream;
using std::iota;
using std::for_each;

using float_avx_t = xs::batch<float, xs::avx>;
using double_avx_t = xs::batch<double, xs::avx>;

using gfac::OscillatorBank;
using gfac::SimpleExactSineOscillator;
using gfac::SineOscillator;
using FloatCosCalc = gfac::ExactCosineCalculator<float, float_avx_t>;
using DoubleCosCalc = gfac::ExactCosineCalculator<double, double_avx_t>;


template <typename GeneratorT>
void do_inharmonic_bench(ankerl::nanobench::Bench* bench, char const* name, size_t chunk_size, size_t n_oscs) {
    using sample_type = typename GeneratorT::sample_type;

    GeneratorT gen(n_oscs);
    vector<sample_type> output(chunk_size);

    vector<sample_type> freqs(n_oscs);
    iota(freqs.begin(), freqs.end(), 0.);
    for_each(freqs.begin(), freqs.end(), [&] (sample_type& freq) {freq /= (2 * n_oscs);});

    bench->run(name, [&]() {
        for (size_t osc_id = 0; osc_id < n_oscs; ++osc_id) {
            gen.reset_osc(osc_id, freqs[osc_id], 1., 0.);
        }
        gen.progress_and_add(output.begin(), output.end());
    });
}

void do_all_inharmonic_benches(size_t chunk_size, size_t n_oscs) {
    ankerl::nanobench::Bench bench;

    ostringstream title_stream;
    title_stream << "All Inharmonic Bench. Chunck Size: " << chunk_size << "; Num of Oscs: " << n_oscs;
    bench.title(title_stream.str());

    do_inharmonic_bench<OscillatorBank<SimpleExactSineOscillator<double>>>(
        &bench, "Phase-to-Amplitude Simple Double", chunk_size, n_oscs
    );

    do_inharmonic_bench<OscillatorBank<SineOscillator<double, double_avx_t, 1,DoubleCosCalc>>>(
        &bench, "Phase-to-Amplitude Exact Double-AVX-1", chunk_size, n_oscs
    );

    do_inharmonic_bench<OscillatorBank<SineOscillator<double, double_avx_t, 2,DoubleCosCalc>>>(
        &bench, "Phase-to-Amplitude Exact Double-AVX-2", chunk_size, n_oscs
    );
}
 
int main() {
    cout << "Hello World!\n";  

    do_all_inharmonic_benches(50000, 1);
    do_all_inharmonic_benches(1024, 1);
    do_all_inharmonic_benches(1, 1);
  
    return 0;
}
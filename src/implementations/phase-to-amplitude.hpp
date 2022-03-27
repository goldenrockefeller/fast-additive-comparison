#ifndef GOLDENROCEKEFELLER_FAST_ADDITIVE_IMPLEMENTATIONS_PHASE_TO_AMPLITUDE_HPP
#define GOLDENROCEKEFELLER_FAST_ADDITIVE_IMPLEMENTATIONS_PHASE_TO_AMPLITUDE_HPP
#include <cstddef>
#include <vector>

#include <cstddef>
#include <vector>
#include <cmath>
#include <iterator>

#include "common.hpp"


namespace goldenrockefeller{ namespace fast_additive_comparison{
    template <typename sample_type>
    class SimpleExactSineOscillator {
        using size_t = std::size_t;
        using vector_t = typename std::vector<sample_type>;
        using vector_iterator_t = typename std::vector<sample_type>::iterator;

        sample_type freq;
        sample_type ampl;
        sample_type phase;

        public:
            SimpleExactSineOscillator(sample_type freq, sample_type ampl, sample_type phase) :
                freq(wrap_phase_offset(freq)),
                ampl(ampl),
                phase(phase)
            {}
            void progress_and_add(vector_iterator_t signal_begin_it, vector_iterator_t signal_end_it) {
                using k = typed_constants<float>;
                for (auto signal_it = signal_begin_it; signal_it < signal_end_it; ++signal_it) {
                    *signal_it += ampl * cos(this->phase);
                    this->phase += k::tau * this->freq;
                    this->phase = wrap_phase_bounded(this->phase);
                }
            }
        // end public
    };

    template <typename sample_type, typename operand_type>
    struct ExactCosineCalculator {
        static_assert(sizeof(operand_type) >= sizeof(sample_type), "The operand type size must be the same size as sample type");
        static_assert((sizeof(operand_type) % sizeof(sample_type)) == 0, "The operand type size must be a multiple of size as sample type");

        static constexpr size_t n_samples_per_operand = sizeof(operand_type) / sizeof(sample_type);

        static inline operand_type cos(operand_type& x) {
            operand_type y;
            sample_type* x_ptr = reinterpret_cast<sample_type*>(&x);
            sample_type* y_ptr = reinterpret_cast<sample_type*>(&y);

            for (size_t i = 0; i < n_samples_per_operand; i++) {
                y_ptr[i] = cos(x_ptr[i]);
            }
            return y;
        }
    };

    template <typename operand_type>
    struct IdentityCalculator {
        static inline operand_type cos(operand_type& x) {
            return x;
        }
    };

    template <typename sample_type, typename operand_type, std::size_t n_operands_per_block, typename CosineCalculatorT> 
    class SineOscillator{
        static_assert(sizeof(operand_type) >= sizeof(sample_type), "The operand type size must be the same size as sample type");
        static_assert((sizeof(operand_type) % sizeof(sample_type)) == 0, "The operand type size must be a multiple of size as sample type");
        static_assert(n_operands_per_block >= 1, "The operand block length must be positive");

        using size_t = std::size_t;
        using vector_t = typename std::vector<sample_type>;
        using vector_iterator_t = typename std::vector<sample_type>::iterator;

        static constexpr size_t n_samples_per_operand = sizeof(operand_type) / sizeof(sample_type);
        static constexpr size_t n_samples_per_block = n_operands_per_block * sizeof(operand_type) / sizeof(sample_type);

        sample_type freq;
        operand_type ampl_operand;

        operand_type delta_phase_per_block;

        vector_t osc_block;
        vector_t phase_block;
        vector_iterator_t osc_block_it;
        vector_iterator_t osc_block_safe_end_it;
        vector_iterator_t osc_block_safe_begin_it;

        static inline void progress_phase_operand(sample_type* phase_ptr, operand_type& delta_phase_per_block) {
            operand_type phase_operand;
            load(phase_ptr, phase_operand);
            phase_operand += delta_phase_per_block;
            phase_operand = wrap_phase_bounded(phase_operand);
            store(phase_ptr, phase_operand);
        }

        template <size_t N>
        struct ProgressPhaseLoop{
            static inline void progress(vector_t& phase_block, operand_type& delta_phase_per_block) {
                SineOscillator::progress_phase_operand(&phase_block[N * n_samples_per_operand], delta_phase_per_block);
                ProgressPhaseLoop<N-1>::progress(phase_block, delta_phase_per_block);
            }
        };

        template <>
        struct ProgressPhaseLoop<0> {
            static inline void progress(vector_t& phase_block, operand_type& delta_phase_per_block) {
                SineOscillator::progress_phase_operand(&phase_block[0], delta_phase_per_block);
            }
        };

        static inline void progress_osc_operand(sample_type* osc_ptr, sample_type* phase_ptr, operand_type& ampl_operand) {
            operand_type osc_operand;
            operand_type phase_operand;
            load(phase_ptr, phase_operand); 
            osc_operand  = ampl_operand * CosineCalculatorT::cos(phase_operand);
            store(osc_ptr, osc_operand);
        }

        template <size_t N>
        struct ProgressOscLoop{
            static inline void progress(vector_t& osc_block, vector_t& phase_block, operand_type& ampl_operand) {
                SineOscillator::progress_osc_operand(&osc_block[(N+1) * n_samples_per_operand], &phase_block[N * n_samples_per_operand], ampl_operand);
                ProgressOscLoop<N-1>::progress(osc_block, phase_block, ampl_operand);
            }
        };

        template <>
        struct ProgressOscLoop<0> {
            static inline void progress(vector_t& osc_block, vector_t& phase_block, operand_type& ampl_operand) {
                SineOscillator::progress_osc_operand(&osc_block[n_samples_per_operand], &phase_block[0], ampl_operand);
            }
        };
    
        public:  
            static vector_t new_phase_block(sample_type freq, sample_type phase) {
                vector_t phase_block(n_samples_per_block, 0.);

                auto delta_phase_per_sample = wrap_phase_offset(typed_constants<sample_type>::tau * freq);

                // Fill the first phase operand
                for (auto it = phase_block.begin(); it < phase_block.begin() + n_samples_per_operand; ++it){
                    *it = phase;
                    phase += delta_phase_per_sample;
                    phase = wrap_phase_bounded(phase);
                }

                operand_type delta_phase_per_operand(wrap_phase_offset(typed_constants<sample_type>::tau * freq * n_samples_per_operand));

                // Fill the first phase block
                operand_type prev_phase_operand;
                load(phase_block.data(), prev_phase_operand);
                for (size_t i = n_samples_per_operand; i < n_samples_per_block; i += n_samples_per_operand) {
                    operand_type phase_operand;
                    load(&phase_block[i], phase_operand);
                    phase_operand = wrap_phase_bounded(prev_phase_operand + delta_phase_per_operand);
                    prev_phase_operand = phase_operand;
                    store(&phase_block[i], phase_operand);
                }  

                return phase_block;
            }

            static vector_t new_osc_block(sample_type freq, sample_type ampl, sample_type phase) {
                auto phase_block = new_phase_block(freq, phase);
                operand_type ampl_operand(ampl);

                vector_t osc_block(phase_block.size() + n_samples_per_operand, 0.);              

                ProgressOscLoop<n_operands_per_block>::progress(osc_block, phase_block, ampl_operand);

                return osc_block;
            }

            SineOscillator(sample_type freq, sample_type ampl, sample_type phase) :
                freq(freq),
                ampl_operand(ampl),
                delta_phase_per_block(wrap_phase_offset(typed_constants<sample_type>::tau * freq * this->n_samples_per_block)),
                phase_block(SineOscillator::new_phase_block(freq, phase)),
                osc_block(SineOscillator::new_osc_block(freq, ampl, phase)),
                osc_block_safe_end_it(osc_block.begin() + n_samples_per_block),
                osc_block_safe_begin_it(osc_block.begin() + n_samples_per_operand),
                osc_block_it(osc_block.begin()+n_samples_per_operand)
            {}

            void progress_phase_block() {
                ProgressPhaseLoop<n_operands_per_block>::progress(this->phase_block, this->delta_phase_per_block);
            }

            void progress_osc_block(size_t sample_offset) {
                operand_type osc_operand_last;
                load(&(*this->osc_block_safe_end_it), osc_operand_last);
                store(this->osc_block.data(), osc_operand_last);

                this->progress_phase_block();
                ProgressOscLoop<n_operands_per_block>::progress(this->osc_block, this->phase_block, this->ampl_operand);
                                
                this->osc_block_it = this->osc_block.begin() + sample_offset;
            }

            void progress_and_add(vector_iterator_t signal_begin_it, vector_iterator_t signal_end_it)    {

                if (signal_end_it < signal_begin_it) {
                    return;
                }

                if (signal_end_it - signal_begin_it < this->n_samples_per_operand) { // it is not safe to vectorize
                    for (auto signal_it = signal_begin_it; signal_it < signal_end_it; ++signal_it) {
                        if (this->osc_block_it >  this->osc_block_safe_end_it) {
                            this->progress_osc_block(size_t(this->osc_block_it - this->osc_block_safe_end_it));
                        }

                        *signal_it += *this->osc_block_it;
                        ++this->osc_block_it;
                    }
                } 

                else { // it is safe to vectorize
                    auto signal_safe_end_it = signal_end_it - this->n_samples_per_operand;
                    operand_type signal_operand_last;
                    load(&(*signal_safe_end_it), signal_operand_last);

                    auto signal_it = signal_begin_it;
                    for (; signal_it < signal_safe_end_it; signal_it += this->n_samples_per_operand) {

                        if (this->osc_block_it >  this->osc_block_safe_end_it) {
                            this->progress_osc_block(size_t(this->osc_block_it - this->osc_block_safe_end_it));
                        }
                        
                        operand_type signal_operand;
                        operand_type osc_operand;
                        load(&(*signal_it), signal_operand);
                        load(&(*this->osc_block_it), osc_operand);

                        signal_operand += osc_operand;

                        store(&(*signal_it), signal_operand);

                        this->osc_block_it += this->n_samples_per_operand;
                    }

                    this->osc_block_it -= signal_it - signal_safe_end_it;

                    if (this->osc_block_it >  this->osc_block_safe_end_it) {
                        this->progress_osc_block(size_t(this->osc_block_it - this->osc_block_safe_end_it));
                    }

                    operand_type osc_operand;
                    load(&(*this->osc_block_it), osc_operand);

                    signal_operand_last = signal_operand_last + osc_operand;

                    store(&(*signal_safe_end_it), signal_operand_last);

                    this->osc_block_it += this->n_samples_per_operand;
                }
            }
        // public
    };
}}

#endif
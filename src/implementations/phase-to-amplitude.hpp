#ifndef GOLDENROCEKEFELLER_FAST_ADDITIVE_IMPLEMENTATIONS_PHASE_TO_AMPLITUDE_HPP
#define GOLDENROCEKEFELLER_FAST_ADDITIVE_IMPLEMENTATIONS_PHASE_TO_AMPLITUDE_HPP
#include <cstddef>
#include <vector>

#include <cstddef>
#include <vector>
#include <cmath>
#include <iterator>

#include "../common.hpp"


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

        static inline operand_type cos(operand_type x) {
            operand_type y;
            sample_type* x_ptr = reinterpret_cast<sample_type*>(&x);
            sample_type* y_ptr = reinterpret_cast<sample_type*>(&y);

            for (size_t i = 0; i < n_samples_per_operand; i++) {
                y_ptr[i] = std::cos(x_ptr[i]);
            }
            return y;
        }
    };

    template <typename sample_type, typename operand_type, size_t n_operands_per_block, typename CosineCalculatorT> 
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
        operand_type ampl;

        operand_type delta_phase_per_block;

        vector_t osc_block;
        vector_t phase_block;
        vector_iterator_t osc_block_it;
        vector_iterator_t osc_block_safe_end_it;
        vector_iterator_t osc_block_safe_begin_it;

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
            operand_type* phase_operands = static_cast<operand_type*>(phase_block.data());
            for (size_t i = 1; i < n_operands_per_block; ++i) {
                phase_operands[i] = wrap_phase_bounded(phase_operands[i-1] + delta_phase_per_operand);
            }  

            return phase_block;
        }

        static vector_t new_osc_block(sample_type freq, sample_type ampl, sample_type phase) {
            auto phase_block = new_phase_block(freq, phase);
            operand_type ampl_operand(ampl);

            vector_t osc_block(phase_block.size() + n_samples_per_operand, 0.);

            auto osc_block_safe_begin_it = osc_block.begin() + n_samples_per_operand;
            
            operand_type* osc_operands = static_cast<operand_type*>(osc_block.data());
            operand_type* phase_operands = static_cast<operand_type*>(phase_block.data());

            for (size_t i = 0; i < n_operands_per_block; ++i) {
                osc_operands[i + 1] = ampl_operand * CosineCalculatorT::cos(phase_operands[i]);
            } 

            return osc_block;
        }
    
        public:  
            SineOscillator(sample_type freq, sample_type ampl, sample_type phase) :
                freq(freq),
                ampl(ampl),
                delta_phase_per_block(wrap_phase_offset(typed_constants<sample_type>::tau * freq * this->n_samples_per_block)),
                phase_block(SineOscillator::new_phase_block(freq, phase)),
                osc_block(SineOscillator::new_osc_block(freq, ampl, phase)),
                osc_block_safe_end_it(osc_block.begin() + n_samples_per_block),
                osc_block_safe_begin_it(osc_block.begin() + n_samples_per_operand),
                osc_block_it(osc_block.begin()+n_samples_per_operand)
            {}

            void progress_phase_block() {
                operand_type* phase_operands = static_cast<operand_type*>(this->phase_block.data());

                for (size_t i = 0; i < n_operands_per_block; ++i) {
                    phase_operands[i] += this->delta_phase_per_block;
                    phase_operands[i] = wrap_phase_bounded(phase_operands[i]);
                } 
            }

            void progress_osc_block(size_t sample_offset) {
                operand_type* osc_operands_last = static_cast<operand_type*>(&(*this->osc_block_safe_end_it));
                operand_type* osc_operands_first = static_cast<operand_type*>(this->osc_block.data());
                *osc_operands_first = *osc_operands_last;

                this->progress_phase_block();

                operand_type* osc_operands = static_cast<operand_type*>(osc_block.data());
                operand_type* phase_operands = static_cast<operand_type*>(phase_block.data());

                for (size_t i = 0; i < n_operands_per_block; ++i) {
                    osc_operands[i + 1] = this->ampl * CosineCalculatorT::cos(phase_operands[i]);
                } 
                
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
                    operand_type signal_operands_last = static_cast<operand_type>(*signal_safe_end_it);

                    auto signal_it = signal_begin_it;
                    for (; signal_it < signal_safe_end_it; signal_it += this->n_samples_per_operand) {

                        if (this->osc_block_it >  this->osc_block_safe_end_it) {
                            this->progress_osc_block(size_t(this->osc_block_it - this->osc_block_safe_end_it));
                        }
                        
                        operand_type* signal_operands = static_cast<operand_type*>(&(*signal_it));
                        operand_type* osc_operands = static_cast<operand_type*>(&(*this->osc_block_it));

                        *signal_operands += *osc_operands;
                        
                        this->osc_block_it += this->n_samples_per_operand;
                    }


                    if (this->osc_block_it >  this->osc_block_safe_end_it) {
                        this->progress_osc_block(size_t(this->osc_block_it - this->osc_block_safe_end_it));
                    }

                    operand_type* osc_operands = static_cast<operand_type*>(&(*this->osc_block_it));
                    operand_type* signal_operands = static_cast<operand_type*>(&(*signal_safe_end_it));
                    *signal_operands = signal_operands_last + *osc_operands;

                    this->osc_block_it += this->n_samples_per_operand - (signal_it - signal_safe_end_it);
                }
            }
        // public
    };
}}

#endif
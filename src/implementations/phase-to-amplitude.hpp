#ifndef GOLDENROCEKEFELLER_FAST_ADDITIVE_IMPLEMENTATIONS_PHASE_TO_AMPLITUDE_HPP
#define GOLDENROCEKEFELLER_FAST_ADDITIVE_IMPLEMENTATIONS_PHASE_TO_AMPLITUDE_HPP
#include <cstddef>
#include <vector>
#include <cmath>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <array>

#include "common.hpp"


namespace goldenrockefeller{ namespace fast_additive_comparison{
    template <typename sample_type>
    class SimpleExactSineOscillator {
        using size_t = std::size_t;
        using vector_type = typename std::vector<sample_type>;

        sample_type freq;
        sample_type ampl;
        sample_type phase;

        public:
            typedef sample_type sample_type;

            SimpleExactSineOscillator() : SimpleExactSineOscillator(sample_type(0), sample_type(0), sample_type(0)) {}

            SimpleExactSineOscillator(sample_type freq, sample_type ampl, sample_type phase) :
                freq(wrap_phase_offset(freq)),
                ampl(ampl),
                phase(phase)
            {}

            void reset(sample_type freq, sample_type ampl, sample_type phase) {
                this->freq = freq;
                this->ampl = ampl;
                this->phase = phase;
            }

            template<typename iterator_type>
            void progress_and_add(iterator_type signal_begin_it, iterator_type signal_end_it) {
                
                for (auto signal_it = signal_begin_it; signal_it < signal_end_it; ++signal_it) {
                    *signal_it += ampl * cos(this->phase);
                    this->phase += tau<sample_type>() * this->freq;
                    this->phase = wrap_phase_bounded(this->phase);
                }
            }
        // end public
    };

    template <typename sample_type>
    struct ExactCosineCalculator {

        template <typename operand_type>
        static inline operand_type cos(const operand_type& x) {
            static_assert(sizeof(operand_type) >= sizeof(sample_type), "The operand type size must be the same size as sample type");
            static_assert((sizeof(operand_type) % sizeof(sample_type)) == 0, "The operand type size must be a multiple of size as sample type");

            static constexpr size_t n_samples_per_operand = sizeof(operand_type) / sizeof(sample_type);

            operand_type y;
            const sample_type* x_ptr = reinterpret_cast<const sample_type*>(&x);
            sample_type* y_ptr = reinterpret_cast<sample_type*>(&y);

            for (size_t i = 0; i < n_samples_per_operand; i++) {
                y_ptr[i] = goldenrockefeller::fast_additive_comparison::cos(sample_type(x_ptr[i]));
            }
            return y;
        }
    };

    
    struct ApproxCos14Calculator {
        template <typename operand_type>
        static inline operand_type cos(const operand_type& x) {
            return approx_cos_deg_14(x);
        }
    };

    struct ApproxCos14CalculatorPre {
        template <typename operand_type>
        static inline operand_type cos(const operand_type& x) {
            return approx_cos_deg_14_pre(x);
        }
    };

    struct ApproxCos10Calculator {
        template <typename operand_type>
        static inline operand_type cos(const operand_type& x) {
            return approx_cos_deg_10(x);
        }
    };

    
    struct IdentityCalculator {
        template <typename operand_type>
        static inline operand_type cos(const operand_type& x) {
            return x;
        }
    };

    template <typename sample_type, typename operand_type, std::size_t n_operands_per_block, typename CosineCalculatorT> 
    class SineOscillator{
        static_assert(sizeof(operand_type) >= sizeof(sample_type), "The operand type size must be the same size as sample type");
        static_assert((sizeof(operand_type) % sizeof(sample_type)) == 0, "The operand type size must be a multiple of size as sample type");
        static_assert(n_operands_per_block >= 1, "The operand block length must be positive");

        using size_t = std::size_t;
        using vector_type = typename std::vector<sample_type>;
        using vector_iterator_type = typename std::vector<sample_type>::iterator;

        static constexpr size_t n_samples_per_operand = sizeof(operand_type) / sizeof(sample_type);
        static constexpr size_t n_samples_per_block = n_operands_per_block * sizeof(operand_type) / sizeof(sample_type);

        sample_type freq;
        operand_type ampl_operand;

        operand_type delta_phase_per_block;

        vector_type osc_block;
        vector_type phase_block;
        vector_iterator_type osc_block_it;
        vector_iterator_type osc_block_safe_end_it;
        vector_iterator_type osc_block_safe_begin_it;

        static inline void progress_phase_operand(sample_type* phase_ptr, const operand_type& delta_phase_per_block) {
            operand_type phase_operand;
            load(phase_ptr, phase_operand);
            phase_operand += delta_phase_per_block;
            phase_operand = wrap_phase_bounded(phase_operand);
            store(phase_ptr, phase_operand);
        }

        template <size_t N>
        struct ProgressPhaseLoop{
            static inline void progress(vector_type& phase_block, const operand_type& delta_phase_per_block) {
                SineOscillator::progress_phase_operand(&phase_block[(N-1) * n_samples_per_operand], delta_phase_per_block);
                ProgressPhaseLoop<N-1>::progress(phase_block, delta_phase_per_block);
            }
        };

        template <>
        struct ProgressPhaseLoop<1> {
            static inline void progress(vector_type& phase_block, const operand_type& delta_phase_per_block) {
                SineOscillator::progress_phase_operand(&phase_block[0], delta_phase_per_block);
            }
        };

        static inline void update_osc_operand(sample_type* osc_ptr, const sample_type* phase_ptr, const operand_type& ampl_operand) {
            operand_type osc_operand;
            operand_type phase_operand;
            load(phase_ptr, phase_operand); 
            osc_operand  = ampl_operand * CosineCalculatorT::cos(phase_operand);
            store(osc_ptr, osc_operand);
        }

        template <size_t N>
        struct UpdateOscLoop{
            static inline void update(vector_type& osc_block, const vector_type& phase_block, const operand_type& ampl_operand) {
                SineOscillator::update_osc_operand(&osc_block[N * n_samples_per_operand], &phase_block[(N-1) * n_samples_per_operand], ampl_operand);
                UpdateOscLoop<N-1>::update(osc_block, phase_block, ampl_operand);
            }
        };

        template <>
        struct UpdateOscLoop<1> {
            static inline void update(vector_type& osc_block, const vector_type& phase_block, const operand_type& ampl_operand) {
                SineOscillator::update_osc_operand(&osc_block[n_samples_per_operand], &phase_block[0], ampl_operand);
            }
        };
    
        public:
            typedef sample_type sample_type;

            static void init_phase_block(vector_type& phase_block, sample_type freq, sample_type phase) {
                if (phase_block.size() !=  n_samples_per_block) {
                    std::ostringstream msg;
                    msg << "The phase block size "
                        << "(phase_block.size() = " << phase_block.size() << ") "
                        << "must be equal to the number of samples per block "
                        << "(n_samples_per_block = " << n_samples_per_block << ") ";
                    throw std::invalid_argument(msg.str());
                }

                auto delta_phase_per_sample = wrap_phase_offset(tau<sample_type>() * freq);

                // Fill the first phase operand
                for (auto it = phase_block.begin(); it < phase_block.begin() + n_samples_per_operand; ++it){
                    *it = phase;
                    phase += delta_phase_per_sample;
                    phase = wrap_phase_bounded(phase);
                }

                operand_type delta_phase_per_operand(wrap_phase_offset(tau<sample_type>() * freq * n_samples_per_operand));

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

            }
            
            // static void update_osc_block(vector_type& osc_block, const vector_type& phase_block, sample_type ampl){
            //     if (phase_block.size() !=  n_samples_per_block) {
            //         std::ostringstream msg;
            //         msg << "The phase block size "
            //             << "(phase_block.size() = " << phase_block.size() << ") "
            //             << "must be equal to the number of samples per block "
            //             << "(n_samples_per_block = " << n_samples_per_block << ") ";
            //         throw std::invalid_argument(msg.str());
            //     }

            //     if (osc_block.size() !=  n_samples_per_block + n_samples_per_operand) {
            //         std::ostringstream msg;
            //         msg << "The oscillator block size "
            //             << "(osc_block.size() = " << osc_block.size() << ") "
            //             << "must be equal to the number of samples per block "
            //             << "(n_samples_per_block = " << n_samples_per_block << ") "
            //             << "plus the offset, which is the number of samples per operand"
            //             << "(n_samples_per_operand = " << n_samples_per_operand << ") ";
            //         throw std::invalid_argument(msg.str());
            //     }
            // }

            static vector_type new_phase_block(sample_type freq, sample_type phase) {

                vector_type phase_block(n_samples_per_block, 0.);                
                SineOscillator::init_phase_block(phase_block, freq, phase);
                return phase_block;
            }

            static vector_type new_osc_block(sample_type freq, sample_type ampl, sample_type phase) {
                auto phase_block = new_phase_block(freq, phase);
                vector_type osc_block(phase_block.size() + n_samples_per_operand, 0.);      

                operand_type ampl_operand(ampl);

                UpdateOscLoop<n_operands_per_block>::update(osc_block, phase_block, ampl_operand);

                return osc_block;
            }

            SineOscillator() : SineOscillator(sample_type(0), sample_type(0), sample_type(0)) {}

            SineOscillator(sample_type freq, sample_type ampl, sample_type phase) :
                freq(freq),
                ampl_operand(ampl),
                delta_phase_per_block(wrap_phase_offset(tau<sample_type>() * freq * this->n_samples_per_block)),
                phase_block(SineOscillator::new_phase_block(freq, phase)),
                osc_block(SineOscillator::new_osc_block(freq, ampl, phase)),
                osc_block_safe_end_it(osc_block.begin() + n_samples_per_block),
                osc_block_safe_begin_it(osc_block.begin() + n_samples_per_operand),
                osc_block_it(osc_block.begin()+n_samples_per_operand)
            {}

            void reset(sample_type freq, sample_type ampl, sample_type phase) {
                this->freq = freq;
                this->ampl_operand = operand_type(ampl);
                this->delta_phase_per_block = operand_type(wrap_phase_offset(tau<sample_type>() * freq * this->n_samples_per_block));
                SineOscillator::init_phase_block(this->phase_block, freq, phase);
                UpdateOscLoop<n_operands_per_block>::update(this->osc_block, this->phase_block, this->ampl_operand);
                this->osc_block_safe_end_it = this->osc_block.begin() + this->n_samples_per_block;
                this->osc_block_safe_begin_it = this->osc_block.begin() + this->n_samples_per_operand;
                this->osc_block_it = this->osc_block.begin() + this->n_samples_per_operand;
            }

            void progress_phase_block() {
                ProgressPhaseLoop<n_operands_per_block>::progress(this->phase_block, this->delta_phase_per_block);
            }

            void update_osc_block(size_t sample_offset) {
                operand_type osc_operand_last;
                
                load(&(*this->osc_block_safe_end_it), osc_operand_last);
                store(this->osc_block.data(), osc_operand_last);

                this->progress_phase_block();
                UpdateOscLoop<n_operands_per_block>::update(this->osc_block, this->phase_block, this->ampl_operand);
                                
                this->osc_block_it = this->osc_block.begin() + sample_offset;
            }

            template<typename iterator_type>
            void progress_and_add(iterator_type signal_begin_it, iterator_type signal_end_it)    {
                
                if (signal_end_it < signal_begin_it) {
                    return;
                }

                if  (signal_end_it - signal_begin_it < this->n_samples_per_operand) { // it is not safe to vectorize
                    for (auto signal_it = signal_begin_it; signal_it < signal_end_it; ++signal_it) {
                        if (this->osc_block_it >  this->osc_block_safe_end_it) {
                            this->update_osc_block(size_t(this->osc_block_it - this->osc_block_safe_end_it));
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
                            this->update_osc_block(size_t(this->osc_block_it - this->osc_block_safe_end_it));
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
                        this->update_osc_block(size_t(this->osc_block_it - this->osc_block_safe_end_it));
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
/*
    template <typename sample_type, typename operand_type, std::size_t n_operands_per_block> 
    struct ApproxSineOscillator : public SineOscillator {
        std::array<operand_type, 8> coefs;

         SineOscillator() : SineOscillator(sample_type(0), sample_type(0), sample_type(0)) {}

        SineOscillator(sample_type freq, sample_type ampl, sample_type phase) :
            freq(freq),
            ampl_operand(ampl),
            delta_phase_per_block(wrap_phase_offset(tau<sample_type>() * freq * this->n_samples_per_block)),
            phase_block(SineOscillator::new_phase_block(freq, phase)),
            osc_block(SineOscillator::new_osc_block(freq, ampl, phase)),
            osc_block_safe_end_it(osc_block.begin() + n_samples_per_block),
            osc_block_safe_begin_it(osc_block.begin() + n_samples_per_operand),
            osc_block_it(osc_block.begin()+n_samples_per_operand)
        {}

        ApproxSineOscillator() : ApproxSineOscillator(sample_type(0), sample_type(0), sample_type(0)) {}

        ApproxSineOscillator(sample_type freq, sample_type ampl, sample_type phase) : SineOscillator(sample_type freq, sample_type ampl, sample_type phase) {
            this-coefs = {
                coefs[0] = operand_type(0x1.ffffffff470fdp-1);
                coefs[1] = operand_type(-0x1.ffffffec1c40dp-2);
                coefs[2] = operand_type(0x1.555553f050eb2p-5);
                coefs[3] = operand_type(-0x1.6c169b776ec06p-10);
                coefs[4] = operand_type(0x1.a0160ea01af9bp-16);
                coefs[5] = operand_type(-0x1.27abf550a036ap-22);
                coefs[6] = operand_type(0x1.1b5c0b8055789p-29);
                coefs[7] = operand_type(-0x1.577f9d3aa99cep-37); 
            }
        }

        static inline void update_osc_operand(sample_type* osc_ptr, const sample_type* phase_ptr, const operand_type& ampl_operand) {
            operand_type osc_operand;
            operand_type phase_operand;
            load(phase_ptr, phase_operand); 
            osc_operand  = ampl_operand * approx_x2_deg_14(this->coefs, phase_operand);
            store(osc_ptr, osc_operand);
        }
    };
*/
}}



#endif
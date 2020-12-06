#include <misc/Exception.h>
#include <psdr/core/sampler.h>

namespace psdr
{

template <typename UInt32>
uint64_array_t<UInt32> sample_tea_64(UInt32 v0, UInt32 v1, int rounds = 4) {
    UIntC sum = 0;

    for ( int i = 0; i < rounds; ++i ) {
        sum += 0x9e3779b9;
        v0 += (sl<4>(v1) + 0xa341316c) ^ (v1 + sum) ^ (sr<5>(v1) + 0xc8013ea4);
        v1 += (sl<4>(v0) + 0xad90777d) ^ (v0 + sum) ^ (sr<5>(v0) + 0x7e95761e);
    }

    return uint64_array_t<UInt32>(v0) + sl<32>(uint64_array_t<UInt32>(v1));
}


// std::shared_ptr<Sampler> Sampler::clone() {
//     std::shared_ptr<Sampler> sampler = std::make_shared<Sampler>();
//     sampler->m_sample_count = m_sample_count;
//     sampler->m_base_seed = m_base_seed;
//     return sampler;
// }


void Sampler::seed(UInt64C seed_value) {
    if ( !m_rng )
        m_rng = std::make_unique<PCG32>();

    seed_value += m_base_seed;

    UInt64C idx = arange<UInt64C>(seed_value.size());
    m_rng->seed(sample_tea_64(seed_value, idx),
                sample_tea_64(idx, seed_value));

    m_sample_count = static_cast<int64_t>(slices(seed_value));
}


template <bool ad>
Float<ad> Sampler::next_1d() {
    if ( m_rng == nullptr )
        throw Exception("Sampler::seed() must be invoked before using this sampler!");
    else {
        if constexpr (ad) {
            return FloatD(m_rng->next_float32());
        } else {
            return m_rng->next_float32();
        }
    }
}


// Explicit instanciations
template FloatC Sampler::next_1d<false>();
template FloatD Sampler::next_1d<true>();

} // namespace psdr

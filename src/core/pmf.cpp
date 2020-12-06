#include <misc/Exception.h>
#include <psdr/core/pmf.h>

namespace psdr
{

void DiscreteDistribution::init(const FloatC &pmf) {
    m_size = static_cast<int>(slices(pmf));
    m_sum = hsum(pmf);
    m_pmf = pmf;
    m_cmf = psum(m_pmf);
    m_pmf_normalized = pmf/m_sum;
    //std::cout << m_cmf << std::endl;
}


std::pair<IntC, FloatC> DiscreteDistribution::sample(const FloatC &_samples) const {
    if ( unlikely(m_size == 1) ) {
        return { zero<IntC>(), full<FloatC>(1.f) };
    }
    FloatC samples = _samples*m_sum;
    IntC idx = binary_search(
        0, m_size - 1, [&](IntC i) { return gather<FloatC>(m_cmf, i) < samples; }
    );
    return { idx, gather<FloatC>(m_pmf, idx)/m_sum };
}


template <bool ad>
std::pair<IntC, FloatC> DiscreteDistribution::sample_reuse(Float<ad> &samples) const {
    if ( unlikely(m_size == 1) ) {
        return { zero<IntC>(), full<FloatC>(1.f) };
    }
    samples *= m_sum;
    IntC idx;
    if constexpr ( ad ) {
        idx = binary_search(
            0, m_size - 1, [&](IntC i) { return gather<FloatC>(m_cmf, i) < detach(samples); }
        );
    } else {
        idx = binary_search(
            0, m_size - 1, [&](IntC i) { return gather<FloatC>(m_cmf, i) < samples; }
        );
    }
    samples -= gather<FloatC>(m_cmf, idx - 1, idx > 0);
    FloatC pmf = gather<FloatC>(m_pmf, idx);
    masked(samples, pmf > 0.f) /= pmf;
    samples = clamp(samples, 0.f, 1.f);
    return { idx, pmf/m_sum };
}

// Explicit instantiations
template std::pair<IntC, FloatC> DiscreteDistribution::sample_reuse<false>(FloatC&) const;
template std::pair<IntC, FloatC> DiscreteDistribution::sample_reuse<true >(FloatD&) const;

} // namespace psdr

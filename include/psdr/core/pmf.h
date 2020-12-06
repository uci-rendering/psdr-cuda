#pragma once

#include <psdr/psdr.h>

namespace psdr
{

struct DiscreteDistribution {
    DiscreteDistribution() = default;

    void init(const FloatC &pmf);

    std::pair<IntC, FloatC> sample(const FloatC &samples) const;

    template <bool ad>
    std::pair<IntC, FloatC> sample_reuse(Float<ad> &samples) const;

    const FloatC &pmf() const { return m_pmf_normalized; }

    int     m_size;
    FloatC  m_sum;

protected:
    FloatC  m_pmf, m_pmf_normalized, m_cmf;
};

} // namespace psdr

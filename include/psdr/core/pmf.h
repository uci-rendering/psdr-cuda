#pragma once

#include <psdr/psdr.h>
#include <psdr/core/records.h>
#include <misc/Exception.h>

namespace psdr
{

struct DiscreteDistribution {
    DiscreteDistribution() = default;

    void init(const FloatC &pmf);

    std::pair<IntC, FloatC> sample(const FloatC &samples) const;

    template <bool ad>
    std::pair<IntC, FloatC> sample_reuse(Float<ad> &samples) const;

    const FloatC &pmf() const { return m_pmf_normalized; }
    const FloatC &cmf() const { return m_cmf_normalized; }

    int     m_size;
    FloatC  m_sum;

// protected:
    FloatC  m_pmf, m_pmf_normalized, m_cmf, m_cmf_normalized;
};

PSDR_CLASS_DECL_BEGIN(DDistribution,, Object)
public:
    virtual ~DDistribution() override {}

    virtual void init(const FloatC &pmf) = 0;
    virtual DiscreteRecordC sample_reuse(FloatC &samples) const = 0;
    virtual DiscreteRecordC sample_reuse(FloatD &samples) const = 0;
    virtual const FloatC &pmf() const = 0;
    virtual const FloatC &cmf() const = 0;

    ENOKI_PINNED_OPERATOR_NEW(FloatD)

PSDR_CLASS_DECL_END(DDistribution)


PSDR_CLASS_DECL_BEGIN(_DiscreteDistribution, final, DDistribution)
public:
    void init(const FloatC &pmf) override {
        m_size = static_cast<int>(slices(pmf));
        m_sum = hsum(pmf);
        m_pmf = pmf;
        m_cmf = psum(m_pmf);
        m_pmf_normalized = pmf/m_sum;
        m_cmf_normalized = m_cmf/m_sum;
    }

    DiscreteRecordC sample_reuse(FloatC &_samples) const override {
        FloatC samples(_samples);
        DiscreteRecordC result;
        if ( unlikely(m_size == 1) ) {
            PSDR_ASSERT(0);
            result.pdf = full<FloatC>(1.f);
            result.idx = zero<IntC>();
            return result;
        }
        samples *= m_sum;
        IntC idx;

        idx = binary_search(
            0, m_size - 1, [&](IntC i) { return gather<FloatC>(m_cmf, i) < samples; }
        );
        samples -= gather<FloatC>(m_cmf, idx - 1, idx > 0);
        FloatC pmf = gather<FloatC>(m_pmf, idx);
        masked(samples, pmf > 0.f) /= pmf;
        samples = clamp(samples, 0.f, 1.f);

        result.rnd = samples;
        result.pdf = pmf/m_sum;
        result.idx = idx;
        return result;
    }

    DiscreteRecordC sample_reuse(FloatD &samples) const override {

        DiscreteRecordC result;
        if ( unlikely(m_size == 1) ) {
            PSDR_ASSERT(0);
            result.pdf = full<FloatC>(1.f);
            result.idx = zero<IntC>();
            return result;
        }
        samples *= m_sum;
        IntC idx;
        idx = binary_search(
            0, m_size - 1, [&](IntC i) { return gather<FloatC>(m_cmf, i) < detach(samples); }
        );
        samples -= gather<FloatC>(m_cmf, idx - 1, idx > 0);
        FloatC pmf = gather<FloatC>(m_pmf, idx);
        masked(samples, pmf > 0.f) /= pmf;
        samples = clamp(samples, 0.f, 1.f);

        result.rnd = detach(samples);
        result.pdf = pmf/m_sum;
        result.idx = idx;
        return result;
    }


    const FloatC &pmf() const override { return m_pmf_normalized; }
    const FloatC &cmf() const override { return m_cmf_normalized; }

    int     m_size;
    FloatC  m_sum;
    FloatC  m_pmf, m_pmf_normalized, m_cmf, m_cmf_normalized;
    ENOKI_PINNED_OPERATOR_NEW(FloatD)

PSDR_CLASS_DECL_END(_DiscreteDistribution)

} // namespace psdr

ENOKI_CALL_SUPPORT_BEGIN(psdr::DDistribution)
    ENOKI_CALL_SUPPORT_METHOD(init)
    ENOKI_CALL_SUPPORT_METHOD(sample_reuse)
    ENOKI_CALL_SUPPORT_METHOD(pmf)
    ENOKI_CALL_SUPPORT_METHOD(cmf)
ENOKI_CALL_SUPPORT_END(psdr::DDistribution)


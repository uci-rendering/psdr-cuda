#pragma once

#include <psdr/psdr.h>
#include "pmf.h"
#include <psdr/core/bitmap.h>

namespace psdr
{

template <int ndim>
struct HyperCubeDistribution {
    static_assert(ndim > 1);

    void set_resolution(const Array<int, ndim> &reso);
    void set_mass(const FloatC &pmf);

    FloatC sample_reuse(Vectorf<ndim, false> &samples) const;
    FloatC pdf(const Vectorf<ndim, false> &p) const;

    bool                    m_ready = false;
    Array<int, ndim>        m_resolution = 0;

    DiscreteDistribution    m_distrb;

    int                     m_num_cells;
    Vectori<ndim, false>    m_cells;
    Array<float, ndim>      m_unit;
};


struct CubeDistribution {
    void set_resolution(const Array<int, 2> &reso);
    void set_mass(const Bitmap3fD &m_radiance);

    FloatC sample_reuse(Vectorf<2, false> &samples) const;
    FloatC pdf(const Vectorf<2, false> &p) const;

    bool                    m_ready = false;
    Array<int, 2>           m_resolution = 0;

    DiscreteDistribution         m_distrb;
    DiscreteDistribution         m_theta_distrb;
    DiscreteDistributionArrayC   m_phi_distrbs;

    int                     m_num_cells;
    Vectori<2, false>    m_cells;
    Array<float, 2>      m_unit;
};


} // namespace psdr

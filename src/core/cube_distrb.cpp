#include <misc/Exception.h>
#include <psdr/core/cube_distrb.h>

namespace psdr
{

template <int ndim>
void HyperCubeDistribution<ndim>::set_resolution(const Array<int, ndim> &reso) {
	if ( m_resolution != reso ) {
		Array<int64_t, ndim> prod_reso;
		prod_reso[ndim - 1] = reso[ndim - 1];
		for ( int i = ndim - 2; i >= 0; --i )
			prod_reso[i] = reso[i]*prod_reso[i + 1];
		PSDR_ASSERT(prod_reso[0] < std::numeric_limits<int>::max());

		m_num_cells = static_cast<int>(prod_reso[0]);
		m_resolution = reso;
		m_unit = rcp(Array<float, ndim>(reso));

		IntC &cur = m_cells[ndim - 1];
		cur = arange<IntC>(m_num_cells);
		for ( int i = 0; i < ndim - 1; ++i ) {
			int denorm = static_cast<int>(prod_reso[i + 1]);
			m_cells[i] = divisor<int>(denorm)(cur);
			cur -= m_cells[i]*denorm;
		}
		m_ready = false;
	}
}


template <int ndim>
void HyperCubeDistribution<ndim>::set_mass(const FloatC &pmf) {
	PSDR_ASSERT(static_cast<int>(slices(pmf)) == m_num_cells);
	m_distrb.init(pmf);
	m_ready = true;
}


template <int ndim>
FloatC HyperCubeDistribution<ndim>::sample_reuse(Vectorf<ndim, false> &samples) const {
	PSDR_ASSERT(m_ready);
	auto [idx, pdf] = m_distrb.sample_reuse<false>(samples[ndim - 1]);
	samples += gather<Vectori<ndim, false>>(m_cells, idx);
	samples *= m_unit;
	return pdf*static_cast<float>(m_num_cells);
}


template <int ndim>
FloatC HyperCubeDistribution<ndim>::pdf(const Vectorf<ndim, false> &p) const {
	PSDR_ASSERT(m_ready);

	auto ip = floor2int<Vectori<ndim, false>, Vectorf<ndim, false>>(p*m_resolution);
	MaskC valid = (ip[0] >= 0 && ip[0] < m_resolution[0]);
	IntC idx = ip[0];
	for ( int i = 1; i < ndim; ++i ) {
		valid &= (ip[i] >= 0 && ip[i] < m_resolution[i]);
		idx = fmadd(idx, m_resolution[i], ip[i]);
	}
	return (gather<FloatC>(m_distrb.pmf(), idx)*static_cast<float>(m_num_cells)) & valid;
}

// Explicit instantiations
template void HyperCubeDistribution<2>::set_resolution(const ScalarVector2i&);
template void HyperCubeDistribution<2>::set_mass(const FloatC&);
template FloatC HyperCubeDistribution<2>::sample_reuse(Vector2fC&) const;
template FloatC HyperCubeDistribution<2>::pdf(const Vector2fC&) const;

template void HyperCubeDistribution<3>::set_resolution(const ScalarVector3i &reso);
template void HyperCubeDistribution<3>::set_mass(const FloatC &pmf);
template FloatC HyperCubeDistribution<3>::sample_reuse(Vector3fC &samples) const;
template FloatC HyperCubeDistribution<3>::pdf(const Vector3fC&) const;

} // namespace psdr

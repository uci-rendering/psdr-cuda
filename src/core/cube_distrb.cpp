#include <misc/Exception.h>
#include <psdr/core/cube_distrb.h>
#include <string> 

// #define USE_2D_DIS

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

void CubeDistribution::set_resolution(const Array<int, 2> &reso) {

	// std::cout << "set reolution for CubeDistribution" << std::endl;
	// std::cout << reso << std::endl;
	if ( m_resolution != reso ) {
		Array<int64_t, 2> prod_reso;
		prod_reso[1] = reso[1];
		for ( int i = 0; i >= 0; --i )
			prod_reso[i] = reso[i]*prod_reso[i + 1];
		PSDR_ASSERT(prod_reso[0] < std::numeric_limits<int>::max());

		m_num_cells = static_cast<int>(prod_reso[0]);
		m_resolution = reso;
		m_unit = rcp(Array<float, 2>(reso));

		IntC &cur = m_cells[1];
		cur = arange<IntC>(m_num_cells);
		for ( int i = 0; i < 1; ++i ) {
			int denorm = static_cast<int>(prod_reso[i + 1]);
			m_cells[i] = divisor<int>(denorm)(cur);
			cur -= m_cells[i]*denorm;
		}
		m_ready = false;
	}
}


void CubeDistribution::set_mass(const Bitmap3fD &m_radiance) {
	int width = m_resolution[0];
	int height = m_resolution[1];

	Vector2fC uv = (m_cells + Vector2fC(.5f, .5f))*m_unit;
	SpectrumC val = m_radiance.eval<false>(uv, false);
	FloatC theta = ((arange<IntC>(width*height) % (height)) + .5f)*(Pi/static_cast<float>(height));
	FloatC radiance_val = rgb2luminance<false>(val)*sin(theta);

	m_distrb.init(radiance_val);


#ifdef USE_2D_DIS
	FloatC theta_pmf = zero<FloatC>(width);
	scatter_add(theta_pmf, radiance_val, m_cells[0]);
	m_theta_distrb.init(theta_pmf);

	std::vector<DDistribution*>      m_DiscreteDistributions;
	for (int i=0; i<width; ++i) {
		std::cout << "\r" << "building envmap: " << i << "/" << width << std::flush;
		FloatC temp_phi_pmf = zero<FloatC>(height);
		IntC phi_id = m_cells[1];
		IntC buf_mask = arange<IntC>(width*height);
		MaskC other_phi = ( i*height <= buf_mask && buf_mask < (i+1)*height);
		scatter_add(temp_phi_pmf, radiance_val, phi_id, other_phi);
		_DiscreteDistribution *distrb = new _DiscreteDistribution();
		distrb->init(temp_phi_pmf);
		distrb->m_id = std::to_string(i);
		m_DiscreteDistributions.push_back(distrb);
	}
	m_phi_distrbs = DiscreteDistributionArrayC::copy(m_DiscreteDistributions.data(), m_DiscreteDistributions.size());
#endif


	m_ready = true;
}


FloatC CubeDistribution::sample_reuse(Vectorf<2, false> &_samples) const {
	PSDR_ASSERT(m_ready);
#ifdef USE_2D_DIS
	Vector2fC samples(_samples);
	auto [idx, pdf_theta] = m_theta_distrb.sample_reuse<false>(samples[1]);
	DiscreteDistributionArrayC arr_phi_distrb = gather<DiscreteDistributionArrayC>(m_phi_distrbs, idx);
	DiscreteRecordC phi_record = arr_phi_distrb->sample_reuse(samples[0]);
	samples[0] = phi_record.rnd;
	_samples[1] = (samples[0] + phi_record.idx)/m_resolution[1];
	_samples[0] = (samples[1] + idx)/m_resolution[0];
	return pdf_theta*phi_record.pdf*static_cast<float>(m_num_cells);
#else
	auto [idx, pdf] = m_distrb.sample_reuse<false>(_samples[1]);
	_samples += gather<Vectori<2, false>>(m_cells, idx);
	_samples *= m_unit;
	return pdf*static_cast<float>(m_num_cells);
#endif
}


FloatC CubeDistribution::pdf(const Vectorf<2, false> &p) const {
	PSDR_ASSERT(m_ready);
	auto ip = floor2int<Vectori<2, false>, Vectorf<2, false>>(p*m_resolution);
	MaskC valid = (ip[0] >= 0 && ip[0] < m_resolution[0]);
	IntC idx = ip[0];
	for ( int i = 1; i < 2; ++i ) {
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

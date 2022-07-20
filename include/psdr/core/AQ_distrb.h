#pragma once

#include <psdr/psdr.h>
#include "pmf.h"

namespace psdr
{

template <int xdim, int ydim, int zdim>
struct Histogram3D
{
    FloatC data;

    Histogram3D() {
        data = zero<FloatC>(xdim*ydim*zdim);
    }

    void update(const Vector3fC &value) {
        IntC idx = IntC(value[0]*FloatC(xdim));
        IntC idy = IntC(value[1]*FloatC(ydim));
        IntC idz = IntC(value[2]*FloatC(zdim));
        IntC idt = idx*xdim*ydim+idy*ydim+idz;
        scatter_add(data, FloatC(1.0f), idt);
    }

    void normalize() {
        data = data * xdim*ydim*zdim / hsum(data) / 1000.f;
    }

    FloatC get_data() {
        return data;
    }

};

template <int ndim>
struct AdaptiveQuadratureDistribution {

    void setup(const Scene &scene, const std::vector<int> &sensor_id, const FloatC &cdfx, const FloatC &cdfy, const FloatC &cdfz, const AQ_Option &option);
	Vector3fC sample(const Vector3fC &rnd, FloatC &pdf);
    FloatC    pdf_mis(const Scene &scene, int sensor_id, const Vector3fC &rnd);

    DiscreteDistribution                       aq_distrb;
    AQLeaf                                     aq_leaf;
    bool                                       aq_edge_direct;
    FloatC                                      aq_sum;

private:
    void __sample_grid(const Scene &scene, const std::vector<int> &sensor_id, int npass, float RMSE_wt);
    FloatC __eval(const AQLeaf &sample_leaf, const Vector3fC &sample);

    FloatC __FZ(const AQLeaf &sample_leaf, const FloatC &rndz);
    FloatC __FY(const AQLeaf &sample_leaf, const FloatC fix_z, const FloatC &rndy);
    FloatC __FX(const AQLeaf &sample_leaf, const FloatC fix_y, const FloatC fix_z, const FloatC &rndx);

};

} // namespace psdr

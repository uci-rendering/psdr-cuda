#include <misc/Exception.h>
#include <psdr/core/AQ_distrb.h>
#include <psdr/core/sampler.h>
#include <psdr/core/ray.h>
#include <psdr/core/intersection.h>
#include <psdr/core/transform.h>
#include <psdr/bsdf/bsdf.h>
#include <psdr/emitter/emitter.h>
#include <psdr/shape/mesh.h>
#include <psdr/scene/scene.h>
#include <psdr/sensor/perspective.h>

#include <cuda/host/NEE.cuh>

namespace psdr
{

inline static FloatC Boundary_NEE(const Scene &scene, int active_sensor, const Vector3fC &sample3, bool m_edge_direct) {
    BoundarySegSampleDirect bss;
    if (m_edge_direct) {
        bss = scene.sample_emitter_ray(sample3);
    } else {
        bss = scene.sample_edge_ray(sample3);
    }

    MaskC valid = bss.is_valid;
    const Vector3fC &_p0    = detach(bss.p0);
    Vector3fC       _p2, _dir;
    _dir = bss.p2;
    IntersectionC _its2 = scene.ray_intersect<false>(RayC(_p0, _dir), valid); 
    IntersectionC _its1 = scene.ray_intersect<false>(RayC(_p0, -_dir), valid);
    _p2 = _its2.p;

    valid &= _its2.is_emitter(valid) && _its2.is_valid() && norm(_its2.p - _p2) < ShadowEpsilon && _its1.is_valid();
    Vector3fC &_p1 = _its1.p;
    SensorDirectSampleC sds = scene.m_sensors[active_sensor]->sample_direct(_p1);
    valid &= sds.is_valid;
    RayC camera_ray = scene.m_sensors[active_sensor]->sample_primary_ray(sds.q);
    IntersectionC its1 = scene.ray_intersect<false>(camera_ray, valid);
    valid &= its1.is_valid() && norm(its1.p - _p1) < ShadowEpsilon;
    FloatC      dist    = norm(_p2 - _p1),
                cos2    = abs(dot(_its2.n, -_dir));
    Vector3fC   e       = cross(bss.edge, _dir);
    FloatC      sinphi  = norm(e);
    Vector3fC   proj    = normalize(cross(e, _its2.n));
    FloatC      sinphi2 = norm(cross(_dir, proj));
    FloatC      base_v  = (_its1.t/dist)*(sinphi/sinphi2)*cos2;
    valid &= (sinphi > Epsilon) && (sinphi2 > Epsilon);
    Vector3fC d0 = -camera_ray.d;
    Vector3fC d0_local = _its1.sh_frame.to_local(d0);
    SpectrumC bsdf_val = _its1.shape->bsdf(valid)->eval(_its1, d0_local, valid);
    FloatC correction = abs((_its1.wi.z()*dot(d0, _its1.n))/(d0_local.z()*dot(_dir, _its1.n)));
    masked(bsdf_val, valid) *= correction;
    SpectrumC value0 = (bsdf_val*_its2.Le(valid)*(base_v*sds.sensor_val/bss.pdf)) & valid;
    masked(value0, ~enoki::isfinite<SpectrumC>(value0)) = 0.f;
    value0 = abs(value0);

    float eps;
    if (m_edge_direct) {
        eps = 0.0f;
    } else {
        eps = 0.0f;
    }
    FloatC buf = max(hmax(value0), eps);
    return buf;
}

static inline FloatC linear_eval(const Vector2fC &data, const FloatC &value) {
    return value * (data[1] - data[0]) + data[0];
}

static inline FloatC linear_int(const Vector2fC &poly, const FloatC &value) {
    return (poly[0] + (value * (poly[1] - poly[0]) + poly[0])) * value / FloatC(2.0f);
}

static inline FloatC LagrangeND(const Vectorf<8, false> &data, const Vector3fC &value) {
    FloatC result = 0.0f;
    Vector2fC data3D;
    data3D[0] = linear_eval(Vector2fC(linear_eval(Vector2fC(data[0], data[1]), value[2]), linear_eval(Vector2fC(data[2], data[3]), value[2])), value[1]);
    data3D[1] = linear_eval(Vector2fC(linear_eval(Vector2fC(data[4], data[5]), value[2]), linear_eval(Vector2fC(data[6], data[7]), value[2])), value[1]);
    return linear_eval(data3D, value[0]);
}

static inline FloatC finverse_norm(const FloatC &aa, const FloatC &bb, const FloatC &yi) {
    Vector3fC abc(aa / (aa+bb), bb / (aa+bb), -yi);
    FloatC delta = abc[1]*abc[1]-4.0f*abc[0]*abc[2];
    return select(abs(abc[0]) > 0.0001, (-abc[1] + sqrt(delta)) / (2.0f*abc[0]), -abc[2]/abc[1]);
}

template <int ndim>
FloatC AdaptiveQuadratureDistribution<ndim>::__FZ(const AQLeaf &sample_leaf, const FloatC &rndz) {
    FloatC val1 = linear_int(Vector2fC(linear_int(Vector2fC(sample_leaf.poly[0], sample_leaf.poly[4]), 1.0f), linear_int(Vector2fC(sample_leaf.poly[2], sample_leaf.poly[6]), 1.0f)), 1.0f);
    FloatC val2 = linear_int(Vector2fC(linear_int(Vector2fC(sample_leaf.poly[1], sample_leaf.poly[5]), 1.0f), linear_int(Vector2fC(sample_leaf.poly[3], sample_leaf.poly[7]), 1.0f)), 1.0f);
    FloatC aa = (val2 - val1) / FloatC(2.0f);
    FloatC bb = val1;
    return finverse_norm(aa, bb, rndz);
}

template <int ndim>
FloatC AdaptiveQuadratureDistribution<ndim>::__FY(const AQLeaf &sample_leaf, const FloatC fix_z, const FloatC &rndy) {
    FloatC val11 = linear_eval(Vector2fC(sample_leaf.poly[0], sample_leaf.poly[1]), fix_z);
    FloatC val12 = linear_eval(Vector2fC(sample_leaf.poly[2], sample_leaf.poly[3]), fix_z);
    FloatC val21 = linear_eval(Vector2fC(sample_leaf.poly[4], sample_leaf.poly[5]), fix_z);
    FloatC val22 = linear_eval(Vector2fC(sample_leaf.poly[6], sample_leaf.poly[7]), fix_z);
    FloatC aa = ((val12 - val11 + val22 - val21) / FloatC(4.0f));
    FloatC bb = (val11 + val21) / FloatC(2.0f);
    return finverse_norm(aa, bb, rndy);
}

template <int ndim>
FloatC AdaptiveQuadratureDistribution<ndim>::__FX(const AQLeaf &sample_leaf, const FloatC fix_y, const FloatC fix_z, const FloatC &rndx) {
    Vector2fC para1(linear_eval(Vector2fC(sample_leaf.poly[0], sample_leaf.poly[1]), fix_z), linear_eval(Vector2fC(sample_leaf.poly[2], sample_leaf.poly[3]), fix_z));
    FloatC val1 = linear_eval(para1, fix_y);
    Vector2fC para2(linear_eval(Vector2fC(sample_leaf.poly[4], sample_leaf.poly[5]), fix_z), linear_eval(Vector2fC(sample_leaf.poly[6], sample_leaf.poly[7]), fix_z));
    FloatC val2 = linear_eval(para2, fix_y);
    FloatC aa = (val2 - val1) / FloatC(2.0f);
    FloatC bb = val1;
    return finverse_norm(aa, bb, rndx);
}

template <int ndim>
void AdaptiveQuadratureDistribution<ndim>::__sample_grid(const Scene &scene, const std::vector<int> &sensor_id, int npass, float RMSE_wt) {
    size_t num_samples = slices(aq_leaf);
    FloatC ymean = zero<FloatC>(num_samples);

#if 1
    size_t sample_size = num_samples*npass;
    IntC idx = arange<IntC>(sample_size) / npass;

    Sampler sampler;
    sampler.seed(arange<UInt64C>(sample_size));

    Vector3fC rnd_base = sampler.next_nd<3, false>();
    Vector3fC mid_base = gather<Vector3fC>(aq_leaf.p0, idx) + gather<Vector3fC>(aq_leaf.p1 - aq_leaf.p0, idx) * rnd_base;
    FloatC fx = Boundary_NEE(scene, sensor_id[0], mid_base, aq_edge_direct);

    scatter_add(ymean, fx, idx);
    ymean  /= static_cast<float>(npass);
    // aq_leaf.poly = Vectorf<8, false>(ymean, ymean, ymean, ymean, ymean, ymean, ymean, ymean);

    FloatC x0mean = zero<FloatC>(num_samples);
    FloatC x1mean = zero<FloatC>(num_samples);
    FloatC x2mean = zero<FloatC>(num_samples);
    FloatC s0 = zero<FloatC>(num_samples);
    FloatC s1 = zero<FloatC>(num_samples);
    FloatC s2 = zero<FloatC>(num_samples);
    FloatC s00 = zero<FloatC>(num_samples);
    FloatC s11 = zero<FloatC>(num_samples);
    FloatC s22 = zero<FloatC>(num_samples);
    FloatC s01 = zero<FloatC>(num_samples);
    FloatC s02 = zero<FloatC>(num_samples);
    FloatC s12 = zero<FloatC>(num_samples);
    FloatC s0y = zero<FloatC>(num_samples);
    FloatC s1y = zero<FloatC>(num_samples);
    FloatC s2y = zero<FloatC>(num_samples);
    FloatC sy = zero<FloatC>(num_samples);

    scatter_add(x0mean, rnd_base[0], idx);
    x0mean /= static_cast<float>(npass);
    scatter_add(x1mean, rnd_base[1], idx);
    x1mean /= static_cast<float>(npass);
    scatter_add(x2mean, rnd_base[2], idx);
    x2mean /= static_cast<float>(npass);

    scatter_add(s0, rnd_base[0], idx);
    scatter_add(s1, rnd_base[1], idx);
    scatter_add(s2, rnd_base[2], idx);

    scatter_add(s00, rnd_base[0]*rnd_base[0], idx);
    scatter_add(s11, rnd_base[1]*rnd_base[1], idx);
    scatter_add(s22, rnd_base[2]*rnd_base[2], idx);
    scatter_add(s01, rnd_base[0]*rnd_base[1], idx);
    scatter_add(s02, rnd_base[0]*rnd_base[2], idx);
    scatter_add(s12, rnd_base[1]*rnd_base[2], idx);
    scatter_add(s0y, rnd_base[0]*fx, idx);
    scatter_add(s1y, rnd_base[1]*fx, idx);
    scatter_add(s2y, rnd_base[2]*fx, idx);
    scatter_add(sy, fx, idx);


    FloatC Sx0x0 = s00 - s0*s0 / static_cast<float>(npass);
    FloatC Sx1x1 = s11 - s1*s1 / static_cast<float>(npass);
    FloatC Sx2x2 = s22 - s2*s2 / static_cast<float>(npass);
    FloatC Sx0x1 = s01 - s0*s1 / static_cast<float>(npass);
    FloatC Sx0x2 = s02 - s0*s2 / static_cast<float>(npass);
    FloatC Sx1x2 = s12 - s1*s2 / static_cast<float>(npass);

    Vector3fC MY(s0y - s0*sy / static_cast<float>(npass), 
                 s1y - s1*sy / static_cast<float>(npass), 
                 s2y - s2*sy / static_cast<float>(npass));

    Matrix3x3fC MX( Vector3fC(Sx0x0, Sx0x1, Sx0x2),
                    Vector3fC(Sx0x1, Sx1x1, Sx1x2),
                    Vector3fC(Sx0x2, Sx1x2, Sx2x2));

    Matrix3fC bb = inverse(MX)*MY;

    FloatC b1 = bb[0][0];
    FloatC b2 = bb[1][0];
    FloatC b3 = bb[2][0];
    FloatC b0 = ymean - b1*x0mean - b2*x1mean - b3*x2mean;

    // RMSE balance
    if (RMSE_wt > 0.00001f) {
        FloatC interp = gather<FloatC>(b1, idx)*mid_base.x() + gather<FloatC>(b2, idx)*mid_base.y() + gather<FloatC>(b3, idx)*mid_base.z() + gather<FloatC>(b0, idx);
        FloatC RMSE = zero<FloatC>(num_samples);
        scatter_add(RMSE, abs(interp - fx) / static_cast<float>(npass), idx);
        b0 += RMSE * RMSE_wt;
    }

    FloatC poly1 = b0;
    FloatC poly2 = b0 + b3;
    FloatC poly3 = b0 + b2;
    FloatC poly4 = b0 + b2 + b3;
    FloatC poly5 = b0 + b1;
    FloatC poly6 = b0 + b1 + b3;
    FloatC poly7 = b0 + b1 + b2;
    FloatC poly8 = b0 + b1 + b2 + b3;

    // y = b1*x + b2*y + b3*z + b0

    MaskC p1_mask = poly1 < 0.0f;
    MaskC p2_mask = poly2 < 0.0f;
    MaskC p3_mask = poly3 < 0.0f;
    MaskC p4_mask = poly4 < 0.0f;

    MaskC p5_mask = poly5 < 0.0f;
    MaskC p6_mask = poly6 < 0.0f;
    MaskC p7_mask = poly7 < 0.0f;
    MaskC p8_mask = poly8 < 0.0f;

    FloatC neg_poly1 = select(p1_mask, poly1, 0.0f);
    FloatC neg_poly2 = select(p2_mask, poly2, 0.0f);
    FloatC neg_poly3 = select(p3_mask, poly3, 0.0f);
    FloatC neg_poly4 = select(p4_mask, poly4, 0.0f);

    FloatC neg_poly5 = select(p5_mask, poly5, 0.0f);
    FloatC neg_poly6 = select(p6_mask, poly6, 0.0f);
    FloatC neg_poly7 = select(p7_mask, poly7, 0.0f);
    FloatC neg_poly8 = select(p8_mask, poly8, 0.0f);

    poly1 = poly1 + neg_poly8;
    poly2 = poly2 + neg_poly7;
    poly3 = poly3 + neg_poly6;
    poly4 = poly4 + neg_poly5;

    poly5 = poly5 + neg_poly4;
    poly6 = poly6 + neg_poly3;
    poly7 = poly7 + neg_poly2;
    poly8 = poly8 + neg_poly1;

    poly1 = max(0.0f, poly1);
    poly2 = max(0.0f, poly2);
    poly3 = max(0.0f, poly3);
    poly4 = max(0.0f, poly4);
    poly5 = max(0.0f, poly5);
    poly6 = max(0.0f, poly6);
    poly7 = max(0.0f, poly7);
    poly8 = max(0.0f, poly8);

    aq_leaf.poly = Vectorf<8, false>(poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8);

#else
    FloatC x0mean = zero<FloatC>(num_samples);
    FloatC x1mean = zero<FloatC>(num_samples);
    FloatC x2mean = zero<FloatC>(num_samples);
    FloatC s0 = zero<FloatC>(num_samples);
    FloatC s1 = zero<FloatC>(num_samples);
    FloatC s2 = zero<FloatC>(num_samples);
    FloatC s00 = zero<FloatC>(num_samples);
    FloatC s11 = zero<FloatC>(num_samples);
    FloatC s22 = zero<FloatC>(num_samples);
    FloatC s01 = zero<FloatC>(num_samples);
    FloatC s02 = zero<FloatC>(num_samples);
    FloatC s12 = zero<FloatC>(num_samples);
    FloatC s0y = zero<FloatC>(num_samples);
    FloatC s1y = zero<FloatC>(num_samples);
    FloatC s2y = zero<FloatC>(num_samples);
    FloatC sy = zero<FloatC>(num_samples);

    Sampler sampler;
    sampler.seed(arange<UInt64C>(num_samples));
    std::vector<Vector3fC> training_data_x;
    std::vector<FloatC>    training_data_y;

    for (int i=0; i<npass; ++i) {
        Vector3fC rnd = sampler.next_nd<3, false>();
        Vector3fC mid = aq_leaf.p0 + (aq_leaf.p1 - aq_leaf.p0) * rnd;
        FloatC fx = Boundary_NEE(scene, sensor_id, mid);

        if (RMSE_wt > 0.00001f) {
            training_data_x.push_back(rnd);
            training_data_y.push_back(fx);
        }

        ymean += fx / static_cast<float>(npass);
        x0mean += rnd[0] / static_cast<float>(npass);
        x1mean += rnd[1] / static_cast<float>(npass);
        x2mean += rnd[2] / static_cast<float>(npass);
        s0  += rnd[0];
        s1  += rnd[1];
        s2  += rnd[2];
        s00 += rnd[0]*rnd[0];
        s11 += rnd[1]*rnd[1];
        s22 += rnd[2]*rnd[2];
        s01 += rnd[0]*rnd[1];
        s02 += rnd[0]*rnd[2];
        s12 += rnd[1]*rnd[2];
        s0y += rnd[0]*fx;
        s1y += rnd[1]*fx;
        s2y += rnd[2]*fx;
        sy  += fx;
    }

    FloatC Sx0x0 = s00 - s0*s0 / static_cast<float>(npass);
    FloatC Sx1x1 = s11 - s1*s1 / static_cast<float>(npass);
    FloatC Sx2x2 = s22 - s2*s2 / static_cast<float>(npass);
    FloatC Sx0x1 = s01 - s0*s1 / static_cast<float>(npass);
    FloatC Sx0x2 = s02 - s0*s2 / static_cast<float>(npass);
    FloatC Sx1x2 = s12 - s1*s2 / static_cast<float>(npass);

    Vector3fC MY(s0y - s0*sy / static_cast<float>(npass), 
                 s1y - s1*sy / static_cast<float>(npass), 
                 s2y - s2*sy / static_cast<float>(npass));

    Matrix3x3fC MX( Vector3fC(Sx0x0, Sx0x1, Sx0x2),
                    Vector3fC(Sx0x1, Sx1x1, Sx1x2),
                    Vector3fC(Sx0x2, Sx1x2, Sx2x2));

    Matrix3fC bb = inverse(MX)*MY;

    FloatC b1 = bb[0][0];
    FloatC b2 = bb[1][0];
    FloatC b3 = bb[2][0];
    FloatC b0 = ymean - b1*x0mean - b2*x1mean - b3*x2mean;

    // RMSE balance
    if (RMSE_wt > 0.00001f) {
        FloatC RMSE = zero<FloatC>(num_samples);    
        for (int i=0; i<npass; ++i) {
            Vector3fC rnd = training_data_x[i];
            Vector3fC mid = aq_leaf.p0 + (aq_leaf.p1 - aq_leaf.p0) * rnd;
            FloatC interp = b1*mid.x() + b2*mid.y() + b3*mid.z() + b0;
            FloatC fx = training_data_y[i];
            RMSE += abs(interp - fx) / static_cast<float>(npass);
        }
        b0 += RMSE * RMSE_wt;
    }

    FloatC poly1 = b0;
    FloatC poly2 = b0 + b3;
    FloatC poly3 = b0 + b2;
    FloatC poly4 = b0 + b2 + b3;
    FloatC poly5 = b0 + b1;
    FloatC poly6 = b0 + b1 + b3;
    FloatC poly7 = b0 + b1 + b2;
    FloatC poly8 = b0 + b1 + b2 + b3;

    MaskC p1_mask = poly1 < 0.0f;
    MaskC p2_mask = poly2 < 0.0f;
    MaskC p3_mask = poly3 < 0.0f;
    MaskC p4_mask = poly4 < 0.0f;

    MaskC p5_mask = poly5 < 0.0f;
    MaskC p6_mask = poly6 < 0.0f;
    MaskC p7_mask = poly7 < 0.0f;
    MaskC p8_mask = poly8 < 0.0f;

    FloatC neg_poly1 = select(p1_mask, poly1, 0.0f);
    FloatC neg_poly2 = select(p2_mask, poly2, 0.0f);
    FloatC neg_poly3 = select(p3_mask, poly3, 0.0f);
    FloatC neg_poly4 = select(p4_mask, poly4, 0.0f);

    FloatC neg_poly5 = select(p5_mask, poly5, 0.0f);
    FloatC neg_poly6 = select(p6_mask, poly6, 0.0f);
    FloatC neg_poly7 = select(p7_mask, poly7, 0.0f);
    FloatC neg_poly8 = select(p8_mask, poly8, 0.0f);

    poly1 = poly1 + neg_poly8;
    poly2 = poly2 + neg_poly7;
    poly3 = poly3 + neg_poly6;
    poly4 = poly4 + neg_poly5;

    poly5 = poly5 + neg_poly4;
    poly6 = poly6 + neg_poly3;
    poly7 = poly7 + neg_poly2;
    poly8 = poly8 + neg_poly1;

    poly1 = max(0.0f, poly1);
    poly2 = max(0.0f, poly2);
    poly3 = max(0.0f, poly3);
    poly4 = max(0.0f, poly4);
    poly5 = max(0.0f, poly5);
    poly6 = max(0.0f, poly6);
    poly7 = max(0.0f, poly7);
    poly8 = max(0.0f, poly8);

    aq_leaf.poly = Vectorf<8, false>(poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8);
    // aq_leaf.poly = Vectorf<8, false>(ymean, ymean, ymean, ymean, ymean, ymean, ymean, ymean);
#endif
}

template <int ndim>
void AdaptiveQuadratureDistribution<ndim>::setup(const Scene &scene, const std::vector<int> &sensor_id, const FloatC &cdfx, const FloatC &cdfy, const FloatC &cdfz, const AQ_Option &option) {
    if (aq_edge_direct) {
        std::cout << "Begin AQ emitter guiding with init grid size =" << slices(cdfx) << " " << slices(cdfy) << " " << slices(cdfz) << std::endl;
    } else {
        std::cout << "Begin AQ direction guiding with init grid size =" << slices(cdfx) << " " << slices(cdfy) << " " << slices(cdfz) << std::endl;
    }

    size_t dimx = slices(cdfx);
    size_t dimy = slices(cdfy);
    size_t dimz = slices(cdfz);
    size_t init_size = dimx*dimy*dimz;
    size_t max_tree_size = option.max_memory;

    if (init_size > max_tree_size){
        PSDR_ASSERT(0);
    }

    tree3D tree_leaf = zero<tree3D>(max_tree_size);
    psdr_cuda::init_tree(cdfx.data(), 
                         cdfy.data(),
                         cdfz.data(),
                         dimx,
                         dimy,
                         dimz,
                         tree_leaf.p0.x().data(),
                         tree_leaf.p0.y().data(),
                         tree_leaf.p0.z().data(),
                         tree_leaf.p1.x().data(),
                         tree_leaf.p1.y().data(),
                         tree_leaf.p1.z().data());

    int fix_size = init_size;
    int cut_dim = 1;

    for (int i=0; i<option.max_depth; ++i) {
        // if (aq_edge_direct) {
        //     std::cout << "we dont cut direct for now" << std::endl;
        //     break;
        // }

        Vector3fC eval_rnd_buf = zero<Vector3fC>(fix_size * 36);
        psdr_cuda::generate_eval_point(fix_size,
                                        tree_leaf.p0.x().data(),
                                        tree_leaf.p0.y().data(),
                                        tree_leaf.p0.z().data(),
                                        tree_leaf.p1.x().data(),
                                        tree_leaf.p1.y().data(),
                                        tree_leaf.p1.z().data(),
                                        eval_rnd_buf.x().data(),
                                        eval_rnd_buf.y().data(),
                                        eval_rnd_buf.z().data(), i);
        FloatC error_value = Boundary_NEE(scene, sensor_id[0], eval_rnd_buf, aq_edge_direct);
        cuda_eval(); cuda_sync();
        int app_size = psdr_cuda::cut_grid( error_value.data(),
                             tree_leaf.p0.x().data(),
                             tree_leaf.p0.y().data(),
                             tree_leaf.p0.z().data(),
                             tree_leaf.p1.x().data(),
                             tree_leaf.p1.y().data(),
                             tree_leaf.p1.z().data(),
                             eval_rnd_buf.x().data(),
                             eval_rnd_buf.y().data(),
                             eval_rnd_buf.z().data(),
                             fix_size, option.thold, option.wt1);
        fix_size += app_size;
        if (fix_size > max_tree_size){
            PSDR_ASSERT(0);
        }
        std::cout << "\r"  << "Depth " << i << " with grid: " << fix_size;
    }
    std::cout << std::endl;
    aq_leaf = zero<AQLeaf>(slices(tree_leaf));
    aq_leaf.p0 = tree_leaf.p0;
    aq_leaf.p1 = tree_leaf.p1;

    FloatC aq_leaf_area = (aq_leaf.p1-aq_leaf.p0).x()*
                          (aq_leaf.p1-aq_leaf.p0).y()*
                          (aq_leaf.p1-aq_leaf.p0).z();
    __sample_grid(scene, sensor_id, option.final_spp, option.RMSE_wt);

    aq_sum = hsum(hsum(aq_leaf.poly) / FloatC(8.0) * aq_leaf_area);
    aq_leaf.poly /= aq_sum;
    aq_leaf.poly = option.eps + (1.f - option.eps)*aq_leaf.poly;
    FloatC pmf = hsum(aq_leaf.poly) / FloatC(8.0) * aq_leaf_area;

    aq_distrb.init(pmf);
    cuda_eval(); cuda_sync();
}

template <int ndim>
Vector3fC AdaptiveQuadratureDistribution<ndim>::sample(const Vector3fC &_rnd, FloatC &pdf) {

    Vector3fC rnd(_rnd);
    auto [idx, _] = aq_distrb.sample_reuse<false>(rnd.x());
    AQLeaf saq = gather<AQLeaf>(aq_leaf, idx);

    FloatC rndz = __FZ(saq, rnd[2]);
    FloatC rndy = __FY(saq, rndz, rnd[1]);
    FloatC rndx = __FX(saq, rndy, rndz, rnd[0]);
    // FloatC rndz = rnd[2];
    // FloatC rndy = rnd[1];
    // FloatC rndx = rnd[0];
    Vector3fC result(saq.p0.x() + rndx*(saq.p1.x()-saq.p0.x()), 
                     saq.p0.y() + rndy*(saq.p1.y()-saq.p0.y()),
                     saq.p0.z() + rndz*(saq.p1.z()-saq.p0.z()));

    pdf = LagrangeND(saq.poly, Vector3fC(rndx, rndy, rndz)) / aq_distrb.m_sum;
    return result;

}

template <int ndim>
FloatC AdaptiveQuadratureDistribution<ndim>::pdf_mis(const Scene &scene, int sensor_id, const Vector3fC &_rnd) {
    FloatC estimate = Boundary_NEE(scene, sensor_id, _rnd, aq_edge_direct);
    std::cout << estimate << std::endl;
    return estimate;
}

template FloatC AdaptiveQuadratureDistribution<3>::pdf_mis(const Scene &scene, int sensor_id, const Vector3fC &rnd);

template Vector3fC AdaptiveQuadratureDistribution<3>::sample(const Vector3fC &rnd, FloatC &pdf);
template void AdaptiveQuadratureDistribution<3>::setup(const Scene &scene, const std::vector<int> &sensor_id, const FloatC &cdfx, const FloatC &cdfy, const FloatC &cdfz, const AQ_Option &option);

} // namespace psdr

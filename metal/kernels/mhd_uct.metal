#include "../metal_config.h"

#if NSUBGRID != 1 && NSUBGRID != 2
#error NSUBGRID must be 1 or 2
#endif

#define MR_NSUBGRIDP2 (NSUBGRID + 2)
#define MR_SUBGRIDSIZE (MR_NSUBGRIDP2 * MR_NSUBGRIDP2 * MR_NSUBGRIDP2)
#define MR_M (2 * NSUBGRID)
#define MR_STENCIL (MR_M + 4)
#define MR_TRACE (MR_M + 2)
#define MR_LOCAL_CELLS (MR_STENCIL * MR_STENCIL * MR_STENCIL)
#define MR_REFINED_CELLS (MR_TRACE * MR_TRACE * MR_TRACE)
#define MR_BF_COMPONENT ((MR_STENCIL - 1) * MR_STENCIL * MR_STENCIL)
#define MR_FACE_CELLS ((MR_M + 1) * MR_M * MR_M)
#define MR_ALL_FACES (3 * MR_FACE_CELLS)
#define MR_SHELL_CELLS (4 * MR_M * (MR_M + 1))
#define MR_ALL_SHELLS (3 * MR_SHELL_CELLS)
#define MR_EDGE_CELLS (MR_M * (MR_M + 1) * (MR_M + 1))
#define MR_ALL_EDGES (3 * MR_EDGE_CELLS)
#define MR_INTERIOR_CELLS (MR_M * MR_M * MR_M)
#define MR_NSUBGRID_CELLS (NSUBGRID * NSUBGRID * NSUBGRID)
#define MR_LOCAL_BASE 0
#define MR_BF_BASE (5 * MR_LOCAL_CELLS)
#define MR_STATE_FLOATS (MR_BF_BASE + 3 * MR_BF_COMPONENT)
#define MR_FACE_BASE (MR_BF_BASE + 3 * MR_BF_COMPONENT)
#define MR_FACE_BUFFER (8 * MR_FACE_CELLS)
#define MR_FLUX_FLOATS (5 * MR_ALL_FACES)
#define MR_RECORD_FLOATS (6 * MR_ALL_FACES)
#define MR_SHELL_FREE (6 * MR_FACE_BUFFER - MR_FLUX_FLOATS - MR_RECORD_FLOATS)
#define MR_SHELL_RECORD_FLOATS (6 * MR_ALL_SHELLS)
#define MR_SHELL_SPILL (MR_SHELL_RECORD_FLOATS - MR_SHELL_FREE)
#define MR_EMF_BASE MR_SHELL_SPILL
#define MR_SMEM_FLOATS (MR_FACE_BASE + 6 * MR_FACE_BUFFER)
#define MR_HLLD_CORNER_FIELDS 7
#define MR_HLLD_CORNER_FLOATS (4 * MR_HLLD_CORNER_FIELDS * MR_EDGE_CELLS)
#define MR_HLLD_EMF_BASE (MR_FACE_BASE + MR_HLLD_CORNER_FLOATS)
struct mr_uct_record_t {
    float aL;
    float dL;
    float dR;
    float vt1;
    float vt2;
    float Bn;
};

struct mr_trace_t {
    primitive_t cell;
    primitive_t sx;
    primitive_t sy;
    primitive_t sz;
    float AL;
    float AR;
    float BL;
    float BR;
    float CL;
    float CR;
};

struct mr_edge_pair_t {
    float left;
    float right;
};

bool mr_switch_to_llf(float rmin, float pmin, float switch_dmin, float switch_pmin) {
    // The pressure threshold takes precedence when both switches are enabled.
    bool use_llf = false;
    if (switch_dmin > 0.0f) use_llf = rmin < switch_dmin;
    if (switch_pmin > 0.0f) use_llf = pmin < switch_pmin;
    return use_llf;
}

int mr_nbor_get(device const int *nbor, int subgrid_idx, int ind_nbor) {
    return nbor[(subgrid_idx - 1) * MR_SUBGRIDSIZE + ind_nbor - 1];
}

int mr_local_index(int i, int j, int k) {
    return i + MR_STENCIL * (j + MR_STENCIL * k);
}

float mr_local_get(threadgroup const float *s, int field, int i, int j, int k) {
    return s[MR_LOCAL_BASE + field * MR_LOCAL_CELLS + mr_local_index(i, j, k)];
}

void mr_local_set(threadgroup float *s, int field, int i, int j, int k, float value) {
    s[MR_LOCAL_BASE + field * MR_LOCAL_CELLS + mr_local_index(i, j, k)] = value;
}

int mr_refined_index(int i, int j, int k) {
    return (i - 1) + MR_TRACE * ((j - 1) + MR_TRACE * (k - 1));
}

bool mr_refined_get(threadgroup const bool *refined, int i, int j, int k) {
    return refined[mr_refined_index(i, j, k)];
}

void mr_refined_set(threadgroup bool *refined, int i, int j, int k, bool value) {
    refined[mr_refined_index(i, j, k)] = value;
}

int mr_bx_index(int i, int j, int k) {
    return i - 1 + (MR_STENCIL - 1) * (j + MR_STENCIL * k);
}

int mr_by_index(int i, int j, int k) {
    return i + MR_STENCIL * (j - 1 + (MR_STENCIL - 1) * k);
}

int mr_bz_index(int i, int j, int k) {
    return i + MR_STENCIL * (j + MR_STENCIL * (k - 1));
}

float mr_bf_get(threadgroup const float *s, int component, int i, int j, int k) {
    int index = component == 0 ? mr_bx_index(i, j, k) : component == 1 ? mr_by_index(i, j, k) : mr_bz_index(i, j, k);
    return s[MR_BF_BASE + component * MR_BF_COMPONENT + index];
}

void mr_bf_set(threadgroup float *s, int component, int i, int j, int k, float value) {
    int index = component == 0 ? mr_bx_index(i, j, k) : component == 1 ? mr_by_index(i, j, k) : mr_bz_index(i, j, k);
    s[MR_BF_BASE + component * MR_BF_COMPONENT + index] = value;
}

int mr_face_index(int orientation, int i, int j, int k) {
    if (orientation == 0) return i + (MR_M + 1) * (j + MR_M * k);
    if (orientation == 1) return i + MR_M * (j + (MR_M + 1) * k);
    return i + MR_M * (j + MR_M * k);
}

float mr_face_get(threadgroup const float *s, int buffer, int field, int index) {
    return s[MR_FACE_BASE + buffer * MR_FACE_BUFFER + field * MR_FACE_CELLS + index];
}

void mr_face_set(threadgroup float *s, int buffer, int field, int index, float value) {
    s[MR_FACE_BASE + buffer * MR_FACE_BUFFER + field * MR_FACE_CELLS + index] = value;
}

primitive_t mr_face_load(threadgroup const float *s, int buffer, int index) {
    primitive_t q;
    q.density = mr_face_get(s, buffer, 0, index);
    q.velocity_x = mr_face_get(s, buffer, 1, index);
    q.velocity_y = mr_face_get(s, buffer, 2, index);
    q.velocity_z = mr_face_get(s, buffer, 3, index);
    q.pressure = mr_face_get(s, buffer, 4, index);
    q.Bx = mr_face_get(s, buffer, 5, index);
    q.By = mr_face_get(s, buffer, 6, index);
    q.Bz = mr_face_get(s, buffer, 7, index);
    return q;
}

void mr_face_store(threadgroup float *s, int buffer, int index, primitive_t q) {
    mr_face_set(s, buffer, 0, index, q.density);
    mr_face_set(s, buffer, 1, index, q.velocity_x);
    mr_face_set(s, buffer, 2, index, q.velocity_y);
    mr_face_set(s, buffer, 3, index, q.velocity_z);
    mr_face_set(s, buffer, 4, index, q.pressure);
    mr_face_set(s, buffer, 5, index, q.Bx);
    mr_face_set(s, buffer, 6, index, q.By);
    mr_face_set(s, buffer, 7, index, q.Bz);
}

float mr_flux_get(threadgroup const float *s, int field, int face) {
    return s[MR_FACE_BASE + field * MR_ALL_FACES + face];
}

void mr_flux_set(threadgroup float *s, int field, int face, float value) {
    s[MR_FACE_BASE + field * MR_ALL_FACES + face] = value;
}

float mr_record_value(mr_uct_record_t r, int field) {
    if (field == 0) return r.aL;
    if (field == 1) return r.dL;
    if (field == 2) return r.dR;
    if (field == 3) return r.vt1;
    if (field == 4) return r.vt2;
    return r.Bn;
}

void mr_record_set_value(thread mr_uct_record_t &r, int field, float value) {
    if (field == 0) r.aL = value;
    else if (field == 1) r.dL = value;
    else if (field == 2) r.dR = value;
    else if (field == 3) r.vt1 = value;
    else if (field == 4) r.vt2 = value;
    else r.Bn = value;
}

int mr_product_face_index(int orientation, int i, int j, int k) {
    return orientation * MR_INTERIOR_CELLS + i + MR_M * (j + MR_M * k);
}

ulong mr_product_index(int subgrid_idx, int head_idx, int field, int face) {
    return (((ulong)(subgrid_idx - head_idx) * MR_UCT_PRODUCT_FIELDS + field) * MR_UCT_PRODUCT_FACES + face);
}

void mr_product_set(device float *product, int subgrid_idx, int head_idx, int field, int face, float value) {
    product[mr_product_index(subgrid_idx, head_idx, field, face)] = value;
}

float mr_product_get(device const float *product, int subgrid_idx, int head_idx, int field, int face) {
    return product[mr_product_index(subgrid_idx, head_idx, field, face)];
}

int mr_positive_mod(int value, int modulus) {
    int result = value % modulus;
    return result < 0 ? result + modulus : result;
}

int mr_product_owner_face(device const int *nbor, int subgrid_idx, int orientation, int normal, int t1, int t2, thread int &owner_subgrid) {
    int i = orientation == 0 ? normal : t1;
    int j = orientation == 1 ? normal : orientation == 0 ? t1 : t2;
    int k = orientation == 2 ? normal : t2;
    int i_sg = (i + 2) / 2;
    int j_sg = (j + 2) / 2;
    int k_sg = (k + 2) / 2;
    int source_idx = mr_nbor_get(nbor, subgrid_idx, 1 + i_sg + MR_NSUBGRIDP2 * j_sg + MR_NSUBGRIDP2 * MR_NSUBGRIDP2 * k_sg);
    owner_subgrid = (source_idx - 1) / MR_NSUBGRID_CELLS + 1;
    return mr_product_face_index(orientation, mr_positive_mod(i, MR_M), mr_positive_mod(j, MR_M), mr_positive_mod(k, MR_M));
}

mr_uct_record_t mr_product_record_get(device const float *product, device const int *nbor, int subgrid_idx, int head_idx, int orientation, int normal, int t1, int t2) {
    int owner_subgrid;
    int face = mr_product_owner_face(nbor, subgrid_idx, orientation, normal, t1, t2, owner_subgrid);
    mr_uct_record_t record;
    for (int field = 0; field < MR_UCT_PRODUCT_RECORD_FIELDS; ++field) mr_record_set_value(record, field, mr_product_get(product, owner_subgrid, head_idx, MR_UCT_PRODUCT_FLUX_FIELDS + field, face));
    record.Bn = 0.0f;
    return record;
}

float mr_product_flux_get(device const float *product, device const int *nbor, int subgrid_idx, int head_idx, int orientation, int normal, int t1, int t2, int field) {
    int owner_subgrid;
    int face = mr_product_owner_face(nbor, subgrid_idx, orientation, normal, t1, t2, owner_subgrid);
    return mr_product_get(product, owner_subgrid, head_idx, field, face);
}

void mr_interior_record_store(threadgroup float *s, int face, mr_uct_record_t r) {
    for (int field = 0; field < 6; ++field) s[MR_FACE_BASE + MR_FLUX_FLOATS + field * MR_ALL_FACES + face] = mr_record_value(r, field);
}

mr_uct_record_t mr_interior_record_load(threadgroup const float *s, int face) {
    mr_uct_record_t r;
    for (int field = 0; field < 6; ++field) mr_record_set_value(r, field, s[MR_FACE_BASE + MR_FLUX_FLOATS + field * MR_ALL_FACES + face]);
    return r;
}

void mr_shell_record_store(threadgroup float *s, int shell, mr_uct_record_t r) {
    for (int field = 0; field < 6; ++field) {
        int index = field * MR_ALL_SHELLS + shell;
        if (index < MR_SHELL_FREE) s[MR_FACE_BASE + MR_FLUX_FLOATS + MR_RECORD_FLOATS + index] = mr_record_value(r, field);
        else s[index - MR_SHELL_FREE] = mr_record_value(r, field);
    }
}

mr_uct_record_t mr_shell_record_load(threadgroup const float *s, int shell) {
    mr_uct_record_t r;
    for (int field = 0; field < 6; ++field) {
        int index = field * MR_ALL_SHELLS + shell;
        float value = index < MR_SHELL_FREE ? s[MR_FACE_BASE + MR_FLUX_FLOATS + MR_RECORD_FLOATS + index] : s[index - MR_SHELL_FREE];
        mr_record_set_value(r, field, value);
    }
    return r;
}

int mr_emfz_index(int i, int j, int k) {
    return i + (MR_M + 1) * (j + (MR_M + 1) * k);
}

int mr_emfy_index(int i, int j, int k) {
    return MR_EDGE_CELLS + i + (MR_M + 1) * (j + MR_M * k);
}

int mr_emfx_index(int i, int j, int k) {
    return 2 * MR_EDGE_CELLS + i + MR_M * (j + (MR_M + 1) * k);
}

float mr_emf_get(threadgroup const float *s, int orientation, int i, int j, int k) {
    int index = orientation == 0 ? mr_emfz_index(i, j, k) : orientation == 1 ? mr_emfy_index(i, j, k) : mr_emfx_index(i, j, k);
    return s[MR_EMF_BASE + index];
}

void mr_emf_set(threadgroup float *s, int orientation, int i, int j, int k, float value) {
    int index = orientation == 0 ? mr_emfz_index(i, j, k) : orientation == 1 ? mr_emfy_index(i, j, k) : mr_emfx_index(i, j, k);
    s[MR_EMF_BASE + index] = value;
}

int mr_hlld_edge_index(int orientation, int i, int j, int k) {
    if (orientation == 0) return i + (MR_M + 1) * (j + (MR_M + 1) * k);
    if (orientation == 1) return i + (MR_M + 1) * (j + MR_M * k);
    return i + MR_M * (j + (MR_M + 1) * k);
}

void mr_hlld_corner_set(threadgroup float *s, int edge, int slot, primitive_t q) {
    int base = MR_FACE_BASE + (slot * MR_HLLD_CORNER_FIELDS) * MR_EDGE_CELLS + edge;
    s[base] = q.density;
    s[base + MR_EDGE_CELLS] = q.velocity_x;
    s[base + 2 * MR_EDGE_CELLS] = q.velocity_y;
    s[base + 3 * MR_EDGE_CELLS] = q.pressure;
    s[base + 4 * MR_EDGE_CELLS] = q.Bx;
    s[base + 5 * MR_EDGE_CELLS] = q.By;
    s[base + 6 * MR_EDGE_CELLS] = q.Bz;
}

primitive_t mr_hlld_corner_get(threadgroup const float *s, int edge, int slot) {
    int base = MR_FACE_BASE + (slot * MR_HLLD_CORNER_FIELDS) * MR_EDGE_CELLS + edge;
    primitive_t q;
    q.density = s[base];
    q.velocity_x = s[base + MR_EDGE_CELLS];
    q.velocity_y = s[base + 2 * MR_EDGE_CELLS];
    q.velocity_z = 0.0f;
    q.pressure = s[base + 3 * MR_EDGE_CELLS];
    q.Bx = s[base + 4 * MR_EDGE_CELLS];
    q.By = s[base + 5 * MR_EDGE_CELLS];
    q.Bz = s[base + 6 * MR_EDGE_CELLS];
    return q;
}

void mr_hlld_emf_set(threadgroup float *s, int orientation, int edge, float value) {
    s[MR_HLLD_EMF_BASE + orientation * MR_EDGE_CELLS + edge] = value;
}

float mr_emag(float x, float y, float z) {
    return 0.5f * magnitude_squared(x, y, z);
}

float mr_face_b_slope(float bm, float b0, float bp, int slope_mag) {
    if (slope_mag == 0) return 0.0f;
    float factor = float(min(slope_mag, 2));
    float dl = factor * (b0 - bm);
    float dr = factor * (bp - b0);
    float dc = 0.5f * (dl + dr) / factor;
    float limit = dl * dr <= 0.0f ? 0.0f : min(abs(dl), abs(dr));
    return copysign(min(limit, abs(dc)), dc);
}

float mr_hll(float sl, float sr, float fl, float fr, float ul, float ur) {
    if (sl >= 0.0f) return fl;
    if (sr <= 0.0f) return fr;
    return (sr * fl - sl * fr + sl * sr * (ur - ul)) / (sr - sl);
}

primitive_t mr_rotate_face(primitive_t q, int orientation) {
    if (orientation == 0) return q;
    primitive_t r = q;
    if (orientation == 1) {
        r.velocity_x = q.velocity_y;
        r.velocity_y = q.velocity_x;
        r.velocity_z = q.velocity_z;
        r.Bx = q.By;
        r.By = q.Bx;
        r.Bz = q.Bz;
    } else {
        r.velocity_x = q.velocity_z;
        r.velocity_y = q.velocity_x;
        r.velocity_z = q.velocity_y;
        r.Bx = q.Bz;
        r.By = q.Bx;
        r.Bz = q.By;
    }
    return r;
}

conserved_t mr_unrotate_flux(conserved_t f, int orientation) {
    if (orientation == 0) return f;
    conserved_t r = f;
    if (orientation == 1) {
        r.momentum_x = f.momentum_y;
        r.momentum_y = f.momentum_x;
        r.momentum_z = f.momentum_z;
    } else {
        r.momentum_x = f.momentum_y;
        r.momentum_y = f.momentum_z;
        r.momentum_z = f.momentum_x;
    }
    return r;
}

void mr_set_uct_record(thread mr_uct_record_t &record, float ul, float ur, float vl, float vr, float wl, float wr, float A, float SL, float SR, float ustar, float SAL, float SAR) {
    float ap = max(SR, 0.0f);
    float am = -min(SL, 0.0f);
    float asum = ap + am;
    record.vt1 = asum > 0.0f ? (ap * vl + am * vr) / asum : 0.5f * (vl + vr);
    record.vt2 = asum > 0.0f ? (ap * wl + am * wr) / asum : 0.5f * (wl + wr);
    float eps = max(1.0e-9f, 16.0f * FLT_EPSILON);
    float den = abs(SAR) + abs(SAL);
    float nustar = abs(SAR - SAL) > eps * abs(SR - SL) && den > 0.0f ? (SAR + SAL) / den : 0.0f;
    den = abs(SAL) + abs(SL);
    float nuL = den > 0.0f ? (SAL + SL) / den : 0.0f;
    den = abs(SAR) + abs(SR);
    float nuR = den > 0.0f ? (SAR + SR) / den : 0.0f;
    den = SAL + SL - 2.0f * ustar;
    float scale = max(abs(SAL), max(abs(SL), abs(ustar)));
    float tchiL = abs(den) > eps * scale ? (ul - ustar) * (SL - ustar) / den : 0.5f * (ul - ustar);
    den = SAR + SR - 2.0f * ustar;
    scale = max(abs(SAR), max(abs(SR), abs(ustar)));
    float tchiR = abs(den) > eps * scale ? (ur - ustar) * (SR - ustar) / den : 0.5f * (ur - ustar);
    record.aL = 0.5f * (1.0f + nustar);
    record.dL = 0.5f * (nuL - nustar) * tchiL + 0.5f * (abs(SAL) - nustar * SAL);
    record.dR = 0.5f * (nuR - nustar) * tchiR + 0.5f * (abs(SAR) - nustar * SAR);
    record.Bn = A;
}

void mr_set_uct_hll_record(thread mr_uct_record_t &record, float vl, float vr, float wl, float wr, float A, float SL, float SR) {
    float den = SR - SL;
    record.aL = SR / den;
    record.dL = -SL * SR / den;
    record.dR = record.dL;
    record.vt1 = (SR * vl - SL * vr) / den;
    record.vt2 = (SR * wl - SL * wr) / den;
    record.Bn = A;
}

bool mr_finite(float x) {
    return fabs(x) <= FLT_MAX;
}

bool mr_uct_degenerate(float SL, float SR, float ustar, float SAL, float SAR) {
    return SL < 0.0f && SR > 0.0f && (SAL - SL < 1.0e-4f * (ustar - SL) || SAR - SR > -1.0e-4f * (SR - ustar));
}

bool mr_uct_record_nonfinite(mr_uct_record_t r, float SL, float SR, float ustar, float SAL, float SAR) {
    return !(mr_finite(SL) && mr_finite(SR) && mr_finite(ustar) && mr_finite(SAL) && mr_finite(SAR) &&
             mr_finite(r.aL) && mr_finite(r.dL) && mr_finite(r.dR) &&
             mr_finite(r.vt1) && mr_finite(r.vt2) && mr_finite(r.Bn));
}

mr_uct_record_t mr_llf_record_lmax(primitive_t left, primitive_t right, float A, float lmax) {
    mr_uct_record_t record;
    record.aL = 0.5f;
    record.dL = 0.5f * lmax;
    record.dR = 0.5f * lmax;
    record.vt1 = 0.5f * (left.velocity_y + right.velocity_y);
    record.vt2 = 0.5f * (left.velocity_z + right.velocity_z);
    record.Bn = A;
    return record;
}

float mr_llf_signal_speed(primitive_t left, primitive_t right, float gamma, float smallr, float smallc2) {
    float smallp = smallc2 / gamma;
    left.density = max(left.density, smallr);
    right.density = max(right.density, smallr);
    left.pressure = max(left.pressure, smallp * left.density);
    right.pressure = max(right.pressure, smallp * right.density);
    float A = 0.5f * (left.Bx + right.Bx);
    float b2l = A * A + left.By * left.By + left.Bz * left.Bz;
    float c2l = gamma * left.pressure / left.density;
    float h2l = 0.5f * (b2l / left.density + c2l);
    float cfl = sqrt(h2l + sqrt(max(h2l * h2l - c2l * A * A / left.density, 0.0f)));
    float b2r = A * A + right.By * right.By + right.Bz * right.Bz;
    float c2r = gamma * right.pressure / right.density;
    float h2r = 0.5f * (b2r / right.density + c2r);
    float cfr = sqrt(h2r + sqrt(max(h2r * h2r - c2r * A * A / right.density, 0.0f)));
    return max(abs(left.velocity_x) + cfl, abs(right.velocity_x) + cfr);
}

mr_uct_record_t mr_llf_record(primitive_t left, primitive_t right, float gamma, float smallr, float smallc2) {
    float A = 0.5f * (left.Bx + right.Bx);
    float lmax = mr_llf_signal_speed(left, right, gamma, smallr, smallc2);
    return mr_llf_record_lmax(left, right, A, lmax);
}

conserved_t mr_hll_mhd_flux_lmax(primitive_t left, primitive_t right, float gamma, float smallr, float smallc2, thread float &lmax_out) {
    float smallp = smallc2 / gamma;
    float A = 0.5f * (left.Bx + right.Bx);
    lmax_out = mr_llf_signal_speed(left, right, gamma, smallr, smallc2);
    left.Bx = A;
    right.Bx = A;
    left.density = max(left.density, smallr);
    right.density = max(right.density, smallr);
    left.pressure = max(left.pressure, smallp * left.density);
    right.pressure = max(right.pressure, smallp * right.density);
    conserved_t ul = primitive_2_conserved(left, gamma);
    conserved_t ur = primitive_2_conserved(right, gamma);
    float ptl = left.pressure + mr_emag(left.Bx, left.By, left.Bz);
    float ptr = right.pressure + mr_emag(right.Bx, right.By, right.Bz);
    float vbl = left.Bx * left.velocity_x + left.By * left.velocity_y + left.Bz * left.velocity_z;
    float vbr = right.Bx * right.velocity_x + right.By * right.velocity_y + right.Bz * right.velocity_z;
    conserved_t fl;
    fl.density = left.density * left.velocity_x;
    fl.momentum_x = left.density * left.velocity_x * left.velocity_x + ptl - A * A;
    fl.momentum_y = left.density * left.velocity_x * left.velocity_y - A * left.By;
    fl.momentum_z = left.density * left.velocity_x * left.velocity_z - A * left.Bz;
    fl.energy = (ul.energy + ptl) * left.velocity_x - A * vbl;
    fl.Bx = 0.0f;
    fl.By = left.By * left.velocity_x - A * left.velocity_y;
    fl.Bz = left.Bz * left.velocity_x - A * left.velocity_z;
    conserved_t fr;
    fr.density = right.density * right.velocity_x;
    fr.momentum_x = right.density * right.velocity_x * right.velocity_x + ptr - A * A;
    fr.momentum_y = right.density * right.velocity_x * right.velocity_y - A * right.By;
    fr.momentum_z = right.density * right.velocity_x * right.velocity_z - A * right.Bz;
    fr.energy = (ur.energy + ptr) * right.velocity_x - A * vbr;
    fr.Bx = 0.0f;
    fr.By = right.By * right.velocity_x - A * right.velocity_y;
    fr.Bz = right.Bz * right.velocity_x - A * right.velocity_z;
    float sl = -lmax_out;
    float sr = lmax_out;
    conserved_t flux;
    flux.density = mr_hll(sl, sr, fl.density, fr.density, ul.density, ur.density);
    flux.momentum_x = mr_hll(sl, sr, fl.momentum_x, fr.momentum_x, ul.momentum_x, ur.momentum_x);
    flux.momentum_y = mr_hll(sl, sr, fl.momentum_y, fr.momentum_y, ul.momentum_y, ur.momentum_y);
    flux.momentum_z = mr_hll(sl, sr, fl.momentum_z, fr.momentum_z, ul.momentum_z, ur.momentum_z);
    flux.energy = mr_hll(sl, sr, fl.energy, fr.energy, ul.energy, ur.energy);
    flux.Bx = mr_hll(sl, sr, fl.Bx, fr.Bx, ul.Bx, ur.Bx);
    flux.By = mr_hll(sl, sr, fl.By, fr.By, ul.By, ur.By);
    flux.Bz = mr_hll(sl, sr, fl.Bz, fr.Bz, ul.Bz, ur.Bz);
    return flux;
}

conserved_t mr_hll_mhd_flux(primitive_t left, primitive_t right, float gamma, float smallr, float smallc2) {
    float lmax;
    return mr_hll_mhd_flux_lmax(left, right, gamma, smallr, smallc2, lmax);
}

conserved_t mr_hlld_mhd_flux(primitive_t left, primitive_t right, float gamma, float smallr, float smallc2, thread mr_uct_record_t &uct, bool make_uct_record) {
    float smallp = smallc2 / gamma;
    float entho = 1.0f / (gamma - 1.0f);
    primitive_t raw_left = left;
    primitive_t raw_right = right;
    left.density = max(left.density, smallr);
    right.density = max(right.density, smallr);
    left.pressure = max(left.pressure, smallp * left.density);
    right.pressure = max(right.pressure, smallp * right.density);
    float A = 0.5f * (left.Bx + right.Bx);
    float sgnm = copysign(1.0f, A);
    float rl = left.density;
    float ul = left.velocity_x;
    float vl = left.velocity_y;
    float wl = left.velocity_z;
    float Pl = left.pressure;
    float Bl = left.By;
    float Cl = left.Bz;
    float ekinl = 0.5f * (ul * ul + vl * vl + wl * wl) * rl;
    float emagl = 0.5f * (A * A + Bl * Bl + Cl * Cl);
    float etotl = Pl * entho + ekinl + emagl;
    float Ptotl = Pl + emagl;
    float vdotBl = ul * A + vl * Bl + wl * Cl;
    float rr = right.density;
    float ur = right.velocity_x;
    float vr = right.velocity_y;
    float wr = right.velocity_z;
    float Pr = right.pressure;
    float Br = right.By;
    float Cr = right.Bz;
    float ekinr = 0.5f * (ur * ur + vr * vr + wr * wr) * rr;
    float emagr = 0.5f * (A * A + Br * Br + Cr * Cr);
    float etotr = Pr * entho + ekinr + emagr;
    float Ptotr = Pr + emagr;
    float vdotBr = ur * A + vr * Br + wr * Cr;
    float b2l = A * A + Bl * Bl + Cl * Cl;
    float c2l = gamma * Pl / rl;
    float d2l = 0.5f * (b2l / rl + c2l);
    float cfastl = sqrt(d2l + sqrt(max(d2l * d2l - c2l * A * A / rl, 0.0f)));
    float b2r = A * A + Br * Br + Cr * Cr;
    float c2r = gamma * Pr / rr;
    float d2r = 0.5f * (b2r / rr + c2r);
    float cfastr = sqrt(d2r + sqrt(max(d2r * d2r - c2r * A * A / rr, 0.0f)));
    float SL = min(ul, ur) - max(cfastl, cfastr);
    float SR = max(ul, ur) + max(cfastl, cfastr);
    float rcl = rl * (ul - SL);
    float rcr = rr * (SR - ur);
    float ustar = (rcr * ur + rcl * ul + Ptotl - Ptotr) / (rcr + rcl);
    float Ptotstar = (rcr * Ptotl + rcl * Ptotr + rcl * rcr * (ul - ur)) / (rcr + rcl);
    float rstarl = max(rl * (SL - ul) / (SL - ustar), smallr);
    float estar = rl * (SL - ul) * (SL - ustar) - A * A;
    float el = rl * (SL - ul) * (SL - ul) - A * A;
    float vstarl;
    float Bstarl;
    float wstarl;
    float Cstarl;
    if (abs(estar) < 1.0e-4f * A * A) {
        vstarl = vl;
        Bstarl = Bl;
        wstarl = wl;
        Cstarl = Cl;
    } else {
        vstarl = vl - A * Bl * (ustar - ul) / estar;
        Bstarl = Bl * el / estar;
        wstarl = wl - A * Cl * (ustar - ul) / estar;
        Cstarl = Cl * el / estar;
    }
    float vdotBstarl = ustar * A + vstarl * Bstarl + wstarl * Cstarl;
    float etotstarl = ((SL - ul) * etotl - Ptotl * ul + Ptotstar * ustar + A * (vdotBl - vdotBstarl)) / (SL - ustar);
    float sqrrstarl = sqrt(rstarl);
    float SAL = ustar - abs(A) / sqrrstarl;
    float rstarr = max(rr * (SR - ur) / (SR - ustar), smallr);
    estar = rr * (SR - ur) * (SR - ustar) - A * A;
    float er = rr * (SR - ur) * (SR - ur) - A * A;
    float vstarr;
    float Bstarr;
    float wstarr;
    float Cstarr;
    if (abs(estar) < 1.0e-4f * A * A) {
        vstarr = vr;
        Bstarr = Br;
        wstarr = wr;
        Cstarr = Cr;
    } else {
        vstarr = vr - A * Br * (ustar - ur) / estar;
        Bstarr = Br * er / estar;
        wstarr = wr - A * Cr * (ustar - ur) / estar;
        Cstarr = Cr * er / estar;
    }
    float vdotBstarr = ustar * A + vstarr * Bstarr + wstarr * Cstarr;
    float etotstarr = ((SR - ur) * etotr - Ptotr * ur + Ptotstar * ustar + A * (vdotBr - vdotBstarr)) / (SR - ustar);
    float sqrrstarr = sqrt(rstarr);
    float SAR = ustar + abs(A) / sqrrstarr;
    // UCT needs a nonsingular paired record and flux when the rotational fan collapses.
    bool use_hllc = make_uct_record && mr_uct_degenerate(SL, SR, ustar, SAL, SAR);
    if (make_uct_record) {
        if (use_hllc) mr_set_uct_hll_record(uct, vl, vr, wl, wr, A, SL, SR);
        else mr_set_uct_record(uct, ul, ur, vl, vr, wl, wr, A, SL, SR, ustar, SAL, SAR);
        if (mr_uct_record_nonfinite(uct, SL, SR, ustar, SAL, SAR)) {
            float llf_lmax;
            conserved_t llf_flux = mr_hll_mhd_flux_lmax(raw_left, raw_right, gamma, smallr, smallc2, llf_lmax);
            uct = mr_llf_record_lmax(left, right, A, llf_lmax);
            return llf_flux;
        }
    }
    if (use_hllc) {
        float inv_span = 1.0f / (SR - SL);
        Bstarl = Bstarr = (SR * Br - SL * Bl + Bl * ul - A * vl - Br * ur + A * vr) * inv_span;
        Cstarl = Cstarr = (SR * Cr - SL * Cl + Cl * ul - A * wl - Cr * ur + A * wr) * inv_span;
        vstarl = vl - A * (Bstarl - Bl) / (rl * (SL - ul));
        wstarl = wl - A * (Cstarl - Cl) / (rl * (SL - ul));
        vstarr = vr - A * (Bstarr - Br) / (rr * (SR - ur));
        wstarr = wr - A * (Cstarr - Cr) / (rr * (SR - ur));
        vdotBstarl = ustar * A + vstarl * Bstarl + wstarl * Cstarl;
        vdotBstarr = ustar * A + vstarr * Bstarr + wstarr * Cstarr;
        etotstarl = ((SL - ul) * etotl - Ptotl * ul + Ptotstar * ustar + A * (vdotBl - vdotBstarl)) / (SL - ustar);
        etotstarr = ((SR - ur) * etotr - Ptotr * ur + Ptotstar * ustar + A * (vdotBr - vdotBstarr)) / (SR - ustar);
        SAL = ustar;
        SAR = ustar;
    }
    float denom = sqrrstarl + sqrrstarr;
    float vstarstar = (sqrrstarl * vstarl + sqrrstarr * vstarr + sgnm * (Bstarr - Bstarl)) / denom;
    float wstarstar = (sqrrstarl * wstarl + sqrrstarr * wstarr + sgnm * (Cstarr - Cstarl)) / denom;
    float Bstarstar = (sqrrstarl * Bstarr + sqrrstarr * Bstarl + sgnm * sqrrstarl * sqrrstarr * (vstarr - vstarl)) / denom;
    float Cstarstar = (sqrrstarl * Cstarr + sqrrstarr * Cstarl + sgnm * sqrrstarl * sqrrstarr * (wstarr - wstarl)) / denom;
    float vdotBstarstar = ustar * A + vstarstar * Bstarstar + wstarstar * Cstarstar;
    float etotstarstarl = etotstarl - sgnm * sqrrstarl * (vdotBstarl - vdotBstarstar);
    float etotstarstarr = etotstarr + sgnm * sqrrstarr * (vdotBstarr - vdotBstarstar);
    float ro;
    float uo;
    float vo;
    float wo;
    float Bo;
    float Co;
    float Ptoto;
    float etoto;
    float vdotBo;
    if (use_hllc && ustar >= 0.0f) {
        ro = rstarl; uo = ustar; vo = vstarl; wo = wstarl; Bo = Bstarl; Co = Cstarl; Ptoto = Ptotstar; etoto = etotstarl; vdotBo = vdotBstarl;
    } else if (use_hllc) {
        ro = rstarr; uo = ustar; vo = vstarr; wo = wstarr; Bo = Bstarr; Co = Cstarr; Ptoto = Ptotstar; etoto = etotstarr; vdotBo = vdotBstarr;
    } else if (SL > 0.0f) {
        ro = rl; uo = ul; vo = vl; wo = wl; Bo = Bl; Co = Cl; Ptoto = Ptotl; etoto = etotl; vdotBo = vdotBl;
    } else if (SAL > 0.0f) {
        ro = rstarl; uo = ustar; vo = vstarl; wo = wstarl; Bo = Bstarl; Co = Cstarl; Ptoto = Ptotstar; etoto = etotstarl; vdotBo = vdotBstarl;
    } else if (ustar > 0.0f) {
        ro = rstarl; uo = ustar; vo = vstarstar; wo = wstarstar; Bo = Bstarstar; Co = Cstarstar; Ptoto = Ptotstar; etoto = etotstarstarl; vdotBo = vdotBstarstar;
    } else if (SAR > 0.0f) {
        ro = rstarr; uo = ustar; vo = vstarstar; wo = wstarstar; Bo = Bstarstar; Co = Cstarstar; Ptoto = Ptotstar; etoto = etotstarstarr; vdotBo = vdotBstarstar;
    } else if (SR > 0.0f) {
        ro = rstarr; uo = ustar; vo = vstarr; wo = wstarr; Bo = Bstarr; Co = Cstarr; Ptoto = Ptotstar; etoto = etotstarr; vdotBo = vdotBstarr;
    } else {
        ro = rr; uo = ur; vo = vr; wo = wr; Bo = Br; Co = Cr; Ptoto = Ptotr; etoto = etotr; vdotBo = vdotBr;
    }
    conserved_t flux;
    flux.density = ro * uo;
    flux.momentum_x = ro * uo * uo + Ptoto - A * A;
    flux.momentum_y = ro * uo * vo - A * Bo;
    flux.momentum_z = ro * uo * wo - A * Co;
    flux.energy = (etoto + Ptoto) * uo - A * vdotBo;
    flux.Bx = 0.0f;
    flux.By = Bo * uo - A * vo;
    flux.Bz = Co * uo - A * wo;
    return flux;
}

float mr_hlld_2d_emf(primitive_t stateLL, primitive_t stateLR, primitive_t stateRL, primitive_t stateRR, float gamma, float smallr, float smallc, float switch_llf_dmin, float switch_llf_pmin) {
    float smallp = smallc * smallc / gamma;
    stateLL.density = max(stateLL.density, smallr);
    stateLR.density = max(stateLR.density, smallr);
    stateRL.density = max(stateRL.density, smallr);
    stateRR.density = max(stateRR.density, smallr);
    stateLL.pressure = max(stateLL.pressure, smallp * stateLL.density);
    stateLR.pressure = max(stateLR.pressure, smallp * stateLR.density);
    stateRL.pressure = max(stateRL.pressure, smallp * stateRL.density);
    stateRR.pressure = max(stateRR.pressure, smallp * stateRR.density);

    float rLL = stateLL.density;
    float rLR = stateLR.density;
    float rRL = stateRL.density;
    float rRR = stateRR.density;
    float pLL = stateLL.pressure;
    float pLR = stateLR.pressure;
    float pRL = stateRL.pressure;
    float pRR = stateRR.pressure;
    float uLL = stateLL.velocity_x;
    float uLR = stateLR.velocity_x;
    float uRL = stateRL.velocity_x;
    float uRR = stateRR.velocity_x;
    float vLL = stateLL.velocity_y;
    float vLR = stateLR.velocity_y;
    float vRL = stateRL.velocity_y;
    float vRR = stateRR.velocity_y;
    float ALL = stateLL.Bx;
    float ALR = stateLR.Bx;
    float ARL = stateRL.Bx;
    float ARR = stateRR.Bx;
    float BLL = stateLL.By;
    float BLR = stateLR.By;
    float BRL = stateRL.By;
    float BRR = stateRR.By;
    float CLL = stateLL.Bz;
    float CLR = stateLR.Bz;
    float CRL = stateRL.Bz;
    float CRR = stateRR.Bz;

    float b2LL = magnitude_squared(ALL, BLL, CLL);
    float c2LL = gamma * pLL / rLL;
    float d2LL = 0.5f * (b2LL / rLL + c2LL);
    float cfastLLx = sqrt(d2LL + sqrt(max(d2LL * d2LL - c2LL * ALL * ALL / rLL, 0.0f)));
    float cfastLLy = sqrt(d2LL + sqrt(max(d2LL * d2LL - c2LL * BLL * BLL / rLL, 0.0f)));
    float b2LR = magnitude_squared(ALR, BLR, CLR);
    float c2LR = gamma * pLR / rLR;
    float d2LR = 0.5f * (b2LR / rLR + c2LR);
    float cfastLRx = sqrt(d2LR + sqrt(max(d2LR * d2LR - c2LR * ALR * ALR / rLR, 0.0f)));
    float cfastLRy = sqrt(d2LR + sqrt(max(d2LR * d2LR - c2LR * BLR * BLR / rLR, 0.0f)));
    float b2RL = magnitude_squared(ARL, BRL, CRL);
    float c2RL = gamma * pRL / rRL;
    float d2RL = 0.5f * (b2RL / rRL + c2RL);
    float cfastRLx = sqrt(d2RL + sqrt(max(d2RL * d2RL - c2RL * ARL * ARL / rRL, 0.0f)));
    float cfastRLy = sqrt(d2RL + sqrt(max(d2RL * d2RL - c2RL * BRL * BRL / rRL, 0.0f)));
    float b2RR = magnitude_squared(ARR, BRR, CRR);
    float c2RR = gamma * pRR / rRR;
    float d2RR = 0.5f * (b2RR / rRR + c2RR);
    float cfastRRx = sqrt(d2RR + sqrt(max(d2RR * d2RR - c2RR * ARR * ARR / rRR, 0.0f)));
    float cfastRRy = sqrt(d2RR + sqrt(max(d2RR * d2RR - c2RR * BRR * BRR / rRR, 0.0f)));

    float umin = min(min(uLL, uLR), min(uRL, uRR));
    float umax = max(max(uLL, uLR), max(uRL, uRR));
    float vmin = min(min(vLL, vLR), min(vRL, vRR));
    float vmax = max(max(vLL, vLR), max(vRL, vRR));
    float cfastx = max(max(cfastLLx, cfastLRx), max(cfastRLx, cfastRRx));
    float cfasty = max(max(cfastLLy, cfastLRy), max(cfastRLy, cfastRRy));
    float SL = umin - cfastx;
    float SR = umax + cfastx;
    float SB = vmin - cfasty;
    float ST = vmax + cfasty;

    float ELL = uLL * BLL - vLL * ALL;
    float ELR = uLR * BLR - vLR * ALR;
    float ERL = uRL * BRL - vRL * ARL;
    float ERR = uRR * BRR - vRR * ARR;
    float rmin = min(min(rLL, rLR), min(rRL, rRR));
    float pmin = min(min(pLL, pLR), min(pRL, pRR));
    bool switch_to_llf = mr_switch_to_llf(rmin, pmin, switch_llf_dmin, switch_llf_pmin);
    if (switch_to_llf) {
        float Smax = max(max(abs(SR), abs(ST)), max(abs(SL), abs(SB)));
        return 0.25f * (ERR + ERL + ELR + ELL) + 0.5f * Smax * (stateRR.Bx - stateLL.Bx) - 0.5f * Smax * (stateRR.By - stateLL.By);
    }

    float PtotLL = pLL + 0.5f * b2LL;
    float PtotLR = pLR + 0.5f * b2LR;
    float PtotRL = pRL + 0.5f * b2RL;
    float PtotRR = pRR + 0.5f * b2RR;
    float rcLLx = rLL * (uLL - SL);
    float rcRLx = rRL * (SR - uRL);
    float rcLRx = rLR * (uLR - SL);
    float rcRRx = rRR * (SR - uRR);
    float rcLLy = rLL * (vLL - SB);
    float rcLRy = rLR * (ST - vLR);
    float rcRLy = rRL * (vRL - SB);
    float rcRRy = rRR * (ST - vRR);
    float ustar = (rcLLx * uLL + rcLRx * uLR + rcRLx * uRL + rcRRx * uRR + PtotLL - PtotRL + PtotLR - PtotRR) / (rcLLx + rcLRx + rcRLx + rcRRx);
    float vstar = (rcLLy * vLL + rcLRy * vLR + rcRLy * vRL + rcRRy * vRR + PtotLL - PtotLR + PtotRL - PtotRR) / (rcLLy + rcLRy + rcRLy + rcRRy);

    float rstarLLx = rLL * (SL - uLL) / (SL - ustar);
    float BstarLL = BLL * (SL - uLL) / (SL - ustar);
    float rstarLLy = rLL * (SB - vLL) / (SB - vstar);
    float AstarLL = ALL * (SB - vLL) / (SB - vstar);
    float rstarLL = rLL * (SL - uLL) / (SL - ustar) * (SB - vLL) / (SB - vstar);
    float EstarLLx = ustar * BstarLL - vLL * ALL;
    float EstarLLy = uLL * BLL - vstar * AstarLL;
    float EstarLL = ustar * BstarLL - vstar * AstarLL;

    float rstarLRx = rLR * (SL - uLR) / (SL - ustar);
    float BstarLR = BLR * (SL - uLR) / (SL - ustar);
    float rstarLRy = rLR * (ST - vLR) / (ST - vstar);
    float AstarLR = ALR * (ST - vLR) / (ST - vstar);
    float rstarLR = rLR * (SL - uLR) / (SL - ustar) * (ST - vLR) / (ST - vstar);
    float EstarLRx = ustar * BstarLR - vLR * ALR;
    float EstarLRy = uLR * BLR - vstar * AstarLR;
    float EstarLR = ustar * BstarLR - vstar * AstarLR;

    float rstarRLx = rRL * (SR - uRL) / (SR - ustar);
    float BstarRL = BRL * (SR - uRL) / (SR - ustar);
    float rstarRLy = rRL * (SB - vRL) / (SB - vstar);
    float AstarRL = ARL * (SB - vRL) / (SB - vstar);
    float rstarRL = rRL * (SR - uRL) / (SR - ustar) * (SB - vRL) / (SB - vstar);
    float EstarRLx = ustar * BstarRL - vRL * ARL;
    float EstarRLy = uRL * BRL - vstar * AstarRL;
    float EstarRL = ustar * BstarRL - vstar * AstarRL;

    float rstarRRx = rRR * (SR - uRR) / (SR - ustar);
    float BstarRR = BRR * (SR - uRR) / (SR - ustar);
    float rstarRRy = rRR * (ST - vRR) / (ST - vstar);
    float AstarRR = ARR * (ST - vRR) / (ST - vstar);
    float rstarRR = rRR * (SR - uRR) / (SR - ustar) * (ST - vRR) / (ST - vstar);
    float EstarRRx = ustar * BstarRR - vRR * ARR;
    float EstarRRy = uRR * BRR - vstar * AstarRR;
    float EstarRR = ustar * BstarRR - vstar * AstarRR;

    rstarLLx = max(rstarLLx, smallr);
    rstarLRx = max(rstarLRx, smallr);
    rstarRLx = max(rstarRLx, smallr);
    rstarRRx = max(rstarRRx, smallr);
    rstarLLy = max(rstarLLy, smallr);
    rstarLRy = max(rstarLRy, smallr);
    rstarRLy = max(rstarRLy, smallr);
    rstarRRy = max(rstarRRy, smallr);
    rstarLL = max(rstarLL, smallr);
    rstarLR = max(rstarLR, smallr);
    rstarRL = max(rstarRL, smallr);
    rstarRR = max(rstarRR, smallr);
    float calfvenL = max(max(abs(ALR) / sqrt(rstarLRx), abs(AstarLR) / sqrt(rstarLR)), max(max(abs(ALL) / sqrt(rstarLLx), abs(AstarLL) / sqrt(rstarLL)), smallc));
    float calfvenR = max(max(abs(ARR) / sqrt(rstarRRx), abs(AstarRR) / sqrt(rstarRR)), max(max(abs(ARL) / sqrt(rstarRLx), abs(AstarRL) / sqrt(rstarRL)), smallc));
    float calfvenB = max(max(abs(BLL) / sqrt(rstarLLy), abs(BstarLL) / sqrt(rstarLL)), max(max(abs(BRL) / sqrt(rstarRLy), abs(BstarRL) / sqrt(rstarRL)), smallc));
    float calfvenT = max(max(abs(BLR) / sqrt(rstarLRy), abs(BstarLR) / sqrt(rstarLR)), max(max(abs(BRR) / sqrt(rstarRRy), abs(BstarRR) / sqrt(rstarRR)), smallc));
    float SAL = min(ustar - calfvenL, 0.0f);
    float SAR = max(ustar + calfvenR, 0.0f);
    float SAB = min(vstar - calfvenB, 0.0f);
    float SAT = max(vstar + calfvenT, 0.0f);
    float AstarT = (SAR * AstarRR - SAL * AstarLR) / (SAR - SAL);
    float AstarB = (SAR * AstarRL - SAL * AstarLL) / (SAR - SAL);
    float BstarR = (SAT * BstarRR - SAB * BstarRL) / (SAT - SAB);
    float BstarL = (SAT * BstarLR - SAB * BstarLL) / (SAT - SAB);

    if (SB > 0.0f) {
        if (SL > 0.0f) return ELL;
        if (SR < 0.0f) return ERL;
        return (SAR * EstarLLx - SAL * EstarRLx + SAR * SAL * (BRL - BLL)) / (SAR - SAL);
    }
    if (ST < 0.0f) {
        if (SL > 0.0f) return ELR;
        if (SR < 0.0f) return ERR;
        return (SAR * EstarLRx - SAL * EstarRRx + SAR * SAL * (BRR - BLR)) / (SAR - SAL);
    }
    if (SL > 0.0f) return (SAT * EstarLLy - SAB * EstarLRy - SAT * SAB * (ALR - ALL)) / (SAT - SAB);
    if (SR < 0.0f) return (SAT * EstarRLy - SAB * EstarRRy - SAT * SAB * (ARR - ARL)) / (SAT - SAB);
    return (SAL * SAB * EstarRR - SAL * SAT * EstarRL - SAR * SAB * EstarLR + SAR * SAT * EstarLL) / (SAR - SAL) / (SAT - SAB) - SAT * SAB / (SAT - SAB) * (AstarT - AstarB) + SAR * SAL / (SAR - SAL) * (BstarR - BstarL);
}

float mr_llf_2d_emf(primitive_t stateLL, primitive_t stateLR, primitive_t stateRL, primitive_t stateRR, float gamma, float smallr, float smallc2) {
    primitive_t xleft;
    primitive_t xright;
    xleft.density = 0.5f * (stateLL.density + stateLR.density);
    xright.density = 0.5f * (stateRR.density + stateRL.density);
    xleft.velocity_x = 0.5f * (stateLL.velocity_x + stateLR.velocity_x);
    xright.velocity_x = 0.5f * (stateRR.velocity_x + stateRL.velocity_x);
    xleft.velocity_y = 0.5f * (stateLL.velocity_y + stateLR.velocity_y);
    xright.velocity_y = 0.5f * (stateRR.velocity_y + stateRL.velocity_y);
    xleft.velocity_z = 0.5f * (stateLL.velocity_z + stateLR.velocity_z);
    xright.velocity_z = 0.5f * (stateRR.velocity_z + stateRL.velocity_z);
    xleft.pressure = 0.5f * (stateLL.pressure + stateLR.pressure);
    xright.pressure = 0.5f * (stateRR.pressure + stateRL.pressure);
    xleft.Bx = 0.5f * (stateLL.Bx + stateLR.Bx);
    xright.Bx = 0.5f * (stateRR.Bx + stateRL.Bx);
    xleft.By = 0.5f * (stateLL.By + stateLR.By);
    xright.By = 0.5f * (stateRR.By + stateRL.By);
    xleft.Bz = 0.5f * (stateLL.Bz + stateLR.Bz);
    xright.Bz = 0.5f * (stateRR.Bz + stateRL.Bz);
    float xdiff = -0.5f * mr_llf_signal_speed(xleft, xright, gamma, smallr, smallc2) * (xright.By - xleft.By);

    primitive_t yleft;
    primitive_t yright;
    yleft.density = 0.5f * (stateLL.density + stateRL.density);
    yright.density = 0.5f * (stateRR.density + stateLR.density);
    yleft.velocity_x = 0.5f * (stateLL.velocity_y + stateRL.velocity_y);
    yright.velocity_x = 0.5f * (stateRR.velocity_y + stateLR.velocity_y);
    yleft.velocity_y = 0.5f * (stateLL.velocity_x + stateRL.velocity_x);
    yright.velocity_y = 0.5f * (stateRR.velocity_x + stateLR.velocity_x);
    yleft.velocity_z = 0.5f * (stateLL.velocity_z + stateRL.velocity_z);
    yright.velocity_z = 0.5f * (stateRR.velocity_z + stateLR.velocity_z);
    yleft.pressure = 0.5f * (stateLL.pressure + stateRL.pressure);
    yright.pressure = 0.5f * (stateRR.pressure + stateLR.pressure);
    yleft.Bx = 0.5f * (stateLL.By + stateRL.By);
    yright.Bx = 0.5f * (stateRR.By + stateLR.By);
    yleft.By = 0.5f * (stateLL.Bx + stateRL.Bx);
    yright.By = 0.5f * (stateRR.Bx + stateLR.Bx);
    yleft.Bz = 0.5f * (stateLL.Bz + stateRL.Bz);
    yright.Bz = 0.5f * (stateRR.Bz + stateLR.Bz);
    float ydiff = -0.5f * mr_llf_signal_speed(yleft, yright, gamma, smallr, smallc2) * (yright.By - yleft.By);

    float ELL = stateLL.velocity_x * stateLL.By - stateLL.velocity_y * stateLL.Bx;
    float ERL = stateRL.velocity_x * stateRL.By - stateRL.velocity_y * stateRL.Bx;
    float ELR = stateLR.velocity_x * stateLR.By - stateLR.velocity_y * stateLR.Bx;
    float ERR = stateRR.velocity_x * stateRR.By - stateRR.velocity_y * stateRR.Bx;
    return 0.25f * (ELL + ERL + ELR + ERR) + xdiff - ydiff;
}

primitive_t mr_cell_load(threadgroup const float *s, int i, int j, int k) {
    primitive_t q;
    q.density = mr_local_get(s, 0, i, j, k);
    q.velocity_x = mr_local_get(s, 1, i, j, k);
    q.velocity_y = mr_local_get(s, 2, i, j, k);
    q.velocity_z = mr_local_get(s, 3, i, j, k);
    q.pressure = mr_local_get(s, 4, i, j, k);
    q.Bx = 0.5f * (mr_bf_get(s, 0, i, j, k) + mr_bf_get(s, 0, i + 1, j, k));
    q.By = 0.5f * (mr_bf_get(s, 1, i, j, k) + mr_bf_get(s, 1, i, j + 1, k));
    q.Bz = 0.5f * (mr_bf_get(s, 2, i, j, k) + mr_bf_get(s, 2, i, j, k + 1));
    return q;
}

void mr_compute_slopes(threadgroup const float *s, int i, int j, int k, int slope, thread primitive_t &sx, thread primitive_t &sy, thread primitive_t &sz) {
    primitive_t q = mr_cell_load(s, i, j, k);
    sx.density = 0.5f * slope_moncen(mr_local_get(s, 0, i - 1, j, k), q.density, mr_local_get(s, 0, i + 1, j, k), slope);
    sx.velocity_x = 0.5f * slope_moncen(mr_local_get(s, 1, i - 1, j, k), q.velocity_x, mr_local_get(s, 1, i + 1, j, k), slope);
    sx.velocity_y = 0.5f * slope_moncen(mr_local_get(s, 2, i - 1, j, k), q.velocity_y, mr_local_get(s, 2, i + 1, j, k), slope);
    sx.velocity_z = 0.5f * slope_moncen(mr_local_get(s, 3, i - 1, j, k), q.velocity_z, mr_local_get(s, 3, i + 1, j, k), slope);
    sx.pressure = 0.5f * slope_moncen(mr_local_get(s, 4, i - 1, j, k), q.pressure, mr_local_get(s, 4, i + 1, j, k), slope);
    float bym = 0.5f * (mr_bf_get(s, 1, i - 1, j, k) + mr_bf_get(s, 1, i - 1, j + 1, k));
    float byp = 0.5f * (mr_bf_get(s, 1, i + 1, j, k) + mr_bf_get(s, 1, i + 1, j + 1, k));
    float bzm = 0.5f * (mr_bf_get(s, 2, i - 1, j, k) + mr_bf_get(s, 2, i - 1, j, k + 1));
    float bzp = 0.5f * (mr_bf_get(s, 2, i + 1, j, k) + mr_bf_get(s, 2, i + 1, j, k + 1));
    sx.Bx = 0.0f;
    sx.By = 0.5f * slope_moncen(bym, q.By, byp, slope);
    sx.Bz = 0.5f * slope_moncen(bzm, q.Bz, bzp, slope);
    sy.density = 0.5f * slope_moncen(mr_local_get(s, 0, i, j - 1, k), q.density, mr_local_get(s, 0, i, j + 1, k), slope);
    sy.velocity_x = 0.5f * slope_moncen(mr_local_get(s, 1, i, j - 1, k), q.velocity_x, mr_local_get(s, 1, i, j + 1, k), slope);
    sy.velocity_y = 0.5f * slope_moncen(mr_local_get(s, 2, i, j - 1, k), q.velocity_y, mr_local_get(s, 2, i, j + 1, k), slope);
    sy.velocity_z = 0.5f * slope_moncen(mr_local_get(s, 3, i, j - 1, k), q.velocity_z, mr_local_get(s, 3, i, j + 1, k), slope);
    sy.pressure = 0.5f * slope_moncen(mr_local_get(s, 4, i, j - 1, k), q.pressure, mr_local_get(s, 4, i, j + 1, k), slope);
    float bxm = 0.5f * (mr_bf_get(s, 0, i, j - 1, k) + mr_bf_get(s, 0, i + 1, j - 1, k));
    float bxp = 0.5f * (mr_bf_get(s, 0, i, j + 1, k) + mr_bf_get(s, 0, i + 1, j + 1, k));
    bzm = 0.5f * (mr_bf_get(s, 2, i, j - 1, k) + mr_bf_get(s, 2, i, j - 1, k + 1));
    bzp = 0.5f * (mr_bf_get(s, 2, i, j + 1, k) + mr_bf_get(s, 2, i, j + 1, k + 1));
    sy.Bx = 0.5f * slope_moncen(bxm, q.Bx, bxp, slope);
    sy.By = 0.0f;
    sy.Bz = 0.5f * slope_moncen(bzm, q.Bz, bzp, slope);
    sz.density = 0.5f * slope_moncen(mr_local_get(s, 0, i, j, k - 1), q.density, mr_local_get(s, 0, i, j, k + 1), slope);
    sz.velocity_x = 0.5f * slope_moncen(mr_local_get(s, 1, i, j, k - 1), q.velocity_x, mr_local_get(s, 1, i, j, k + 1), slope);
    sz.velocity_y = 0.5f * slope_moncen(mr_local_get(s, 2, i, j, k - 1), q.velocity_y, mr_local_get(s, 2, i, j, k + 1), slope);
    sz.velocity_z = 0.5f * slope_moncen(mr_local_get(s, 3, i, j, k - 1), q.velocity_z, mr_local_get(s, 3, i, j, k + 1), slope);
    sz.pressure = 0.5f * slope_moncen(mr_local_get(s, 4, i, j, k - 1), q.pressure, mr_local_get(s, 4, i, j, k + 1), slope);
    bxm = 0.5f * (mr_bf_get(s, 0, i, j, k - 1) + mr_bf_get(s, 0, i + 1, j, k - 1));
    bxp = 0.5f * (mr_bf_get(s, 0, i, j, k + 1) + mr_bf_get(s, 0, i + 1, j, k + 1));
    bym = 0.5f * (mr_bf_get(s, 1, i, j, k - 1) + mr_bf_get(s, 1, i, j + 1, k - 1));
    byp = 0.5f * (mr_bf_get(s, 1, i, j, k + 1) + mr_bf_get(s, 1, i, j + 1, k + 1));
    sz.Bx = 0.5f * slope_moncen(bxm, q.Bx, bxp, slope);
    sz.By = 0.5f * slope_moncen(bym, q.By, byp, slope);
    sz.Bz = 0.0f;
}

float mr_edge_x(threadgroup const float *s, int i, int j, int k) {
    float v = 0.25f * (mr_local_get(s, 2, i, j - 1, k - 1) + mr_local_get(s, 2, i, j - 1, k) + mr_local_get(s, 2, i, j, k - 1) + mr_local_get(s, 2, i, j, k));
    float w = 0.25f * (mr_local_get(s, 3, i, j - 1, k - 1) + mr_local_get(s, 3, i, j - 1, k) + mr_local_get(s, 3, i, j, k - 1) + mr_local_get(s, 3, i, j, k));
    float B = 0.5f * (mr_bf_get(s, 1, i, j, k - 1) + mr_bf_get(s, 1, i, j, k));
    float C = 0.5f * (mr_bf_get(s, 2, i, j - 1, k) + mr_bf_get(s, 2, i, j, k));
    return v * C - w * B;
}

float mr_edge_y(threadgroup const float *s, int i, int j, int k) {
    float u = 0.25f * (mr_local_get(s, 1, i - 1, j, k - 1) + mr_local_get(s, 1, i - 1, j, k) + mr_local_get(s, 1, i, j, k - 1) + mr_local_get(s, 1, i, j, k));
    float w = 0.25f * (mr_local_get(s, 3, i - 1, j, k - 1) + mr_local_get(s, 3, i - 1, j, k) + mr_local_get(s, 3, i, j, k - 1) + mr_local_get(s, 3, i, j, k));
    float A = 0.5f * (mr_bf_get(s, 0, i, j, k - 1) + mr_bf_get(s, 0, i, j, k));
    float C = 0.5f * (mr_bf_get(s, 2, i - 1, j, k) + mr_bf_get(s, 2, i, j, k));
    return w * A - u * C;
}

float mr_edge_z(threadgroup const float *s, int i, int j, int k) {
    float u = 0.25f * (mr_local_get(s, 1, i - 1, j - 1, k) + mr_local_get(s, 1, i - 1, j, k) + mr_local_get(s, 1, i, j - 1, k) + mr_local_get(s, 1, i, j, k));
    float v = 0.25f * (mr_local_get(s, 2, i - 1, j - 1, k) + mr_local_get(s, 2, i - 1, j, k) + mr_local_get(s, 2, i, j - 1, k) + mr_local_get(s, 2, i, j, k));
    float A = 0.5f * (mr_bf_get(s, 0, i, j - 1, k) + mr_bf_get(s, 0, i, j, k));
    float B = 0.5f * (mr_bf_get(s, 1, i - 1, j, k) + mr_bf_get(s, 1, i, j, k));
    return u * B - v * A;
}

mr_trace_t mr_predict_cell(threadgroup const float *s, int i, int j, int k, float gamma, float dtdx, int slope, int slope_mag, int induction) {
    mr_trace_t t;
    t.AL = mr_bf_get(s, 0, i, j, k);
    t.AR = mr_bf_get(s, 0, i + 1, j, k);
    t.BL = mr_bf_get(s, 1, i, j, k);
    t.BR = mr_bf_get(s, 1, i, j + 1, k);
    t.CL = mr_bf_get(s, 2, i, j, k);
    t.CR = mr_bf_get(s, 2, i, j, k + 1);
    float A = 0.5f * (t.AL + t.AR);
    float B = 0.5f * (t.BL + t.BR);
    float C = 0.5f * (t.CL + t.CR);
    float ELL = mr_edge_x(s, i, j, k);
    float ELR = mr_edge_x(s, i, j, k + 1);
    float ERL = mr_edge_x(s, i, j + 1, k);
    float ERR = mr_edge_x(s, i, j + 1, k + 1);
    float FLL = mr_edge_y(s, i, j, k);
    float FLR = mr_edge_y(s, i, j, k + 1);
    float FRL = mr_edge_y(s, i + 1, j, k);
    float FRR = mr_edge_y(s, i + 1, j, k + 1);
    float GLL = mr_edge_z(s, i, j, k);
    float GLR = mr_edge_z(s, i, j + 1, k);
    float GRL = mr_edge_z(s, i + 1, j, k);
    float GRR = mr_edge_z(s, i + 1, j + 1, k);
    t.AL += ((GLR - GLL) - (FLR - FLL)) * dtdx * 0.5f;
    t.AR += ((GRR - GRL) - (FRR - FRL)) * dtdx * 0.5f;
    t.BL += (-(GRL - GLL) + (ELR - ELL)) * dtdx * 0.5f;
    t.BR += (-(GRR - GLR) + (ERR - ERL)) * dtdx * 0.5f;
    t.CL += ((FRL - FLL) - (ERL - ELL)) * dtdx * 0.5f;
    t.CR += ((FRR - FLR) - (ERR - ELR)) * dtdx * 0.5f;
    t.cell = mr_cell_load(s, i, j, k);
    mr_compute_slopes(s, i, j, k, slope, t.sx, t.sy, t.sz);
    float divv = t.sx.velocity_x + t.sy.velocity_y + t.sz.velocity_z;
    primitive_t source;
    source.density = -t.cell.velocity_x * t.sx.density - t.cell.velocity_y * t.sy.density - t.cell.velocity_z * t.sz.density - divv * t.cell.density;
    source.velocity_x = -t.cell.velocity_x * t.sx.velocity_x - t.cell.velocity_y * t.sy.velocity_x - t.cell.velocity_z * t.sz.velocity_x - t.sx.pressure / t.cell.density;
    source.velocity_y = -t.cell.velocity_x * t.sx.velocity_y - t.cell.velocity_y * t.sy.velocity_y - t.cell.velocity_z * t.sz.velocity_y - t.sy.pressure / t.cell.density;
    source.velocity_z = -t.cell.velocity_x * t.sx.velocity_z - t.cell.velocity_y * t.sy.velocity_z - t.cell.velocity_z * t.sz.velocity_z - t.sz.pressure / t.cell.density;
    source.pressure = -t.cell.velocity_x * t.sx.pressure - t.cell.velocity_y * t.sy.pressure - t.cell.velocity_z * t.sz.pressure - divv * gamma * t.cell.pressure;
    source.velocity_x += (-(B * t.sx.By + C * t.sx.Bz) + B * t.sy.Bx + C * t.sz.Bx) / t.cell.density;
    source.velocity_y += (A * t.sx.By - (A * t.sy.Bx + C * t.sy.Bz) + C * t.sz.By) / t.cell.density;
    source.velocity_z += (A * t.sx.Bz + B * t.sy.Bz - (A * t.sz.Bx + B * t.sz.By)) / t.cell.density;
    if (induction != 0) {
        source.velocity_x = 0.0f;
        source.velocity_y = 0.0f;
        source.velocity_z = 0.0f;
    }
    t.cell.density += dtdx * source.density;
    t.cell.velocity_x += dtdx * source.velocity_x;
    t.cell.velocity_y += dtdx * source.velocity_y;
    t.cell.velocity_z += dtdx * source.velocity_z;
    t.cell.pressure += dtdx * source.pressure;
    t.cell.Bx = 0.5f * (t.AL + t.AR);
    t.cell.By = 0.5f * (t.BL + t.BR);
    t.cell.Bz = 0.5f * (t.CL + t.CR);
    return t;
}

primitive_t mr_traced_face(threadgroup const float *s, mr_trace_t t, int orientation, int side, int i, int j, int k, float smallr, float smallc2) {
    primitive_t q = t.cell;
    primitive_t ds = orientation == 0 ? t.sx : orientation == 1 ? t.sy : t.sz;
    float sign = side < 0 ? -1.0f : 1.0f;
    q.density += sign * ds.density;
    q.velocity_x += sign * ds.velocity_x;
    q.velocity_y += sign * ds.velocity_y;
    q.velocity_z += sign * ds.velocity_z;
    q.pressure += sign * ds.pressure;
    if (orientation == 0) {
        q.Bx = side < 0 ? t.AL : t.AR;
        q.By = t.cell.By + sign * t.sx.By;
        q.Bz = t.cell.Bz + sign * t.sx.Bz;
    } else if (orientation == 1) {
        q.Bx = t.cell.Bx + sign * t.sy.Bx;
        q.By = side < 0 ? t.BL : t.BR;
        q.Bz = t.cell.Bz + sign * t.sy.Bz;
    } else {
        q.Bx = t.cell.Bx + sign * t.sz.Bx;
        q.By = t.cell.By + sign * t.sz.By;
        q.Bz = side < 0 ? t.CL : t.CR;
    }
    if (q.density < smallr) q.density = mr_local_get(s, 0, i, j, k);
    if (q.pressure < smallr * smallc2) q.pressure = mr_local_get(s, 4, i, j, k);
    return q;
}

int mr_shell_index(int orientation, int normal, int t1, int t2) {
    int perimeter;
    if (t2 == -1) perimeter = t1;
    else if (t2 == MR_M) perimeter = MR_M + t1;
    else if (t1 == -1) perimeter = 2 * MR_M + t2;
    else perimeter = 3 * MR_M + t2;
    return orientation * MR_SHELL_CELLS + normal + (MR_M + 1) * perimeter;
}

void mr_decode_shell(int shell, thread int &orientation, thread int &normal, thread int &t1, thread int &t2) {
    orientation = shell / MR_SHELL_CELLS;
    int local = shell - orientation * MR_SHELL_CELLS;
    normal = local % (MR_M + 1);
    int perimeter = local / (MR_M + 1);
    if (perimeter < MR_M) {
        t1 = perimeter;
        t2 = -1;
    } else if (perimeter < 2 * MR_M) {
        t1 = perimeter - MR_M;
        t2 = MR_M;
    } else if (perimeter < 3 * MR_M) {
        t1 = -1;
        t2 = perimeter - 2 * MR_M;
    } else {
        t1 = MR_M;
        t2 = perimeter - 3 * MR_M;
    }
}

mr_uct_record_t mr_get_record(threadgroup const float *s, int orientation, int normal, int t1, int t2) {
    bool interior = t1 >= 0 && t1 < MR_M && t2 >= 0 && t2 < MR_M;
    if (interior) {
        int i = orientation == 0 ? normal : t1;
        int j = orientation == 1 ? normal : orientation == 0 ? t1 : t2;
        int k = orientation == 2 ? normal : t2;
        int face = orientation * MR_FACE_CELLS + mr_face_index(orientation, i, j, k);
        return mr_interior_record_load(s, face);
    }
    return mr_shell_record_load(s, mr_shell_index(orientation, normal, t1, t2));
}

float mr_uct_edge(mr_uct_record_t f1lo, mr_uct_record_t f1hi, mr_uct_record_t f2lo, mr_uct_record_t f2hi, mr_edge_pair_t v1, mr_edge_pair_t b2, mr_edge_pair_t v2, mr_edge_pair_t b1) {
    float a1L = 0.5f * (f1lo.aL + f1hi.aL);
    float a1R = 1.0f - a1L;
    float d1L = 0.5f * (f1lo.dL + f1hi.dL);
    float d1R = 0.5f * (f1lo.dR + f1hi.dR);
    float a2L = 0.5f * (f2lo.aL + f2hi.aL);
    float a2R = 1.0f - a2L;
    float d2L = 0.5f * (f2lo.dL + f2hi.dL);
    float d2R = 0.5f * (f2lo.dR + f2hi.dR);
    return a1L * v1.left * b2.left + a1R * v1.right * b2.right - a2L * v2.left * b1.left - a2R * v2.right * b1.right + d1L * b2.left - d1R * b2.right - d2L * b1.left + d2R * b1.right;
}

mr_uct_record_t mr_shell_solve(threadgroup const float *s, int shell, float gamma, float smallr, float smallc2, float dtdx, int slope, int slope_mag, int induction, float switch_dmin, float switch_pmin) {
    int orientation;
    int normal;
    int t1;
    int t2;
    mr_decode_shell(shell, orientation, normal, t1, t2);
    int il;
    int jl;
    int kl;
    int ir;
    int jr;
    int kr;
    if (orientation == 0) {
        il = normal + 1; jl = t1 + 2; kl = t2 + 2;
        ir = normal + 2; jr = jl; kr = kl;
    } else if (orientation == 1) {
        il = t1 + 2; jl = normal + 1; kl = t2 + 2;
        ir = il; jr = normal + 2; kr = kl;
    } else {
        il = t1 + 2; jl = t2 + 2; kl = normal + 1;
        ir = il; jr = jl; kr = normal + 2;
    }
    mr_trace_t tl = mr_predict_cell(s, il, jl, kl, gamma, dtdx, slope, slope_mag, induction);
    mr_trace_t tr = mr_predict_cell(s, ir, jr, kr, gamma, dtdx, slope, slope_mag, induction);
    primitive_t left = mr_rotate_face(mr_traced_face(s, tl, orientation, 1, il, jl, kl, smallr, smallc2), orientation);
    primitive_t right = mr_rotate_face(mr_traced_face(s, tr, orientation, -1, ir, jr, kr, smallr, smallc2), orientation);
    bool use_llf = mr_switch_to_llf(min(left.density, right.density), min(left.pressure, right.pressure), switch_dmin, switch_pmin);
    if (use_llf) return mr_llf_record(left, right, gamma, smallr, smallc2);
    mr_uct_record_t record;
    mr_hlld_mhd_flux(left, right, gamma, smallr, smallc2, record, true);
    return record;
}

int mr_staged_velocity_index(int oct_idx, int orientation, int cell_idx, int slot) {
    return (((oct_idx - 1) * 3 + orientation) * TWOTONDIM + cell_idx - 1) * 2 + slot;
}

void mr_staged_velocity_set(device float *velocity, int oct_idx, int orientation, int cell_idx, mr_uct_record_t record) {
    int index = mr_staged_velocity_index(oct_idx, orientation, cell_idx, 0);
    velocity[index] = record.vt1;
    velocity[index + 1] = record.vt2;
}

bool mr_staged_velocity_get(device const float *velocity, device const int *nbor, int subgrid_idx, int first_oct, int num_octs, int orientation, int normal, int t1, int t2, int slot, thread float &value) {
    int i = orientation == 0 ? normal : t1;
    int j = orientation == 1 ? normal : orientation == 0 ? t1 : t2;
    int k = orientation == 2 ? normal : t2;
    i += 2;
    j += 2;
    k += 2;
    int i_sg = i / 2;
    int j_sg = j / 2;
    int k_sg = k / 2;
    int source_idx = mr_nbor_get(nbor, subgrid_idx, 1 + i_sg + MR_NSUBGRIDP2 * j_sg + MR_NSUBGRIDP2 * MR_NSUBGRIDP2 * k_sg);
    if (source_idx < first_oct || source_idx >= first_oct + num_octs) return false;
    int cell_idx = 1 + i % 2 + 2 * (j % 2) + 4 * (k % 2);
    value = velocity[mr_staged_velocity_index(source_idx, orientation, cell_idx, slot)];
    return true;
}

float mr_face_b_get(threadgroup const float *s, int orientation, int normal, int t1, int t2) {
    int i = orientation == 0 ? normal : t1;
    int j = orientation == 1 ? normal : orientation == 0 ? t1 : t2;
    int k = orientation == 2 ? normal : t2;
    return mr_bf_get(s, orientation, i + 2, j + 2, k + 2);
}

float mr_vanleer_slope(float qm, float q0, float qp) {
    float dl = q0 - qm;
    float dr = qp - q0;
    return dl * dr > 0.0f ? 2.0f * dl * dr / (dl + dr) : 0.0f;
}

mr_edge_pair_t mr_reconstruct_velocity(float qm, float q0, float qp, float qpp, bool has_m, bool has_p) {
    mr_edge_pair_t pair;
    pair.left = q0 + 0.5f * (has_m ? mr_vanleer_slope(qm, q0, qp) : 0.0f);
    pair.right = qp - 0.5f * (has_p ? mr_vanleer_slope(q0, qp, qpp) : 0.0f);
    return pair;
}

mr_edge_pair_t mr_reconstruct_b(float qm, float q0, float qp, float qpp, int slope_mag) {
    mr_edge_pair_t pair;
    pair.left = q0 + 0.5f * mr_face_b_slope(qm, q0, qp, slope_mag);
    pair.right = qp - 0.5f * mr_face_b_slope(q0, qp, qpp, slope_mag);
    return pair;
}

mr_edge_pair_t mr_staged_velocity_pair(device const float *velocity, device const int *nbor, int subgrid_idx, int first_oct, int num_octs, int orientation, int normal, int edge, int fixed, bool vary_t1, int slot, float raw_left, float raw_right) {
    float qm = raw_left;
    float q0 = raw_left;
    float qp = raw_right;
    float qpp = raw_right;
    int t1 = vary_t1 ? edge - 2 : fixed;
    int t2 = vary_t1 ? fixed : edge - 2;
    bool has_m = mr_staged_velocity_get(velocity, nbor, subgrid_idx, first_oct, num_octs, orientation, normal, t1, t2, slot, qm);
    t1 = vary_t1 ? edge - 1 : fixed;
    t2 = vary_t1 ? fixed : edge - 1;
    mr_staged_velocity_get(velocity, nbor, subgrid_idx, first_oct, num_octs, orientation, normal, t1, t2, slot, q0);
    t1 = vary_t1 ? edge : fixed;
    t2 = vary_t1 ? fixed : edge;
    mr_staged_velocity_get(velocity, nbor, subgrid_idx, first_oct, num_octs, orientation, normal, t1, t2, slot, qp);
    t1 = vary_t1 ? edge + 1 : fixed;
    t2 = vary_t1 ? fixed : edge + 1;
    bool has_p = mr_staged_velocity_get(velocity, nbor, subgrid_idx, first_oct, num_octs, orientation, normal, t1, t2, slot, qpp);
    return mr_reconstruct_velocity(qm, q0, qp, qpp, has_m, has_p);
}

mr_edge_pair_t mr_product_velocity_pair(device const float *product, device const int *nbor, int subgrid_idx, int head_idx, int orientation, int normal, int edge, int fixed, bool vary_t1, int slot) {
    int t1 = vary_t1 ? edge - 2 : fixed;
    int t2 = vary_t1 ? fixed : edge - 2;
    mr_uct_record_t qm_record = mr_product_record_get(product, nbor, subgrid_idx, head_idx, orientation, normal, t1, t2);
    t1 = vary_t1 ? edge - 1 : fixed;
    t2 = vary_t1 ? fixed : edge - 1;
    mr_uct_record_t q0_record = mr_product_record_get(product, nbor, subgrid_idx, head_idx, orientation, normal, t1, t2);
    t1 = vary_t1 ? edge : fixed;
    t2 = vary_t1 ? fixed : edge;
    mr_uct_record_t qp_record = mr_product_record_get(product, nbor, subgrid_idx, head_idx, orientation, normal, t1, t2);
    t1 = vary_t1 ? edge + 1 : fixed;
    t2 = vary_t1 ? fixed : edge + 1;
    mr_uct_record_t qpp_record = mr_product_record_get(product, nbor, subgrid_idx, head_idx, orientation, normal, t1, t2);
    float qm = slot == 0 ? qm_record.vt1 : qm_record.vt2;
    float q0 = slot == 0 ? q0_record.vt1 : q0_record.vt2;
    float qp = slot == 0 ? qp_record.vt1 : qp_record.vt2;
    float qpp = slot == 0 ? qpp_record.vt1 : qpp_record.vt2;
    return mr_reconstruct_velocity(qm, q0, qp, qpp, true, true);
}

mr_edge_pair_t mr_face_b_pair(threadgroup const float *s, int orientation, int normal, int edge, int fixed, bool vary_t1, int slope_mag) {
    int t1 = vary_t1 ? edge - 2 : fixed;
    int t2 = vary_t1 ? fixed : edge - 2;
    float qm = mr_face_b_get(s, orientation, normal, t1, t2);
    t1 = vary_t1 ? edge - 1 : fixed;
    t2 = vary_t1 ? fixed : edge - 1;
    float q0 = mr_face_b_get(s, orientation, normal, t1, t2);
    t1 = vary_t1 ? edge : fixed;
    t2 = vary_t1 ? fixed : edge;
    float qp = mr_face_b_get(s, orientation, normal, t1, t2);
    t1 = vary_t1 ? edge + 1 : fixed;
    t2 = vary_t1 ? fixed : edge + 1;
    float qpp = mr_face_b_get(s, orientation, normal, t1, t2);
    return mr_reconstruct_b(qm, q0, qp, qpp, slope_mag);
}

void subgrid_conserved_2_primitive_mhd(device const oct_t *grid, device const float *uold, device const float *bold, device const float *f, device const int *nbor, constant float *constant_gravity, int subgrid_idx, float gamma, float smallr, float smallc2, float dt, threadgroup float *s, threadgroup bool *refined, int tid, int threads_per_group) {
    for (int work = tid; work < MR_LOCAL_CELLS; work += threads_per_group) {
        int oct_lattice = work / TWOTONDIM;
        int i_sg;
        int j_sg;
        int k_sg;
        index_1Dto3D(oct_lattice, MR_NSUBGRIDP2, MR_NSUBGRIDP2, i_sg, j_sg, k_sg);
        int ind_nbor = oct_lattice + 1;
        int source_idx = mr_nbor_get(nbor, subgrid_idx, ind_nbor);
        int cell_idx = work % TWOTONDIM + 1;
        int ib;
        int jb;
        int kb;
        index_1Dto3D(cell_idx - 1, 2, 2, ib, jb, kb);
        int i = ib + 2 * i_sg;
        int j = jb + 2 * j_sg;
        int k = kb + 2 * k_sg;
        float b0[6];
        for (int component = 0; component < 6; ++component) b0[component] = b_get(bold, source_idx, component + 1, cell_idx);
        if (i >= 1) {
            int face_source = ib == 0 ? mr_nbor_get(nbor, subgrid_idx, ind_nbor - 1) : source_idx;
            int face_cell = cell_idx + 1 - 2 * ib;
            mr_bf_set(s, 0, i, j, k, 0.5f * (b0[0] + b_get(bold, face_source, 4, face_cell)));
        }
        if (j >= 1) {
            int face_source = jb == 0 ? mr_nbor_get(nbor, subgrid_idx, ind_nbor - MR_NSUBGRIDP2) : source_idx;
            int face_cell = cell_idx + 2 * (1 - 2 * jb);
            mr_bf_set(s, 1, i, j, k, 0.5f * (b0[1] + b_get(bold, face_source, 5, face_cell)));
        }
        if (k >= 1) {
            int face_source = kb == 0 ? mr_nbor_get(nbor, subgrid_idx, ind_nbor - MR_NSUBGRIDP2 * MR_NSUBGRIDP2) : source_idx;
            int face_cell = cell_idx + 4 * (1 - 2 * kb);
            mr_bf_set(s, 2, i, j, k, 0.5f * (b0[2] + b_get(bold, face_source, 6, face_cell)));
        }
        conserved_t u;
        u.density = u_get(uold, source_idx, 1, cell_idx);
        u.momentum_x = u_get(uold, source_idx, 2, cell_idx);
        u.momentum_y = u_get(uold, source_idx, 3, cell_idx);
        u.momentum_z = u_get(uold, source_idx, 4, cell_idx);
        u.energy = u_get(uold, source_idx, 5, cell_idx);
        u.Bx = 0.5f * (b0[0] + b0[3]);
        u.By = 0.5f * (b0[1] + b0[4]);
        u.Bz = 0.5f * (b0[2] + b0[5]);
        primitive_t q = conserved_2_primitive(u, gamma, smallr, smallc2);
#ifdef GRAV
        int grav_base = (source_idx - 1) * 3 * TWOTONDIM + cell_idx - 1;
        q.velocity_x += 0.5f * dt * f[grav_base];
        q.velocity_y += 0.5f * dt * f[grav_base + TWOTONDIM];
        q.velocity_z += 0.5f * dt * f[grav_base + 2 * TWOTONDIM];
#else
        q.velocity_x += 0.5f * dt * constant_gravity[0];
        q.velocity_y += 0.5f * dt * constant_gravity[1];
        q.velocity_z += 0.5f * dt * constant_gravity[2];
#endif
        mr_local_set(s, 0, i, j, k, q.density);
        mr_local_set(s, 1, i, j, k, q.velocity_x);
        mr_local_set(s, 2, i, j, k, q.velocity_y);
        mr_local_set(s, 3, i, j, k, q.velocity_z);
        mr_local_set(s, 4, i, j, k, q.pressure);
        if (i >= 1 && i <= MR_TRACE && j >= 1 && j <= MR_TRACE && k >= 1 && k <= MR_TRACE) mr_refined_set(refined, i, j, k, grid[source_idx - 1].refined[cell_idx - 1] != 0);
    }
}

void mr_load_b_stencil(device const float *bold, device const int *nbor, int subgrid_idx, threadgroup float *s, int tid, int threads_per_group) {
    for (int work = tid; work < MR_LOCAL_CELLS; work += threads_per_group) {
        int oct_lattice = work / TWOTONDIM;
        int i_sg;
        int j_sg;
        int k_sg;
        index_1Dto3D(oct_lattice, MR_NSUBGRIDP2, MR_NSUBGRIDP2, i_sg, j_sg, k_sg);
        int ind_nbor = oct_lattice + 1;
        int source_idx = mr_nbor_get(nbor, subgrid_idx, ind_nbor);
        int cell_idx = work % TWOTONDIM + 1;
        int ib;
        int jb;
        int kb;
        index_1Dto3D(cell_idx - 1, 2, 2, ib, jb, kb);
        int i = ib + 2 * i_sg;
        int j = jb + 2 * j_sg;
        int k = kb + 2 * k_sg;
        float b0[6];
        for (int component = 0; component < 6; ++component) b0[component] = b_get(bold, source_idx, component + 1, cell_idx);
        if (i >= 1) {
            int face_source = ib == 0 ? mr_nbor_get(nbor, subgrid_idx, ind_nbor - 1) : source_idx;
            int face_cell = cell_idx + 1 - 2 * ib;
            mr_bf_set(s, 0, i, j, k, 0.5f * (b0[0] + b_get(bold, face_source, 4, face_cell)));
        }
        if (j >= 1) {
            int face_source = jb == 0 ? mr_nbor_get(nbor, subgrid_idx, ind_nbor - MR_NSUBGRIDP2) : source_idx;
            int face_cell = cell_idx + 2 * (1 - 2 * jb);
            mr_bf_set(s, 1, i, j, k, 0.5f * (b0[1] + b_get(bold, face_source, 5, face_cell)));
        }
        if (k >= 1) {
            int face_source = kb == 0 ? mr_nbor_get(nbor, subgrid_idx, ind_nbor - MR_NSUBGRIDP2 * MR_NSUBGRIDP2) : source_idx;
            int face_cell = cell_idx + 4 * (1 - 2 * kb);
            mr_bf_set(s, 2, i, j, k, 0.5f * (b0[2] + b_get(bold, face_source, 6, face_cell)));
        }
    }
}

struct mr_cf_cell_t {
    bool coarse;
    int oct_idx;
    int cell_idx;
};

mr_cf_cell_t mr_get_cf_cell(device const oct_t *grid, device const int *father, device const int *nbor, int subgrid_idx, int ngridmax, int i, int j, int k) {
    mr_cf_cell_t cell;
    cell.coarse = false;
    cell.oct_idx = 0;
    cell.cell_idx = 1;
    int stencil_idx = 1 + i + MR_NSUBGRIDP2 * j + MR_NSUBGRIDP2 * MR_NSUBGRIDP2 * k;
    int nbor_idx = mr_nbor_get(nbor, subgrid_idx, stencil_idx);
    if (nbor_idx > ngridmax) {
        int father_idx = father[nbor_idx - 1];
        int ci = grid[nbor_idx - 1].ckey[0] - 2 * grid[father_idx - 1].ckey[0];
        int cj = grid[nbor_idx - 1].ckey[1] - 2 * grid[father_idx - 1].ckey[1];
        int ck = grid[nbor_idx - 1].ckey[2] - 2 * grid[father_idx - 1].ckey[2];
        cell.coarse = true;
        cell.oct_idx = father_idx;
        cell.cell_idx = 1 + ci + 2 * cj + 4 * ck;
    }
    return cell;
}

void mr_u_atomic_add(device float *u, int oct_idx, int ivar, int cell_idx, float addend) {
    atomic_add_float((device atomic_uint *)&u[u_flat(oct_idx, ivar, cell_idx)], addend);
}

void mr_cf_b_add(device float *b, mr_cf_cell_t cell, int ivar, float addend) {
    if (cell.coarse) atomic_add_float((device atomic_uint *)&b[b_flat(cell.oct_idx, ivar, cell.cell_idx)], addend);
}

float mr_cf_weight(mr_cf_cell_t c1, mr_cf_cell_t c2, mr_cf_cell_t c3) {
    return c1.coarse && c2.coarse && c3.coarse ? 1.0f : 0.5f;
}

void mr_coarse_hydro_update(threadgroup const float *s, device const oct_t *grid, device const int *father, device const int *nbor, device float *unew, int subgrid_idx, int ngridmax, int tid, float dtdx) {
    int interface_size = NSUBGRID * NSUBGRID;
    if (tid >= 6 * interface_size) return;
    int boundary = tid / interface_size;
    int work = tid - boundary * interface_size;
    int first = work % NSUBGRID;
    int second = work / NSUBGRID;
    int orientation = boundary / 2;
    int side = boundary % 2;
    int normal = side == 0 ? 0 : MR_M;
    int i_sg = orientation == 0 ? (side == 0 ? 0 : NSUBGRID + 1) : first + 1;
    int j_sg = orientation == 1 ? (side == 0 ? 0 : NSUBGRID + 1) : orientation == 0 ? first + 1 : second + 1;
    int k_sg = orientation == 2 ? (side == 0 ? 0 : NSUBGRID + 1) : second + 1;
    float flux[5] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    float scale = (side == 0 ? -0.125f : 0.125f) * dtdx;
    for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
            int t1 = 2 * first + a;
            int t2 = 2 * second + b;
            int i = orientation == 0 ? normal : t1;
            int j = orientation == 1 ? normal : orientation == 0 ? t1 : t2;
            int k = orientation == 2 ? normal : t2;
            int face = orientation * MR_FACE_CELLS + mr_face_index(orientation, i, j, k);
            for (int field = 0; field < 5; ++field) flux[field] += mr_flux_get(s, field, face) * scale;
        }
    }
    mr_cf_cell_t cell = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg, k_sg);
    if (cell.coarse) for (int field = 0; field < 5; ++field) mr_u_atomic_add(unew, cell.oct_idx, field + 1, cell.cell_idx, flux[field]);
}

void mr_coarse_ct_update(threadgroup const float *s, device const oct_t *grid, device const int *father, device const int *nbor, device float *bnew, int subgrid_idx, int ngridmax, float dtdx) {
    for (int k_sg = 1; k_sg <= NSUBGRID; ++k_sg) {
        int k = 2 * k_sg - 2;
        for (int j_sg = 1; j_sg <= NSUBGRID; ++j_sg) {
            int j = 2 * j_sg - 2;
            for (int i_sg = 1; i_sg <= NSUBGRID; ++i_sg) {
                int i = 2 * i_sg - 2;
                mr_cf_cell_t c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg - 1, k_sg);
                mr_cf_cell_t c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg - 1, j_sg - 1, k_sg);
                mr_cf_cell_t c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg - 1, j_sg, k_sg);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 0, i, j, k) + mr_emf_get(s, 0, i, j, k + 1)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 5, dflux); mr_cf_b_add(bnew, c1, 1, dflux);
                    mr_cf_b_add(bnew, c2, 4, dflux); mr_cf_b_add(bnew, c2, 5, -dflux);
                    mr_cf_b_add(bnew, c3, 2, -dflux); mr_cf_b_add(bnew, c3, 4, -dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg - 1, j_sg, k_sg);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg - 1, j_sg + 1, k_sg);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg + 1, k_sg);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 0, i, j + 2, k) + mr_emf_get(s, 0, i, j + 2, k + 1)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 4, dflux); mr_cf_b_add(bnew, c1, 5, -dflux);
                    mr_cf_b_add(bnew, c2, 2, -dflux); mr_cf_b_add(bnew, c2, 4, -dflux);
                    mr_cf_b_add(bnew, c3, 1, -dflux); mr_cf_b_add(bnew, c3, 2, dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg + 1, k_sg);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg + 1, j_sg + 1, k_sg);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg + 1, j_sg, k_sg);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 0, i + 2, j + 2, k) + mr_emf_get(s, 0, i + 2, j + 2, k + 1)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 2, -dflux); mr_cf_b_add(bnew, c1, 4, -dflux);
                    mr_cf_b_add(bnew, c2, 1, -dflux); mr_cf_b_add(bnew, c2, 2, dflux);
                    mr_cf_b_add(bnew, c3, 5, dflux); mr_cf_b_add(bnew, c3, 1, dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg + 1, j_sg, k_sg);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg + 1, j_sg - 1, k_sg);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg - 1, k_sg);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 0, i + 2, j, k) + mr_emf_get(s, 0, i + 2, j, k + 1)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 1, -dflux); mr_cf_b_add(bnew, c1, 2, dflux);
                    mr_cf_b_add(bnew, c2, 5, dflux); mr_cf_b_add(bnew, c2, 1, dflux);
                    mr_cf_b_add(bnew, c3, 4, dflux); mr_cf_b_add(bnew, c3, 5, -dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg, k_sg - 1);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg - 1, k_sg - 1);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg - 1, k_sg);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 2, i, j, k) + mr_emf_get(s, 2, i + 1, j, k)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 6, dflux); mr_cf_b_add(bnew, c1, 2, dflux);
                    mr_cf_b_add(bnew, c2, 5, dflux); mr_cf_b_add(bnew, c2, 6, -dflux);
                    mr_cf_b_add(bnew, c3, 3, -dflux); mr_cf_b_add(bnew, c3, 5, -dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg - 1, k_sg);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg - 1, k_sg + 1);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg, k_sg + 1);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 2, i, j, k + 2) + mr_emf_get(s, 2, i + 1, j, k + 2)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 5, dflux); mr_cf_b_add(bnew, c1, 6, -dflux);
                    mr_cf_b_add(bnew, c2, 3, -dflux); mr_cf_b_add(bnew, c2, 5, -dflux);
                    mr_cf_b_add(bnew, c3, 2, -dflux); mr_cf_b_add(bnew, c3, 3, dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg, k_sg + 1);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg + 1, k_sg + 1);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg + 1, k_sg);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 2, i, j + 2, k + 2) + mr_emf_get(s, 2, i + 1, j + 2, k + 2)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 3, -dflux); mr_cf_b_add(bnew, c1, 5, -dflux);
                    mr_cf_b_add(bnew, c2, 2, -dflux); mr_cf_b_add(bnew, c2, 3, dflux);
                    mr_cf_b_add(bnew, c3, 6, dflux); mr_cf_b_add(bnew, c3, 2, dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg + 1, k_sg);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg + 1, k_sg - 1);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg, k_sg - 1);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 2, i, j + 2, k) + mr_emf_get(s, 2, i + 1, j + 2, k)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 2, -dflux); mr_cf_b_add(bnew, c1, 3, dflux);
                    mr_cf_b_add(bnew, c2, 6, dflux); mr_cf_b_add(bnew, c2, 2, dflux);
                    mr_cf_b_add(bnew, c3, 5, dflux); mr_cf_b_add(bnew, c3, 6, -dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg, k_sg - 1);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg - 1, j_sg, k_sg - 1);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg - 1, j_sg, k_sg);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 1, i, j, k) + mr_emf_get(s, 1, i, j + 1, k)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 6, -dflux); mr_cf_b_add(bnew, c1, 1, -dflux);
                    mr_cf_b_add(bnew, c2, 4, -dflux); mr_cf_b_add(bnew, c2, 6, dflux);
                    mr_cf_b_add(bnew, c3, 3, dflux); mr_cf_b_add(bnew, c3, 4, dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg - 1, j_sg, k_sg);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg - 1, j_sg, k_sg + 1);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg, k_sg + 1);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 1, i, j, k + 2) + mr_emf_get(s, 1, i, j + 1, k + 2)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 4, -dflux); mr_cf_b_add(bnew, c1, 6, dflux);
                    mr_cf_b_add(bnew, c2, 3, dflux); mr_cf_b_add(bnew, c2, 4, dflux);
                    mr_cf_b_add(bnew, c3, 1, dflux); mr_cf_b_add(bnew, c3, 3, -dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg, k_sg + 1);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg + 1, j_sg, k_sg + 1);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg + 1, j_sg, k_sg);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 1, i + 2, j, k + 2) + mr_emf_get(s, 1, i + 2, j + 1, k + 2)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 3, dflux); mr_cf_b_add(bnew, c1, 4, dflux);
                    mr_cf_b_add(bnew, c2, 1, dflux); mr_cf_b_add(bnew, c2, 3, -dflux);
                    mr_cf_b_add(bnew, c3, 6, -dflux); mr_cf_b_add(bnew, c3, 1, -dflux);
                }
                c1 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg + 1, j_sg, k_sg);
                c2 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg + 1, j_sg, k_sg - 1);
                c3 = mr_get_cf_cell(grid, father, nbor, subgrid_idx, ngridmax, i_sg, j_sg, k_sg - 1);
                if (c1.coarse || c3.coarse) {
                    float dflux = (mr_emf_get(s, 1, i + 2, j, k) + mr_emf_get(s, 1, i + 2, j + 1, k)) * 0.25f * mr_cf_weight(c1, c2, c3) * dtdx;
                    mr_cf_b_add(bnew, c1, 1, dflux); mr_cf_b_add(bnew, c1, 3, -dflux);
                    mr_cf_b_add(bnew, c2, 6, -dflux); mr_cf_b_add(bnew, c2, 1, -dflux);
                    mr_cf_b_add(bnew, c3, 4, -dflux); mr_cf_b_add(bnew, c3, 6, dflux);
                }
            }
        }
    }
}

kernel void uct_velocity_kernel(
    device const oct_t *grid [[buffer(0)]],
    device const float *uold [[buffer(1)]],
    device const float *bold [[buffer(2)]],
    device const int *nbor [[buffer(3)]],
    constant int &head_idx [[buffer(4)]],
    constant int &num_subgrids [[buffer(5)]],
    constant float &gamma [[buffer(6)]],
    constant float &smallr [[buffer(7)]],
    constant float &smallc2 [[buffer(8)]],
    constant float &dt [[buffer(9)]],
    constant float &dx [[buffer(10)]],
    constant int &slope [[buffer(11)]],
    constant int &slope_mag [[buffer(12)]],
    constant float &switch_llf_dmin [[buffer(13)]],
    constant float &switch_llf_pmin [[buffer(14)]],
    constant int &induction [[buffer(15)]],
    constant float *constant_gravity [[buffer(16)]],
    device const float *f [[buffer(17)]],
    device float *velocity [[buffer(18)]],
    uint block_idx [[threadgroup_position_in_grid]],
    uint tid [[thread_position_in_threadgroup]],
    uint threads_per_group [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_subgrids) return;
    threadgroup float state[MR_STATE_FLOATS];
    threadgroup bool refined[MR_REFINED_CELLS];
    int subgrid_idx = head_idx + int(block_idx);
    float dtdx = dt / dx;
    subgrid_conserved_2_primitive_mhd(grid, uold, bold, f, nbor, constant_gravity, subgrid_idx, gamma, smallr, smallc2, dt, state, refined, int(tid), int(threads_per_group));
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (int work = int(tid); work < 3 * MR_INTERIOR_CELLS; work += int(threads_per_group)) {
        int orientation = work / MR_INTERIOR_CELLS;
        int i;
        int j;
        int k;
        index_1Dto3D(work % MR_INTERIOR_CELLS, MR_M, MR_M, i, j, k);
        int ir = i + 2;
        int jr = j + 2;
        int kr = k + 2;
        int il = ir;
        int jl = jr;
        int kl = kr;
        if (orientation == 0) --il;
        else if (orientation == 1) --jl;
        else --kl;
        mr_trace_t tl = mr_predict_cell(state, il, jl, kl, gamma, dtdx, slope, slope_mag, induction);
        mr_trace_t tr = mr_predict_cell(state, ir, jr, kr, gamma, dtdx, slope, slope_mag, induction);
        primitive_t left = mr_rotate_face(mr_traced_face(state, tl, orientation, 1, il, jl, kl, smallr, smallc2), orientation);
        primitive_t right = mr_rotate_face(mr_traced_face(state, tr, orientation, -1, ir, jr, kr, smallr, smallc2), orientation);
        bool use_llf = mr_switch_to_llf(min(left.density, right.density), min(left.pressure, right.pressure), switch_llf_dmin, switch_llf_pmin);
        mr_uct_record_t record;
        if (use_llf) record = mr_llf_record(left, right, gamma, smallr, smallc2);
        else mr_hlld_mhd_flux(left, right, gamma, smallr, smallc2, record, true);
        int source_idx = mr_nbor_get(nbor, subgrid_idx, 1 + 1 + i / 2 + MR_NSUBGRIDP2 * (1 + j / 2) + MR_NSUBGRIDP2 * MR_NSUBGRIDP2 * (1 + k / 2));
        int cell_idx = 1 + i % 2 + 2 * (j % 2) + 4 * (k % 2);
        mr_staged_velocity_set(velocity, source_idx, orientation, cell_idx, record);
    }
}

kernel void mhd_uct_face_product_kernel(
    device const oct_t *grid [[buffer(0)]],
    device const float *uold [[buffer(1)]],
    device const float *bold [[buffer(2)]],
    device const int *nbor [[buffer(3)]],
    constant int &head_idx [[buffer(4)]],
    constant int &num_subgrids [[buffer(5)]],
    constant float &gamma [[buffer(6)]],
    constant float &smallr [[buffer(7)]],
    constant float &smallc2 [[buffer(8)]],
    constant float &dt [[buffer(9)]],
    constant float &dx [[buffer(10)]],
    constant int &slope [[buffer(11)]],
    constant int &slope_mag [[buffer(12)]],
    constant float &switch_llf_dmin [[buffer(13)]],
    constant float &switch_llf_pmin [[buffer(14)]],
    constant int &induction [[buffer(15)]],
    constant float *constant_gravity [[buffer(16)]],
    device const float *f [[buffer(17)]],
    device float *product [[buffer(18)]],
    uint block_idx [[threadgroup_position_in_grid]],
    uint tid [[thread_position_in_threadgroup]],
    uint threads_per_group [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_subgrids) return;
    threadgroup float state[MR_STATE_FLOATS];
    threadgroup bool refined[MR_REFINED_CELLS];
    int subgrid_idx = head_idx + int(block_idx);
    float dtdx = dt / dx;
    subgrid_conserved_2_primitive_mhd(grid, uold, bold, f, nbor, constant_gravity, subgrid_idx, gamma, smallr, smallc2, dt, state, refined, int(tid), int(threads_per_group));
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (int work = int(tid); work < MR_UCT_PRODUCT_FACES; work += int(threads_per_group)) {
        int orientation = work / MR_INTERIOR_CELLS;
        int i;
        int j;
        int k;
        index_1Dto3D(work % MR_INTERIOR_CELLS, MR_M, MR_M, i, j, k);
        int ir = i + 2;
        int jr = j + 2;
        int kr = k + 2;
        int il = ir;
        int jl = jr;
        int kl = kr;
        if (orientation == 0) --il;
        else if (orientation == 1) --jl;
        else --kl;
        mr_trace_t tl = mr_predict_cell(state, il, jl, kl, gamma, dtdx, slope, slope_mag, induction);
        mr_trace_t tr = mr_predict_cell(state, ir, jr, kr, gamma, dtdx, slope, slope_mag, induction);
        primitive_t left = mr_rotate_face(mr_traced_face(state, tl, orientation, 1, il, jl, kl, smallr, smallc2), orientation);
        primitive_t right = mr_rotate_face(mr_traced_face(state, tr, orientation, -1, ir, jr, kr, smallr, smallc2), orientation);
        bool use_llf = mr_switch_to_llf(min(left.density, right.density), min(left.pressure, right.pressure), switch_llf_dmin, switch_llf_pmin);
        mr_uct_record_t record;
        conserved_t flux;
        if (use_llf) {
            float llf_lmax;
            flux = mr_hll_mhd_flux_lmax(left, right, gamma, smallr, smallc2, llf_lmax);
            record = mr_llf_record_lmax(left, right, 0.5f * (left.Bx + right.Bx), llf_lmax);
        } else {
            flux = mr_hlld_mhd_flux(left, right, gamma, smallr, smallc2, record, true);
        }
        flux = mr_unrotate_flux(flux, orientation);
        if (induction != 0) {
            flux.density = 0.0f;
            flux.momentum_x = 0.0f;
            flux.momentum_y = 0.0f;
            flux.momentum_z = 0.0f;
            flux.energy = 0.0f;
        }
        mr_product_set(product, subgrid_idx, head_idx, 0, work, flux.density);
        mr_product_set(product, subgrid_idx, head_idx, 1, work, flux.momentum_x);
        mr_product_set(product, subgrid_idx, head_idx, 2, work, flux.momentum_y);
        mr_product_set(product, subgrid_idx, head_idx, 3, work, flux.momentum_z);
        mr_product_set(product, subgrid_idx, head_idx, 4, work, flux.energy);
        for (int field = 0; field < MR_UCT_PRODUCT_RECORD_FIELDS; ++field) mr_product_set(product, subgrid_idx, head_idx, MR_UCT_PRODUCT_FLUX_FIELDS + field, work, mr_record_value(record, field));
    }
}

void trace_3d_mhd(threadgroup float *s, int tid, float gamma, float smallr, float smallc2, float dtdx, int slope, int slope_mag, int induction) {
    if (tid < MR_REFINED_CELLS) {
        int i;
        int j;
        int k;
        index_1Dto3D(tid, MR_TRACE, MR_TRACE, i, j, k);
        ++i;
        ++j;
        ++k;
        mr_trace_t trace = mr_predict_cell(s, i, j, k, gamma, dtdx, slope, slope_mag, induction);
        if (i > 1 && j > 1 && j < MR_TRACE && k > 1 && k < MR_TRACE) {
            int index = mr_face_index(0, i - 2, j - 2, k - 2);
            mr_face_store(s, 1, index, mr_traced_face(s, trace, 0, -1, i, j, k, smallr, smallc2));
        }
        if (i < MR_TRACE && j > 1 && j < MR_TRACE && k > 1 && k < MR_TRACE) {
            int index = mr_face_index(0, i - 1, j - 2, k - 2);
            mr_face_store(s, 0, index, mr_traced_face(s, trace, 0, 1, i, j, k, smallr, smallc2));
        }
        if (i > 1 && i < MR_TRACE && j > 1 && k > 1 && k < MR_TRACE) {
            int index = mr_face_index(1, i - 2, j - 2, k - 2);
            mr_face_store(s, 3, index, mr_traced_face(s, trace, 1, -1, i, j, k, smallr, smallc2));
        }
        if (i > 1 && i < MR_TRACE && j < MR_TRACE && k > 1 && k < MR_TRACE) {
            int index = mr_face_index(1, i - 2, j - 1, k - 2);
            mr_face_store(s, 2, index, mr_traced_face(s, trace, 1, 1, i, j, k, smallr, smallc2));
        }
        if (i > 1 && i < MR_TRACE && j > 1 && j < MR_TRACE && k > 1) {
            int index = mr_face_index(2, i - 2, j - 2, k - 2);
            mr_face_store(s, 5, index, mr_traced_face(s, trace, 2, -1, i, j, k, smallr, smallc2));
        }
        if (i > 1 && i < MR_TRACE && j > 1 && j < MR_TRACE && k < MR_TRACE) {
            int index = mr_face_index(2, i - 2, j - 2, k - 1);
            mr_face_store(s, 4, index, mr_traced_face(s, trace, 2, 1, i, j, k, smallr, smallc2));
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

float mr_face_transverse_slope(threadgroup const float *s, int component, int i, int j, int k, int axis, int slope_mag) {
    int im = i;
    int jm = j;
    int km = k;
    int ip = i;
    int jp = j;
    int kp = k;
    if (axis == 0) {
        --im;
        ++ip;
    } else if (axis == 1) {
        --jm;
        ++jp;
    } else {
        --km;
        ++kp;
    }
    return 0.5f * mr_face_b_slope(mr_bf_get(s, component, im, jm, km), mr_bf_get(s, component, i, j, k), mr_bf_get(s, component, ip, jp, kp), slope_mag);
}

primitive_t mr_rotate_corner(primitive_t q, int orientation) {
    primitive_t r = q;
    r.velocity_z = 0.0f;
    if (orientation == 1) {
        r.velocity_x = q.velocity_z;
        r.velocity_y = q.velocity_x;
        r.Bx = q.Bz;
        r.By = q.Bx;
        r.Bz = q.By;
    } else if (orientation == 2) {
        r.velocity_x = q.velocity_y;
        r.velocity_y = q.velocity_z;
        r.Bx = q.By;
        r.By = q.Bz;
        r.Bz = q.Bx;
    }
    return r;
}

primitive_t mr_hlld_corner_state(threadgroup const float *s, mr_trace_t t, int orientation, float sign1, float sign2, int i, int j, int k, int slope_mag, float smallr, float smallc2) {
    primitive_t q = t.cell;
    primitive_t d1 = orientation == 0 ? t.sx : orientation == 1 ? t.sx : t.sy;
    primitive_t d2 = orientation == 0 ? t.sy : t.sz;
    q.density += sign1 * d1.density + sign2 * d2.density;
    q.velocity_x += sign1 * d1.velocity_x + sign2 * d2.velocity_x;
    q.velocity_y += sign1 * d1.velocity_y + sign2 * d2.velocity_y;
    q.velocity_z += sign1 * d1.velocity_z + sign2 * d2.velocity_z;
    q.pressure += sign1 * d1.pressure + sign2 * d2.pressure;
    if (orientation == 0) {
        int ix = i + (sign1 > 0.0f ? 1 : 0);
        int jy = j + (sign2 > 0.0f ? 1 : 0);
        q.Bx = (sign1 > 0.0f ? t.AR : t.AL) + sign2 * mr_face_transverse_slope(s, 0, ix, j, k, 1, slope_mag);
        q.By = (sign2 > 0.0f ? t.BR : t.BL) + sign1 * mr_face_transverse_slope(s, 1, i, jy, k, 0, slope_mag);
        q.Bz = t.cell.Bz + sign1 * t.sx.Bz + sign2 * t.sy.Bz;
    } else if (orientation == 1) {
        int ix = i + (sign1 > 0.0f ? 1 : 0);
        int kz = k + (sign2 > 0.0f ? 1 : 0);
        q.Bx = (sign1 > 0.0f ? t.AR : t.AL) + sign2 * mr_face_transverse_slope(s, 0, ix, j, k, 2, slope_mag);
        q.By = t.cell.By + sign1 * t.sx.By + sign2 * t.sz.By;
        q.Bz = (sign2 > 0.0f ? t.CR : t.CL) + sign1 * mr_face_transverse_slope(s, 2, i, j, kz, 0, slope_mag);
    } else {
        int jy = j + (sign1 > 0.0f ? 1 : 0);
        int kz = k + (sign2 > 0.0f ? 1 : 0);
        q.Bx = t.cell.Bx + sign1 * t.sy.Bx + sign2 * t.sz.Bx;
        q.By = (sign1 > 0.0f ? t.BR : t.BL) + sign2 * mr_face_transverse_slope(s, 1, i, jy, k, 2, slope_mag);
        q.Bz = (sign2 > 0.0f ? t.CR : t.CL) + sign1 * mr_face_transverse_slope(s, 2, i, j, kz, 1, slope_mag);
    }
    if (q.density < smallr) q.density = mr_local_get(s, 0, i, j, k);
    if (q.pressure < smallr * smallc2) q.pressure = mr_local_get(s, 4, i, j, k);
    return mr_rotate_corner(q, orientation);
}

void mr_hlld_corner_stage(threadgroup float *s, int orientation, mr_trace_t t, int i, int j, int k, int slope_mag, float smallr, float smallc2) {
    if (orientation == 0 && k > 1 && k < MR_TRACE) {
        if (i < MR_TRACE && j < MR_TRACE) mr_hlld_corner_set(s, mr_hlld_edge_index(0, i - 1, j - 1, k - 2), 0, mr_hlld_corner_state(s, t, 0, 1.0f, 1.0f, i, j, k, slope_mag, smallr, smallc2));
        if (i < MR_TRACE && j > 1) mr_hlld_corner_set(s, mr_hlld_edge_index(0, i - 1, j - 2, k - 2), 2, mr_hlld_corner_state(s, t, 0, 1.0f, -1.0f, i, j, k, slope_mag, smallr, smallc2));
        if (i > 1 && j < MR_TRACE) mr_hlld_corner_set(s, mr_hlld_edge_index(0, i - 2, j - 1, k - 2), 1, mr_hlld_corner_state(s, t, 0, -1.0f, 1.0f, i, j, k, slope_mag, smallr, smallc2));
        if (i > 1 && j > 1) mr_hlld_corner_set(s, mr_hlld_edge_index(0, i - 2, j - 2, k - 2), 3, mr_hlld_corner_state(s, t, 0, -1.0f, -1.0f, i, j, k, slope_mag, smallr, smallc2));
    } else if (orientation == 1 && j > 1 && j < MR_TRACE) {
        if (i < MR_TRACE && k < MR_TRACE) mr_hlld_corner_set(s, mr_hlld_edge_index(1, i - 1, j - 2, k - 1), 0, mr_hlld_corner_state(s, t, 1, 1.0f, 1.0f, i, j, k, slope_mag, smallr, smallc2));
        if (i < MR_TRACE && k > 1) mr_hlld_corner_set(s, mr_hlld_edge_index(1, i - 1, j - 2, k - 2), 1, mr_hlld_corner_state(s, t, 1, 1.0f, -1.0f, i, j, k, slope_mag, smallr, smallc2));
        if (i > 1 && k < MR_TRACE) mr_hlld_corner_set(s, mr_hlld_edge_index(1, i - 2, j - 2, k - 1), 2, mr_hlld_corner_state(s, t, 1, -1.0f, 1.0f, i, j, k, slope_mag, smallr, smallc2));
        if (i > 1 && k > 1) mr_hlld_corner_set(s, mr_hlld_edge_index(1, i - 2, j - 2, k - 2), 3, mr_hlld_corner_state(s, t, 1, -1.0f, -1.0f, i, j, k, slope_mag, smallr, smallc2));
    } else if (orientation == 2 && i > 1 && i < MR_TRACE) {
        if (j < MR_TRACE && k < MR_TRACE) mr_hlld_corner_set(s, mr_hlld_edge_index(2, i - 2, j - 1, k - 1), 0, mr_hlld_corner_state(s, t, 2, 1.0f, 1.0f, i, j, k, slope_mag, smallr, smallc2));
        if (j < MR_TRACE && k > 1) mr_hlld_corner_set(s, mr_hlld_edge_index(2, i - 2, j - 1, k - 2), 2, mr_hlld_corner_state(s, t, 2, 1.0f, -1.0f, i, j, k, slope_mag, smallr, smallc2));
        if (j > 1 && k < MR_TRACE) mr_hlld_corner_set(s, mr_hlld_edge_index(2, i - 2, j - 2, k - 1), 1, mr_hlld_corner_state(s, t, 2, -1.0f, 1.0f, i, j, k, slope_mag, smallr, smallc2));
        if (j > 1 && k > 1) mr_hlld_corner_set(s, mr_hlld_edge_index(2, i - 2, j - 2, k - 2), 3, mr_hlld_corner_state(s, t, 2, -1.0f, -1.0f, i, j, k, slope_mag, smallr, smallc2));
    }
}

float mr_hlld_edge_solve(threadgroup const float *s, int edge, float gamma, float smallr, float smallc2, float switch_llf_dmin, float switch_llf_pmin, int riemann2d) {
    primitive_t qLL = mr_hlld_corner_get(s, edge, 0);
    primitive_t qRL = mr_hlld_corner_get(s, edge, 1);
    primitive_t qLR = mr_hlld_corner_get(s, edge, 2);
    primitive_t qRR = mr_hlld_corner_get(s, edge, 3);
    float b1LL = 0.5f * (qLL.Bx + qRL.Bx);
    float b1LR = 0.5f * (qLR.Bx + qRR.Bx);
    qLL.Bx = b1LL;
    qRL.Bx = b1LL;
    qLR.Bx = b1LR;
    qRR.Bx = b1LR;
    float b2LL = 0.5f * (qLL.By + qLR.By);
    float b2RL = 0.5f * (qRL.By + qRR.By);
    qLL.By = b2LL;
    qLR.By = b2LL;
    qRL.By = b2RL;
    qRR.By = b2RL;
    if (riemann2d == MR_SOLVER2D_HLLD) return mr_hlld_2d_emf(qLL, qLR, qRL, qRR, gamma, smallr, sqrt(smallc2), switch_llf_dmin, switch_llf_pmin);
    return mr_llf_2d_emf(qLL, qLR, qRL, qRR, gamma, smallr, smallc2);
}

void hlld_emf_driver(threadgroup float *s, threadgroup const bool *refined, int ilevel, int levelmax, float gamma, float smallr, float smallc2, float dtdx, float dx, float etamag, int slope, int slope_mag, int induction, float switch_llf_dmin, float switch_llf_pmin, int riemann2d, int tid, int threads_per_group) {
    int i = 0;
    int j = 0;
    int k = 0;
    mr_trace_t trace;
    if (tid < MR_REFINED_CELLS) {
        index_1Dto3D(tid, MR_TRACE, MR_TRACE, i, j, k);
        ++i;
        ++j;
        ++k;
        trace = mr_predict_cell(s, i, j, k, gamma, dtdx, slope, slope_mag, induction);
    }
    for (int orientation = 0; orientation < 3; ++orientation) {
        if (tid < MR_REFINED_CELLS) mr_hlld_corner_stage(s, orientation, trace, i, j, k, slope_mag, smallr, smallc2);
        threadgroup_barrier(mem_flags::mem_threadgroup);
        if (tid < MR_EDGE_CELLS) {
            int I;
            int J;
            int K;
            float emf = mr_hlld_edge_solve(s, tid, gamma, smallr, smallc2, switch_llf_dmin, switch_llf_pmin, riemann2d);
            if (orientation == 0) {
                index_1Dto3D(tid, MR_M + 1, MR_M + 1, I, J, K);
                if (etamag > 0.0f) {
                    int X = I + 2;
                    int Y = J + 2;
                    int Z = K + 2;
                    emf -= etamag / dx * ((mr_bf_get(s, 1, X, Y, Z) - mr_bf_get(s, 1, X - 1, Y, Z)) - (mr_bf_get(s, 0, X, Y, Z) - mr_bf_get(s, 0, X, Y - 1, Z)));
                }
                if (ilevel < levelmax) {
                    int X = I + 2;
                    int Y = J + 2;
                    int Z = K + 2;
                    if (mr_refined_get(refined, X - 1, Y - 1, Z) || mr_refined_get(refined, X, Y - 1, Z) || mr_refined_get(refined, X - 1, Y, Z) || mr_refined_get(refined, X, Y, Z)) emf = 0.0f;
                }
            } else if (orientation == 1) {
                index_1Dto3D(tid, MR_M + 1, MR_M, I, J, K);
                if (etamag > 0.0f) {
                    int X = I + 2;
                    int Y = J + 2;
                    int Z = K + 2;
                    emf -= etamag / dx * ((mr_bf_get(s, 0, X, Y, Z) - mr_bf_get(s, 0, X, Y, Z - 1)) - (mr_bf_get(s, 2, X, Y, Z) - mr_bf_get(s, 2, X - 1, Y, Z)));
                }
                if (ilevel < levelmax) {
                    int X = I + 2;
                    int Y = J + 2;
                    int Z = K + 2;
                    if (mr_refined_get(refined, X - 1, Y, Z - 1) || mr_refined_get(refined, X, Y, Z - 1) || mr_refined_get(refined, X - 1, Y, Z) || mr_refined_get(refined, X, Y, Z)) emf = 0.0f;
                }
            } else {
                index_1Dto3D(tid, MR_M, MR_M + 1, I, J, K);
                if (etamag > 0.0f) {
                    int X = I + 2;
                    int Y = J + 2;
                    int Z = K + 2;
                    emf -= etamag / dx * ((mr_bf_get(s, 2, X, Y, Z) - mr_bf_get(s, 2, X, Y - 1, Z)) - (mr_bf_get(s, 1, X, Y, Z) - mr_bf_get(s, 1, X, Y, Z - 1)));
                }
                if (ilevel < levelmax) {
                    int X = I + 2;
                    int Y = J + 2;
                    int Z = K + 2;
                    if (mr_refined_get(refined, X, Y - 1, Z - 1) || mr_refined_get(refined, X, Y, Z - 1) || mr_refined_get(refined, X, Y - 1, Z) || mr_refined_get(refined, X, Y, Z)) emf = 0.0f;
                }
            }
            mr_hlld_emf_set(s, orientation, tid, emf);
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    for (int edge = tid; edge < MR_ALL_EDGES; edge += threads_per_group) s[MR_EMF_BASE + edge] = s[MR_HLLD_EMF_BASE + edge];
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

void riemann_driver_mhd(threadgroup float *s, int tid, float gamma, float smallr, float smallc2, float dtdx, int slope, int slope_mag, int induction, int riemann, float switch_llf_dmin, float switch_llf_pmin) {
    bool use_uct = riemann == MR_SOLVER_UCT_HLLD;
    primitive_t left_face;
    primitive_t right_face;
    int orientation = tid / MR_FACE_CELLS;
    int face_local = tid - orientation * MR_FACE_CELLS;
    if (tid < MR_ALL_FACES) {
        left_face = mr_face_load(s, 2 * orientation, face_local);
        right_face = mr_face_load(s, 2 * orientation + 1, face_local);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (tid < MR_ALL_FACES) {
        left_face = mr_rotate_face(left_face, orientation);
        right_face = mr_rotate_face(right_face, orientation);
        bool use_llf = riemann == MR_SOLVER_LLF || mr_switch_to_llf(min(left_face.density, right_face.density), min(left_face.pressure, right_face.pressure), switch_llf_dmin, switch_llf_pmin);
        mr_uct_record_t record;
        conserved_t flux;
        if (use_llf) {
            float llf_lmax;
            flux = mr_hll_mhd_flux_lmax(left_face, right_face, gamma, smallr, smallc2, llf_lmax);
            if (use_uct) record = mr_llf_record_lmax(left_face, right_face, 0.5f * (left_face.Bx + right_face.Bx), llf_lmax);
        } else {
            flux = mr_hlld_mhd_flux(left_face, right_face, gamma, smallr, smallc2, record, use_uct);
        }
        flux = mr_unrotate_flux(flux, orientation);
        if (induction != 0) {
            flux.density = 0.0f;
            flux.momentum_x = 0.0f;
            flux.momentum_y = 0.0f;
            flux.momentum_z = 0.0f;
            flux.energy = 0.0f;
        }
        mr_flux_set(s, 0, tid, flux.density);
        mr_flux_set(s, 1, tid, flux.momentum_x);
        mr_flux_set(s, 2, tid, flux.momentum_y);
        mr_flux_set(s, 3, tid, flux.momentum_z);
        mr_flux_set(s, 4, tid, flux.energy);
        if (use_uct) mr_interior_record_store(s, tid, record);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (use_uct) {
        mr_uct_record_t shell_record;
        if (tid < MR_ALL_SHELLS) shell_record = mr_shell_solve(s, tid, gamma, smallr, smallc2, dtdx, slope, slope_mag, induction, switch_llf_dmin, switch_llf_pmin);
        threadgroup_barrier(mem_flags::mem_threadgroup);
        if (tid < MR_ALL_SHELLS) mr_shell_record_store(s, tid, shell_record);
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
}

void uct_emf_driver(threadgroup float *s, threadgroup const bool *refined, device const float *velocity, device const int *nbor, int subgrid_idx, int head_idx, int num_subgrids, int ilevel, int levelmax, float etamag, float dx, int slope_mag, int tid, int threads_per_group) {
    int first_oct = (head_idx - 1) * MR_NSUBGRID_CELLS + 1;
    int num_octs = num_subgrids * MR_NSUBGRID_CELLS;
    for (int edge = tid; edge < MR_ALL_EDGES; edge += threads_per_group) {
        int I;
        int J;
        int K;
        float emf;
        int edge_orientation;
        if (edge < MR_EDGE_CELLS) {
            index_1Dto3D(edge, MR_M + 1, MR_M + 1, I, J, K);
            mr_uct_record_t f1lo = mr_get_record(s, 0, I, J - 1, K);
            mr_uct_record_t f1hi = mr_get_record(s, 0, I, J, K);
            mr_uct_record_t f2lo = mr_get_record(s, 1, J, I - 1, K);
            mr_uct_record_t f2hi = mr_get_record(s, 1, J, I, K);
            mr_edge_pair_t v1 = mr_staged_velocity_pair(velocity, nbor, subgrid_idx, first_oct, num_octs, 1, J, I, K, true, 0, f2lo.vt1, f2hi.vt1);
            mr_edge_pair_t b2 = mr_face_b_pair(s, 1, J, I, K, true, slope_mag);
            mr_edge_pair_t v2 = mr_staged_velocity_pair(velocity, nbor, subgrid_idx, first_oct, num_octs, 0, I, J, K, true, 0, f1lo.vt1, f1hi.vt1);
            mr_edge_pair_t b1 = mr_face_b_pair(s, 0, I, J, K, true, slope_mag);
            emf = mr_uct_edge(f1lo, f1hi, f2lo, f2hi, v1, b2, v2, b1);
            edge_orientation = 0;
            if (etamag > 0.0f) {
                int X = I + 2;
                int Y = J + 2;
                int Z = K + 2;
                emf -= etamag / dx * ((mr_bf_get(s, 1, X, Y, Z) - mr_bf_get(s, 1, X - 1, Y, Z)) - (mr_bf_get(s, 0, X, Y, Z) - mr_bf_get(s, 0, X, Y - 1, Z)));
            }
            if (ilevel < levelmax) {
                int X = I + 2;
                int Y = J + 2;
                int Z = K + 2;
                if (mr_refined_get(refined, X - 1, Y - 1, Z) || mr_refined_get(refined, X, Y - 1, Z) || mr_refined_get(refined, X - 1, Y, Z) || mr_refined_get(refined, X, Y, Z)) emf = 0.0f;
            }
        } else if (edge < 2 * MR_EDGE_CELLS) {
            index_1Dto3D(edge - MR_EDGE_CELLS, MR_M + 1, MR_M, I, J, K);
            mr_uct_record_t f1lo = mr_get_record(s, 2, K, I - 1, J);
            mr_uct_record_t f1hi = mr_get_record(s, 2, K, I, J);
            mr_uct_record_t f2lo = mr_get_record(s, 0, I, J, K - 1);
            mr_uct_record_t f2hi = mr_get_record(s, 0, I, J, K);
            mr_edge_pair_t v1 = mr_staged_velocity_pair(velocity, nbor, subgrid_idx, first_oct, num_octs, 0, I, K, J, false, 1, f2lo.vt2, f2hi.vt2);
            mr_edge_pair_t b2 = mr_face_b_pair(s, 0, I, K, J, false, slope_mag);
            mr_edge_pair_t v2 = mr_staged_velocity_pair(velocity, nbor, subgrid_idx, first_oct, num_octs, 2, K, I, J, true, 0, f1lo.vt1, f1hi.vt1);
            mr_edge_pair_t b1 = mr_face_b_pair(s, 2, K, I, J, true, slope_mag);
            emf = mr_uct_edge(f1lo, f1hi, f2lo, f2hi, v1, b2, v2, b1);
            edge_orientation = 1;
            if (etamag > 0.0f) {
                int X = I + 2;
                int Y = J + 2;
                int Z = K + 2;
                emf -= etamag / dx * ((mr_bf_get(s, 0, X, Y, Z) - mr_bf_get(s, 0, X, Y, Z - 1)) - (mr_bf_get(s, 2, X, Y, Z) - mr_bf_get(s, 2, X - 1, Y, Z)));
            }
            if (ilevel < levelmax) {
                int X = I + 2;
                int Y = J + 2;
                int Z = K + 2;
                if (mr_refined_get(refined, X - 1, Y, Z - 1) || mr_refined_get(refined, X, Y, Z - 1) || mr_refined_get(refined, X - 1, Y, Z) || mr_refined_get(refined, X, Y, Z)) emf = 0.0f;
            }
        } else {
            index_1Dto3D(edge - 2 * MR_EDGE_CELLS, MR_M, MR_M + 1, I, J, K);
            mr_uct_record_t f1lo = mr_get_record(s, 1, J, I, K - 1);
            mr_uct_record_t f1hi = mr_get_record(s, 1, J, I, K);
            mr_uct_record_t f2lo = mr_get_record(s, 2, K, I, J - 1);
            mr_uct_record_t f2hi = mr_get_record(s, 2, K, I, J);
            mr_edge_pair_t v1 = mr_staged_velocity_pair(velocity, nbor, subgrid_idx, first_oct, num_octs, 2, K, J, I, false, 1, f2lo.vt2, f2hi.vt2);
            mr_edge_pair_t b2 = mr_face_b_pair(s, 2, K, J, I, false, slope_mag);
            mr_edge_pair_t v2 = mr_staged_velocity_pair(velocity, nbor, subgrid_idx, first_oct, num_octs, 1, J, K, I, false, 1, f1lo.vt2, f1hi.vt2);
            mr_edge_pair_t b1 = mr_face_b_pair(s, 1, J, K, I, false, slope_mag);
            emf = mr_uct_edge(f1lo, f1hi, f2lo, f2hi, v1, b2, v2, b1);
            edge_orientation = 2;
            if (etamag > 0.0f) {
                int X = I + 2;
                int Y = J + 2;
                int Z = K + 2;
                emf -= etamag / dx * ((mr_bf_get(s, 2, X, Y, Z) - mr_bf_get(s, 2, X, Y - 1, Z)) - (mr_bf_get(s, 1, X, Y, Z) - mr_bf_get(s, 1, X, Y, Z - 1)));
            }
            if (ilevel < levelmax) {
                int X = I + 2;
                int Y = J + 2;
                int Z = K + 2;
                if (mr_refined_get(refined, X, Y - 1, Z - 1) || mr_refined_get(refined, X, Y, Z - 1) || mr_refined_get(refined, X, Y - 1, Z) || mr_refined_get(refined, X, Y, Z)) emf = 0.0f;
            }
        }
        mr_emf_set(s, edge_orientation, I, J, K, emf);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

void zero_fine_fluxes_mhd(threadgroup float *s, threadgroup const bool *refined, int ilevel, int levelmax, int tid) {
    int orientation = tid / MR_FACE_CELLS;
    int face_local = tid - orientation * MR_FACE_CELLS;
    if (tid < MR_ALL_FACES && ilevel < levelmax) {
        int i;
        int j;
        int k;
        if (orientation == 0) {
            index_1Dto3D(face_local, MR_M + 1, MR_M, i, j, k);
            if (mr_refined_get(refined, i + 1, j + 2, k + 2) || mr_refined_get(refined, i + 2, j + 2, k + 2)) for (int field = 0; field < 5; ++field) mr_flux_set(s, field, tid, 0.0f);
        } else if (orientation == 1) {
            index_1Dto3D(face_local, MR_M, MR_M + 1, i, j, k);
            if (mr_refined_get(refined, i + 2, j + 1, k + 2) || mr_refined_get(refined, i + 2, j + 2, k + 2)) for (int field = 0; field < 5; ++field) mr_flux_set(s, field, tid, 0.0f);
        } else {
            index_1Dto3D(face_local, MR_M, MR_M, i, j, k);
            if (mr_refined_get(refined, i + 2, j + 2, k + 1) || mr_refined_get(refined, i + 2, j + 2, k + 2)) for (int field = 0; field < 5; ++field) mr_flux_set(s, field, tid, 0.0f);
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

void conservative_update_mhd(device float *unew, device float *bnew, device const int *nbor, threadgroup const float *s, int subgrid_idx, int tid, float dtdx, bool update_hydro, bool update_ct) {
    if (tid < MR_INTERIOR_CELLS) {
        int i_sg;
        int j_sg;
        int k_sg;
        index_1Dto3D(tid / TWOTONDIM, NSUBGRID, NSUBGRID, i_sg, j_sg, k_sg);
        ++i_sg;
        ++j_sg;
        ++k_sg;
        int ind_nbor = 1 + i_sg + MR_NSUBGRIDP2 * j_sg + MR_NSUBGRIDP2 * MR_NSUBGRIDP2 * k_sg;
        int oct_idx = mr_nbor_get(nbor, subgrid_idx, ind_nbor);
        int cell_idx = tid % TWOTONDIM + 1;
        int i;
        int j;
        int k;
        index_1Dto3D(cell_idx - 1, 2, 2, i, j, k);
        i += 2 * (i_sg - 1);
        j += 2 * (j_sg - 1);
        k += 2 * (k_sg - 1);
        if (update_hydro) {
            int fx0 = mr_face_index(0, i, j, k);
            int fx1 = mr_face_index(0, i + 1, j, k);
            int fy0 = MR_FACE_CELLS + mr_face_index(1, i, j, k);
            int fy1 = MR_FACE_CELLS + mr_face_index(1, i, j + 1, k);
            int fz0 = 2 * MR_FACE_CELLS + mr_face_index(2, i, j, k);
            int fz1 = 2 * MR_FACE_CELLS + mr_face_index(2, i, j, k + 1);
            for (int field = 0; field < 5; ++field) {
                float update = (mr_flux_get(s, field, fx0) - mr_flux_get(s, field, fx1) + mr_flux_get(s, field, fy0) - mr_flux_get(s, field, fy1) + mr_flux_get(s, field, fz0) - mr_flux_get(s, field, fz1)) * dtdx;
                u_set(unew, oct_idx, field + 1, cell_idx, u_get(unew, oct_idx, field + 1, cell_idx) + update);
            }
        }
        if (update_ct) {
            float dflux = ((mr_emf_get(s, 1, i, j, k) - mr_emf_get(s, 1, i, j, k + 1)) - (mr_emf_get(s, 0, i, j, k) - mr_emf_get(s, 0, i, j + 1, k))) * dtdx;
            b_set(bnew, oct_idx, 1, cell_idx, b_get(bnew, oct_idx, 1, cell_idx) + dflux);
            dflux = ((mr_emf_get(s, 1, i + 1, j, k) - mr_emf_get(s, 1, i + 1, j, k + 1)) - (mr_emf_get(s, 0, i + 1, j, k) - mr_emf_get(s, 0, i + 1, j + 1, k))) * dtdx;
            b_set(bnew, oct_idx, 4, cell_idx, b_get(bnew, oct_idx, 4, cell_idx) + dflux);
            dflux = ((mr_emf_get(s, 0, i, j, k) - mr_emf_get(s, 0, i + 1, j, k)) - (mr_emf_get(s, 2, i, j, k) - mr_emf_get(s, 2, i, j, k + 1))) * dtdx;
            b_set(bnew, oct_idx, 2, cell_idx, b_get(bnew, oct_idx, 2, cell_idx) + dflux);
            dflux = ((mr_emf_get(s, 0, i, j + 1, k) - mr_emf_get(s, 0, i + 1, j + 1, k)) - (mr_emf_get(s, 2, i, j + 1, k) - mr_emf_get(s, 2, i, j + 1, k + 1))) * dtdx;
            b_set(bnew, oct_idx, 5, cell_idx, b_get(bnew, oct_idx, 5, cell_idx) + dflux);
            dflux = ((mr_emf_get(s, 2, i, j, k) - mr_emf_get(s, 2, i, j + 1, k)) - (mr_emf_get(s, 1, i, j, k) - mr_emf_get(s, 1, i + 1, j, k))) * dtdx;
            b_set(bnew, oct_idx, 3, cell_idx, b_get(bnew, oct_idx, 3, cell_idx) + dflux);
            dflux = ((mr_emf_get(s, 2, i, j, k + 1) - mr_emf_get(s, 2, i, j + 1, k + 1)) - (mr_emf_get(s, 1, i, j, k + 1) - mr_emf_get(s, 1, i + 1, j, k + 1))) * dtdx;
            b_set(bnew, oct_idx, 6, cell_idx, b_get(bnew, oct_idx, 6, cell_idx) + dflux);
        }
    }
}

void coarse_cell_update_mhd(threadgroup const float *s, device const oct_t *grid, device const int *father, device const int *nbor, device float *unew, device float *bnew, int subgrid_idx, int ngridmax, int tid, float dtdx) {
    mr_coarse_hydro_update(s, grid, father, nbor, unew, subgrid_idx, ngridmax, tid, dtdx);
    if (tid == 0) mr_coarse_ct_update(s, grid, father, nbor, bnew, subgrid_idx, ngridmax, dtdx);
}

kernel void hydro_integrator_kernel(
    device const oct_t *grid [[buffer(0)]],
    device const float *uold [[buffer(1)]],
    device float *unew [[buffer(2)]],
    device const float *bold [[buffer(3)]],
    device float *bnew [[buffer(4)]],
    device const int *father [[buffer(5)]],
    device const int *nbor [[buffer(6)]],
    constant int &head_idx [[buffer(7)]],
    constant int &num_subgrids [[buffer(8)]],
    constant int &ngridmax [[buffer(9)]],
    constant int &ilevel [[buffer(10)]],
    constant int &levelmin [[buffer(11)]],
    constant int &levelmax [[buffer(12)]],
    constant float &gamma [[buffer(13)]],
    constant float &smallr [[buffer(14)]],
    constant float &smallc2 [[buffer(15)]],
    constant float &dt [[buffer(16)]],
    constant float &dx [[buffer(17)]],
    constant int &slope [[buffer(18)]],
    constant int &slope_mag [[buffer(19)]],
    constant int &riemann [[buffer(20)]],
    constant int &riemann2d [[buffer(21)]],
    constant float &switch_llf_dmin [[buffer(22)]],
    constant float &switch_llf_pmin [[buffer(23)]],
    constant int &induction [[buffer(24)]],
    constant float &etamag [[buffer(25)]],
    constant float *constant_gravity [[buffer(26)]],
    device const float *f [[buffer(27)]],
    device const float *velocity [[buffer(28)]],
    uint block_idx [[threadgroup_position_in_grid]],
    uint tid [[thread_position_in_threadgroup]],
    uint threads_per_group [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_subgrids) return;
    threadgroup float smem[MR_SMEM_FLOATS];
    threadgroup bool refined[MR_REFINED_CELLS];
    int subgrid_idx = head_idx + int(block_idx);
    float dtdx = dt / dx;
    subgrid_conserved_2_primitive_mhd(grid, uold, bold, f, nbor, constant_gravity, subgrid_idx, gamma, smallr, smallc2, dt, smem, refined, int(tid), int(threads_per_group));
    threadgroup_barrier(mem_flags::mem_threadgroup);
    trace_3d_mhd(smem, int(tid), gamma, smallr, smallc2, dtdx, slope, slope_mag, induction);
    riemann_driver_mhd(smem, int(tid), gamma, smallr, smallc2, dtdx, slope, slope_mag, induction, riemann, switch_llf_dmin, switch_llf_pmin);
    if (riemann == MR_SOLVER_UCT_HLLD) uct_emf_driver(smem, refined, velocity, nbor, subgrid_idx, head_idx, num_subgrids, ilevel, levelmax, etamag, dx, slope_mag, int(tid), int(threads_per_group));
    zero_fine_fluxes_mhd(smem, refined, ilevel, levelmax, int(tid));
    if (riemann == MR_SOLVER_UCT_HLLD) {
        conservative_update_mhd(unew, bnew, nbor, smem, subgrid_idx, int(tid), dtdx, true, true);
        if (ilevel > levelmin) coarse_cell_update_mhd(smem, grid, father, nbor, unew, bnew, subgrid_idx, ngridmax, int(tid), dtdx);
    } else if (riemann == MR_SOLVER_HLLD || riemann == MR_SOLVER_LLF) {
        conservative_update_mhd(unew, bnew, nbor, smem, subgrid_idx, int(tid), dtdx, true, false);
        if (ilevel > levelmin) mr_coarse_hydro_update(smem, grid, father, nbor, unew, subgrid_idx, ngridmax, int(tid), dtdx);
        threadgroup_barrier(mem_flags::mem_threadgroup);
        hlld_emf_driver(smem, refined, ilevel, levelmax, gamma, smallr, smallc2, dtdx, dx, etamag, slope, slope_mag, induction, switch_llf_dmin, switch_llf_pmin, riemann2d, int(tid), int(threads_per_group));
        conservative_update_mhd(unew, bnew, nbor, smem, subgrid_idx, int(tid), dtdx, false, true);
        if (ilevel > levelmin && tid == 0) mr_coarse_ct_update(smem, grid, father, nbor, bnew, subgrid_idx, ngridmax, dtdx);
    }
}

kernel void hydro_integrator_uct_reuse_kernel(
    device float *unew [[buffer(0)]],
    device const float *bold [[buffer(1)]],
    device float *bnew [[buffer(2)]],
    device const int *nbor [[buffer(3)]],
    constant int &head_idx [[buffer(4)]],
    constant int &num_subgrids [[buffer(5)]],
    constant float &dt [[buffer(6)]],
    constant float &dx [[buffer(7)]],
    constant int &slope_mag [[buffer(8)]],
    constant float &etamag [[buffer(9)]],
    device const float *product [[buffer(10)]],
    uint block_idx [[threadgroup_position_in_grid]],
    uint tid [[thread_position_in_threadgroup]],
    uint threads_per_group [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_subgrids) return;
    threadgroup float smem[MR_STATE_FLOATS];
    int subgrid_idx = head_idx + int(block_idx);
    float dtdx = dt / dx;
    mr_load_b_stencil(bold, nbor, subgrid_idx, smem, int(tid), int(threads_per_group));
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (int edge = int(tid); edge < MR_ALL_EDGES; edge += int(threads_per_group)) {
        int I;
        int J;
        int K;
        float emf;
        int edge_orientation;
        if (edge < MR_EDGE_CELLS) {
            index_1Dto3D(edge, MR_M + 1, MR_M + 1, I, J, K);
            mr_uct_record_t f1lo = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 0, I, J - 1, K);
            mr_uct_record_t f1hi = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 0, I, J, K);
            mr_uct_record_t f2lo = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 1, J, I - 1, K);
            mr_uct_record_t f2hi = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 1, J, I, K);
            mr_edge_pair_t v1 = mr_product_velocity_pair(product, nbor, subgrid_idx, head_idx, 1, J, I, K, true, 0);
            mr_edge_pair_t b2 = mr_face_b_pair(smem, 1, J, I, K, true, slope_mag);
            mr_edge_pair_t v2 = mr_product_velocity_pair(product, nbor, subgrid_idx, head_idx, 0, I, J, K, true, 0);
            mr_edge_pair_t b1 = mr_face_b_pair(smem, 0, I, J, K, true, slope_mag);
            emf = mr_uct_edge(f1lo, f1hi, f2lo, f2hi, v1, b2, v2, b1);
            edge_orientation = 0;
            if (etamag > 0.0f) {
                int X = I + 2;
                int Y = J + 2;
                int Z = K + 2;
                emf -= etamag / dx * ((mr_bf_get(smem, 1, X, Y, Z) - mr_bf_get(smem, 1, X - 1, Y, Z)) - (mr_bf_get(smem, 0, X, Y, Z) - mr_bf_get(smem, 0, X, Y - 1, Z)));
            }
        } else if (edge < 2 * MR_EDGE_CELLS) {
            index_1Dto3D(edge - MR_EDGE_CELLS, MR_M + 1, MR_M, I, J, K);
            mr_uct_record_t f1lo = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 2, K, I - 1, J);
            mr_uct_record_t f1hi = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 2, K, I, J);
            mr_uct_record_t f2lo = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 0, I, J, K - 1);
            mr_uct_record_t f2hi = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 0, I, J, K);
            mr_edge_pair_t v1 = mr_product_velocity_pair(product, nbor, subgrid_idx, head_idx, 0, I, K, J, false, 1);
            mr_edge_pair_t b2 = mr_face_b_pair(smem, 0, I, K, J, false, slope_mag);
            mr_edge_pair_t v2 = mr_product_velocity_pair(product, nbor, subgrid_idx, head_idx, 2, K, I, J, true, 0);
            mr_edge_pair_t b1 = mr_face_b_pair(smem, 2, K, I, J, true, slope_mag);
            emf = mr_uct_edge(f1lo, f1hi, f2lo, f2hi, v1, b2, v2, b1);
            edge_orientation = 1;
            if (etamag > 0.0f) {
                int X = I + 2;
                int Y = J + 2;
                int Z = K + 2;
                emf -= etamag / dx * ((mr_bf_get(smem, 0, X, Y, Z) - mr_bf_get(smem, 0, X, Y, Z - 1)) - (mr_bf_get(smem, 2, X, Y, Z) - mr_bf_get(smem, 2, X - 1, Y, Z)));
            }
        } else {
            index_1Dto3D(edge - 2 * MR_EDGE_CELLS, MR_M, MR_M + 1, I, J, K);
            mr_uct_record_t f1lo = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 1, J, I, K - 1);
            mr_uct_record_t f1hi = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 1, J, I, K);
            mr_uct_record_t f2lo = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 2, K, I, J - 1);
            mr_uct_record_t f2hi = mr_product_record_get(product, nbor, subgrid_idx, head_idx, 2, K, I, J);
            mr_edge_pair_t v1 = mr_product_velocity_pair(product, nbor, subgrid_idx, head_idx, 2, K, J, I, false, 1);
            mr_edge_pair_t b2 = mr_face_b_pair(smem, 2, K, J, I, false, slope_mag);
            mr_edge_pair_t v2 = mr_product_velocity_pair(product, nbor, subgrid_idx, head_idx, 1, J, K, I, false, 1);
            mr_edge_pair_t b1 = mr_face_b_pair(smem, 1, J, K, I, false, slope_mag);
            emf = mr_uct_edge(f1lo, f1hi, f2lo, f2hi, v1, b2, v2, b1);
            edge_orientation = 2;
            if (etamag > 0.0f) {
                int X = I + 2;
                int Y = J + 2;
                int Z = K + 2;
                emf -= etamag / dx * ((mr_bf_get(smem, 2, X, Y, Z) - mr_bf_get(smem, 2, X, Y - 1, Z)) - (mr_bf_get(smem, 1, X, Y, Z) - mr_bf_get(smem, 1, X, Y, Z - 1)));
            }
        }
        mr_emf_set(smem, edge_orientation, I, J, K, emf);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (tid < MR_INTERIOR_CELLS) {
        int i_sg;
        int j_sg;
        int k_sg;
        index_1Dto3D(int(tid) / TWOTONDIM, NSUBGRID, NSUBGRID, i_sg, j_sg, k_sg);
        ++i_sg;
        ++j_sg;
        ++k_sg;
        int ind_nbor = 1 + i_sg + MR_NSUBGRIDP2 * j_sg + MR_NSUBGRIDP2 * MR_NSUBGRIDP2 * k_sg;
        int oct_idx = mr_nbor_get(nbor, subgrid_idx, ind_nbor);
        int cell_idx = int(tid) % TWOTONDIM + 1;
        int i;
        int j;
        int k;
        index_1Dto3D(cell_idx - 1, 2, 2, i, j, k);
        i += 2 * (i_sg - 1);
        j += 2 * (j_sg - 1);
        k += 2 * (k_sg - 1);
        for (int field = 0; field < 5; ++field) {
            float update = (mr_product_flux_get(product, nbor, subgrid_idx, head_idx, 0, i, j, k, field) - mr_product_flux_get(product, nbor, subgrid_idx, head_idx, 0, i + 1, j, k, field) + mr_product_flux_get(product, nbor, subgrid_idx, head_idx, 1, j, i, k, field) - mr_product_flux_get(product, nbor, subgrid_idx, head_idx, 1, j + 1, i, k, field) + mr_product_flux_get(product, nbor, subgrid_idx, head_idx, 2, k, i, j, field) - mr_product_flux_get(product, nbor, subgrid_idx, head_idx, 2, k + 1, i, j, field)) * dtdx;
            u_set(unew, oct_idx, field + 1, cell_idx, u_get(unew, oct_idx, field + 1, cell_idx) + update);
        }
        float dflux = ((mr_emf_get(smem, 1, i, j, k) - mr_emf_get(smem, 1, i, j, k + 1)) - (mr_emf_get(smem, 0, i, j, k) - mr_emf_get(smem, 0, i, j + 1, k))) * dtdx;
        b_set(bnew, oct_idx, 1, cell_idx, b_get(bnew, oct_idx, 1, cell_idx) + dflux);
        dflux = ((mr_emf_get(smem, 1, i + 1, j, k) - mr_emf_get(smem, 1, i + 1, j, k + 1)) - (mr_emf_get(smem, 0, i + 1, j, k) - mr_emf_get(smem, 0, i + 1, j + 1, k))) * dtdx;
        b_set(bnew, oct_idx, 4, cell_idx, b_get(bnew, oct_idx, 4, cell_idx) + dflux);
        dflux = ((mr_emf_get(smem, 0, i, j, k) - mr_emf_get(smem, 0, i + 1, j, k)) - (mr_emf_get(smem, 2, i, j, k) - mr_emf_get(smem, 2, i, j, k + 1))) * dtdx;
        b_set(bnew, oct_idx, 2, cell_idx, b_get(bnew, oct_idx, 2, cell_idx) + dflux);
        dflux = ((mr_emf_get(smem, 0, i, j + 1, k) - mr_emf_get(smem, 0, i + 1, j + 1, k)) - (mr_emf_get(smem, 2, i, j + 1, k) - mr_emf_get(smem, 2, i, j + 1, k + 1))) * dtdx;
        b_set(bnew, oct_idx, 5, cell_idx, b_get(bnew, oct_idx, 5, cell_idx) + dflux);
        dflux = ((mr_emf_get(smem, 2, i, j, k) - mr_emf_get(smem, 2, i, j + 1, k)) - (mr_emf_get(smem, 1, i, j, k) - mr_emf_get(smem, 1, i + 1, j, k))) * dtdx;
        b_set(bnew, oct_idx, 3, cell_idx, b_get(bnew, oct_idx, 3, cell_idx) + dflux);
        dflux = ((mr_emf_get(smem, 2, i, j, k + 1) - mr_emf_get(smem, 2, i, j + 1, k + 1)) - (mr_emf_get(smem, 1, i, j, k + 1) - mr_emf_get(smem, 1, i + 1, j, k + 1))) * dtdx;
        b_set(bnew, oct_idx, 6, cell_idx, b_get(bnew, oct_idx, 6, cell_idx) + dflux);
    }
}

/*
 * metal/kernels/hydro.metal
 *
 * Metal Shading Language port of gpu/gpu_hydro.cuf.
 * Scope: HYDRO=1, NDIM=3, NPRE=4 (float32 throughout); MHD extends this file below.
 *
 * Structure intentionally mirrors gpu_hydro.cuf so the CUDA and Metal kernels
 * are recognisable side-by-side.  Each device function corresponds 1:1 to its
 * attributes(device) counterpart and each kernel to an attributes(global) one.
 *
 * Global buffer layout (hydro kernels):
 *   set_unew / set_uold:  [0]=uold  [1]=unew  [2]=head_idx  [3]=num_octs
 *   cmpdt:  [0]=grid [1]=uold [2]=data_buf [3]=head_idx [4]=num_octs
 *           [5]=dx [6]=gamma [7]=smallr [8]=smallc2 [9]=courant_factor
 *           [10]=constant_gravity(3 floats)
 *   hydro_integrator:  [0]=grid [1]=uold [2]=unew [3]=nbor
 *           [4]=head_idx [5]=num_subgrids [6]=ngridmax
 *           [7]=ilevel [8]=levelmin [9]=levelmax
 *           [10]=gamma [11]=smallr [12]=smallc2 [13]=dt [14]=dx
 *           [15]=slope [16]=riemann [17]=constant_gravity(3 floats)
 */

#include <metal_stdlib>
#include <metal_atomic>
using namespace metal;

#include "../metal_types.h"
#include "metal_utils.h"

/* ---------------------------------------------------------------------------
 * Compile-time constants — injected via xcrun metal -DNDIM=3 etc.
 * --------------------------------------------------------------------------*/
#ifndef NVAR
#define NVAR 5
#endif

constant int TWOTONDIM   = 8;    /* 2^NDIM                      */
constant int NSUBGRID    = 1;    /* subgrid size for PoC        */
constant int NSUBGRIDP2  = 3;    /* NSUBGRID + 2                */
constant int NTHREADS_Y  = 16;   /* threadgroup y-dim for set_* */

/* Riemann solver IDs (mirrors hydro_parameters.f90) */
constant int SOLVER_LLF  = 1;
constant int SOLVER_HLL  = 2;
constant int SOLVER_HLLC = 3;

/* ===========================================================================
 * Types — mirror gpu_hydro.cuf
 * ========================================================================= */

struct conserved_t {
    float density, momentum_x, momentum_y, momentum_z, energy;
#ifdef MHD
    float Bx, By, Bz;
#endif
};
struct primitive_t {
    float density, velocity_x, velocity_y, velocity_z, pressure;
#ifdef MHD
    float Bx, By, Bz;
#endif
};

/* 6×6×6 stencil subgrid (nsubgrid=1 → indices 0..5 per axis).
 * Matches subgrid_6x6x6cell_primitive in gpu_utils.cuf.
 * C row-major layout differs from Fortran column-major, but all accesses
 * are threadgroup-local so the ordering does not matter for correctness. */
struct local_subgrid_t {
    float density   [6][6][6];
    float velocity_x[6][6][6];
    float velocity_y[6][6][6];
    float velocity_z[6][6][6];
    float pressure  [6][6][6];
    /* Fortran: refined(1:2*nsubgrid+2, 1:2*nsubgrid+2, 1:2*nsubgrid+2) = 4×4×4
     * for nsubgrid=1.  Stored 0-based here. */
    bool  refined   [4][4][4];
};

/* Interface state arrays (nsubgrid=1):
 *   x: Fortran (0:2, 0:1, 0:1) → C [3][2][2]   (subgrid_3x2x2cell_primitive)
 *   y: Fortran (0:1, 0:2, 0:1) → C [2][3][2]   (subgrid_2x3x2cell_primitive)
 *   z: Fortran (0:1, 0:1, 0:2) → C [2][2][3]   (subgrid_2x2x3cell_primitive)
 * After riemann_driver these arrays are overwritten with the computed fluxes
 * (conserved) while keeping the primitive member names — same convention as CUDA. */
struct interfaces_x_t {
    float density   [3][2][2];
    float velocity_x[3][2][2];
    float velocity_y[3][2][2];
    float velocity_z[3][2][2];
    float pressure  [3][2][2];
};
struct interfaces_y_t {
    float density   [2][3][2];
    float velocity_x[2][3][2];
    float velocity_y[2][3][2];
    float velocity_z[2][3][2];
    float pressure  [2][3][2];
};
struct interfaces_z_t {
    float density   [2][2][3];
    float velocity_x[2][2][3];
    float velocity_y[2][2][3];
    float velocity_z[2][2][3];
    float pressure  [2][2][3];
};

/* ===========================================================================
 * Global-buffer access helpers
 *
 * Fortran column-major: uold(cell_idx, ivar, oct_idx) stores cell_idx fastest.
 * Flat C index: (oct_0)*NVAR*8 + (ivar_0)*8 + (cell_0)   (all 0-based).
 * ========================================================================= */
inline float u_get(device const float *u, int oct_1, int ivar_1, int cell_1) {
    return u[(oct_1-1)*(NVAR)*TWOTONDIM + (ivar_1-1)*TWOTONDIM + (cell_1-1)];
}
inline void u_set(device float *u, int oct_1, int ivar_1, int cell_1, float v) {
    u[(oct_1-1)*(NVAR)*TWOTONDIM + (ivar_1-1)*TWOTONDIM + (cell_1-1)] = v;
}
/* Flat 0-based index into u — for atomic operations that need a raw pointer. */
inline int u_flat(int oct_1, int ivar_1, int cell_1) {
    return (oct_1-1)*(NVAR)*TWOTONDIM + (ivar_1-1)*TWOTONDIM + (cell_1-1);
}
#ifdef MHD
inline float b_get(device const float *b, int oct_1, int ivar_1, int cell_1) {
    return b[(oct_1-1)*6*TWOTONDIM + (ivar_1-1)*TWOTONDIM + (cell_1-1)];
}
inline void b_set(device float *b, int oct_1, int ivar_1, int cell_1, float v) {
    b[(oct_1-1)*6*TWOTONDIM + (ivar_1-1)*TWOTONDIM + (cell_1-1)] = v;
}
inline int b_flat(int oct_1, int ivar_1, int cell_1) {
    return (oct_1-1)*6*TWOTONDIM + (ivar_1-1)*TWOTONDIM + (cell_1-1);
}
#endif

/* ===========================================================================
 * Scalar helpers — mirror gpu_hydro.cuf
 * ========================================================================= */
float magnitude_squared(float x, float y, float z) { return x*x + (y*y + z*z); }

float compute_pressure(conserved_t c, float gamma, float smallr, float smallc2) {
    float smallp = smallc2 / gamma;
    float density = max(c.density, smallr);
    float eint    = c.energy - 0.5f * magnitude_squared(c.momentum_x, c.momentum_y, c.momentum_z) / density;
#ifdef MHD
    eint -= 0.5f * magnitude_squared(c.Bx, c.By, c.Bz);
#endif
    return max((gamma - 1.0f) * eint, density * smallp);
}

float compute_energy(primitive_t p, float gamma) {
    float energy = p.pressure / (gamma - 1.0f) + 0.5f * p.density * magnitude_squared(p.velocity_x, p.velocity_y, p.velocity_z);
#ifdef MHD
    energy += 0.5f * magnitude_squared(p.Bx, p.By, p.Bz);
#endif
    return energy;
}

float sound_speed(primitive_t p, float gamma) {
    return sqrt(gamma * p.pressure / p.density);
}

primitive_t conserved_2_primitive(conserved_t c, float gamma, float smallr, float smallc2) {
    primitive_t p;
    p.density    = max(c.density, smallr);
    p.velocity_x = c.momentum_x / p.density;
    p.velocity_y = c.momentum_y / p.density;
    p.velocity_z = c.momentum_z / p.density;
    p.pressure   = compute_pressure(c, gamma, smallr, smallc2);
#ifdef MHD
    p.Bx = c.Bx;
    p.By = c.By;
    p.Bz = c.Bz;
#endif
    return p;
}

conserved_t primitive_2_conserved(primitive_t p, float gamma) {
    conserved_t c;
    c.density    = p.density;
    c.momentum_x = p.velocity_x * p.density;
    c.momentum_y = p.velocity_y * p.density;
    c.momentum_z = p.velocity_z * p.density;
    c.energy     = compute_energy(p, gamma);
#ifdef MHD
    c.Bx = p.Bx;
    c.By = p.By;
    c.Bz = p.Bz;
#endif
    return c;
}

/* slope_moncen — mirrors gpu_hydro.cuf
 * slope=0: 1st-order, slope=1: minmod, slope=2: moncen */
float slope_moncen(float left, float middle, float right, int slope) {
    float sl = middle - left;
    float sr = right  - middle;
    float sc = 0.5f * (sl + sr);
    float f  = float(slope);
    if (sl * sr <= 0.0f) return 0.0f;
    if (sl > 0.0f) return min(f * min(sl, sr), sc);
    else           return max(f * max(sl, sr), sc);
}

/* ===========================================================================
 * compute_tvd_slopes — mirrors gpu_hydro.cuf (no MHD, no scalars)
 * ========================================================================= */
void compute_tvd_slopes(threadgroup local_subgrid_t &ls, int i, int j, int k, int slope,
                        thread primitive_t &sx, thread primitive_t &sy, thread primitive_t &sz) {
    primitive_t cell;
    cell.density    = ls.density   [i][j][k];
    cell.velocity_x = ls.velocity_x[i][j][k];
    cell.velocity_y = ls.velocity_y[i][j][k];
    cell.velocity_z = ls.velocity_z[i][j][k];
    cell.pressure   = ls.pressure  [i][j][k];

    sx.density    = 0.5f * slope_moncen(ls.density   [i-1][j][k], cell.density,    ls.density   [i+1][j][k], slope);
    sx.velocity_x = 0.5f * slope_moncen(ls.velocity_x[i-1][j][k], cell.velocity_x, ls.velocity_x[i+1][j][k], slope);
    sx.velocity_y = 0.5f * slope_moncen(ls.velocity_y[i-1][j][k], cell.velocity_y, ls.velocity_y[i+1][j][k], slope);
    sx.velocity_z = 0.5f * slope_moncen(ls.velocity_z[i-1][j][k], cell.velocity_z, ls.velocity_z[i+1][j][k], slope);
    sx.pressure   = 0.5f * slope_moncen(ls.pressure  [i-1][j][k], cell.pressure,   ls.pressure  [i+1][j][k], slope);

    sy.density    = 0.5f * slope_moncen(ls.density   [i][j-1][k], cell.density,    ls.density   [i][j+1][k], slope);
    sy.velocity_x = 0.5f * slope_moncen(ls.velocity_x[i][j-1][k], cell.velocity_x, ls.velocity_x[i][j+1][k], slope);
    sy.velocity_y = 0.5f * slope_moncen(ls.velocity_y[i][j-1][k], cell.velocity_y, ls.velocity_y[i][j+1][k], slope);
    sy.velocity_z = 0.5f * slope_moncen(ls.velocity_z[i][j-1][k], cell.velocity_z, ls.velocity_z[i][j+1][k], slope);
    sy.pressure   = 0.5f * slope_moncen(ls.pressure  [i][j-1][k], cell.pressure,   ls.pressure  [i][j+1][k], slope);

    sz.density    = 0.5f * slope_moncen(ls.density   [i][j][k-1], cell.density,    ls.density   [i][j][k+1], slope);
    sz.velocity_x = 0.5f * slope_moncen(ls.velocity_x[i][j][k-1], cell.velocity_x, ls.velocity_x[i][j][k+1], slope);
    sz.velocity_y = 0.5f * slope_moncen(ls.velocity_y[i][j][k-1], cell.velocity_y, ls.velocity_y[i][j][k+1], slope);
    sz.velocity_z = 0.5f * slope_moncen(ls.velocity_z[i][j][k-1], cell.velocity_z, ls.velocity_z[i][j][k+1], slope);
    sz.pressure   = 0.5f * slope_moncen(ls.pressure  [i][j][k-1], cell.pressure,   ls.pressure  [i][j][k+1], slope);
}

/* ===========================================================================
 * subgrid_conserved_2_primitive — mirrors gpu_hydro.cuf (no MHD, no GRAV, no scalars)
 *
 * Loads a 6×6×6 stencil of primitive variables from the global uold buffer into
 * the threadgroup local_subgrid, using the nbor array for neighbour oct lookup.
 * Gravity predictor uses constant_gravity (GRAV=0 path).
 * ========================================================================= */
void subgrid_conserved_2_primitive(
    device const oct_t  *grid,
    device const float  *uold,
    device const int    *nbor,
    constant float      *constant_gravity, /* [3] */
    device const float  *f,
    int head_idx, int block_idx, int thread_idx,
    float gamma, float smallr, float smallc2, float dt,
    uint threads_per_tg,
    threadgroup local_subgrid_t &local_subgrid)
{
    const int work_size  = 2*NSUBGRID + 4; /* = 6 */
    const int total_work = work_size * work_size * work_size; /* = 216 */

    for (int work_idx = thread_idx; work_idx < total_work; work_idx += int(threads_per_tg)) {

        /* Determine which neighbour oct to load (oct-level index in 3×3×3 nbor grid) */
        int i_sg, j_sg, k_sg;
        index_1Dto3D(work_idx / 8, work_size / 2, work_size / 2, i_sg, j_sg, k_sg);

        int subgrid_idx = head_idx + block_idx;       /* 1-based oct index  */
        int ind_nbor    = work_idx / 8 + 1;           /* 1-based (1..27)    */
        int source_idx  = nbor_get(nbor, subgrid_idx, ind_nbor); /* 1-based */

        /* Which cell within the oct */
        int cell_idx = work_idx % 8 + 1;              /* 1-based (1..8)     */

        /* Position in the 6×6×6 stencil */
        int ib, jb, kb;
        index_1Dto3D(cell_idx - 1, 2, 2, ib, jb, kb);
        int i = ib + 2 * i_sg;
        int j = jb + 2 * j_sg;
        int k = kb + 2 * k_sg;

        /* Load conserved variables from global uold */
        conserved_t cv;
        cv.density    = u_get(uold, source_idx, 1, cell_idx);
        cv.momentum_x = u_get(uold, source_idx, 2, cell_idx);
        cv.momentum_y = u_get(uold, source_idx, 3, cell_idx);
        cv.momentum_z = u_get(uold, source_idx, 4, cell_idx);
        cv.energy     = u_get(uold, source_idx, 5, cell_idx);

        /* Convert to primitives */
        primitive_t pv = conserved_2_primitive(cv, gamma, smallr, smallc2);

        /* Gravity predictor step */
#ifdef GRAV
        pv.velocity_x += f[(source_idx - 1)*3*8 + 0*8 + (cell_idx - 1)] * 0.5f * dt;
        pv.velocity_y += f[(source_idx - 1)*3*8 + 1*8 + (cell_idx - 1)] * 0.5f * dt;
        pv.velocity_z += f[(source_idx - 1)*3*8 + 2*8 + (cell_idx - 1)] * 0.5f * dt;
#else
        pv.velocity_x += constant_gravity[0] * 0.5f * dt;
        pv.velocity_y += constant_gravity[1] * 0.5f * dt;
        pv.velocity_z += constant_gravity[2] * 0.5f * dt;
#endif

        /* Write to threadgroup stencil */
        local_subgrid.density   [i][j][k] = pv.density;
        local_subgrid.velocity_x[i][j][k] = pv.velocity_x;
        local_subgrid.velocity_y[i][j][k] = pv.velocity_y;
        local_subgrid.velocity_z[i][j][k] = pv.velocity_z;
        local_subgrid.pressure  [i][j][k] = pv.pressure;

        /* Store refinement flag for inner cells (Fortran 1:4 → 0-based [0..3]) */
        if (i >= 1 && i <= 2*NSUBGRID+2 && j >= 1 && j <= 2*NSUBGRID+2
                                         && k >= 1 && k <= 2*NSUBGRID+2) {
            local_subgrid.refined[i-1][j-1][k-1] =
                (grid[source_idx - 1].refined[cell_idx - 1] != 0);
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

/* ===========================================================================
 * trace_3d — mirrors gpu_hydro.cuf (no MHD, no scalars)
 *
 * MUSCL half-step predictor: computes limited slopes, applies the linearised
 * Euler source-term predictor, and writes left/right interface states.
 * ========================================================================= */
void trace_3d(
    threadgroup local_subgrid_t  &local_subgrid,
    int thread_idx, uint threads_per_tg,
    threadgroup interfaces_x_t   &left_x,  threadgroup interfaces_x_t &right_x,
    threadgroup interfaces_y_t   &left_y,  threadgroup interfaces_y_t &right_y,
    threadgroup interfaces_z_t   &left_z,  threadgroup interfaces_z_t &right_z,
    float gamma, float smallr, float smallc2, float dtdx, int slope)
{
    const int work_size  = 2*NSUBGRID + 2; /* = 4 */
    const int total_work = work_size * work_size * work_size; /* = 64 */
    float smallp = smallr * smallc2;

    for (int work_idx = thread_idx; work_idx < total_work; work_idx += int(threads_per_tg)) {
        int i, j, k;
        index_1Dto3D(work_idx, work_size, work_size, i, j, k);
        /* Shift: work_idx in 0..work_size³-1 maps to 1-based stencil indices 1..4 */
        i += 1; j += 1; k += 1;

        /* Cell-centred values */
        primitive_t cell;
        cell.density    = local_subgrid.density   [i][j][k];
        cell.velocity_x = local_subgrid.velocity_x[i][j][k];
        cell.velocity_y = local_subgrid.velocity_y[i][j][k];
        cell.velocity_z = local_subgrid.velocity_z[i][j][k];
        cell.pressure   = local_subgrid.pressure  [i][j][k];

        /* TVD slopes */
        primitive_t sx, sy, sz;
        compute_tvd_slopes(local_subgrid, i, j, k, slope, sx, sy, sz);

        /* Source terms (linearised Euler, no MHD) */
        primitive_t src;
        src.density    = -(cell.velocity_x*sx.density    + cell.velocity_y*sy.density    + cell.velocity_z*sz.density)
                         - (sx.velocity_x + sy.velocity_y + sz.velocity_z)*cell.density;
        src.velocity_x = -(cell.velocity_x*sx.velocity_x + cell.velocity_y*sy.velocity_x + cell.velocity_z*sz.velocity_x)
                         - sx.pressure / cell.density;
        src.velocity_y = -(cell.velocity_x*sx.velocity_y + cell.velocity_y*sy.velocity_y + cell.velocity_z*sz.velocity_y)
                         - sy.pressure / cell.density;
        src.velocity_z = -(cell.velocity_x*sx.velocity_z + cell.velocity_y*sy.velocity_z + cell.velocity_z*sz.velocity_z)
                         - sz.pressure / cell.density;
        src.pressure   = -(cell.velocity_x*sx.pressure   + cell.velocity_y*sy.pressure   + cell.velocity_z*sz.pressure)
                         - (sx.velocity_x + sy.velocity_y + sz.velocity_z)*gamma*cell.pressure;

        /* Half-step predicted state */
        cell.density    += dtdx * src.density;
        cell.velocity_x += dtdx * src.velocity_x;
        cell.velocity_y += dtdx * src.velocity_y;
        cell.velocity_z += dtdx * src.velocity_z;
        cell.pressure   += dtdx * src.pressure;

        /* Write interface states.  Index arithmetic matches Fortran (i-2, j-2, k-2) → 0-based. */

        /* Right state at i-1/2 */
        if (i > 1 && (j > 1 && j < work_size) && (k > 1 && k < work_size)) {
            int ri=i-2, rj=j-2, rk=k-2;
            right_x.density   [ri][rj][rk] = cell.density    - sx.density;
            right_x.velocity_x[ri][rj][rk] = cell.velocity_x - sx.velocity_x;
            right_x.velocity_y[ri][rj][rk] = cell.velocity_y - sx.velocity_y;
            right_x.velocity_z[ri][rj][rk] = cell.velocity_z - sx.velocity_z;
            right_x.pressure  [ri][rj][rk] = cell.pressure   - sx.pressure;
            if (right_x.density [ri][rj][rk] < smallr) right_x.density [ri][rj][rk] = local_subgrid.density [i][j][k];
            if (right_x.pressure[ri][rj][rk] < smallp) right_x.pressure[ri][rj][rk] = local_subgrid.pressure[i][j][k];
        }
        /* Left state at i+1/2 */
        if (i < work_size && (j > 1 && j < work_size) && (k > 1 && k < work_size)) {
            int li=i-1, lj=j-2, lk=k-2;
            left_x.density   [li][lj][lk] = cell.density    + sx.density;
            left_x.velocity_x[li][lj][lk] = cell.velocity_x + sx.velocity_x;
            left_x.velocity_y[li][lj][lk] = cell.velocity_y + sx.velocity_y;
            left_x.velocity_z[li][lj][lk] = cell.velocity_z + sx.velocity_z;
            left_x.pressure  [li][lj][lk] = cell.pressure   + sx.pressure;
            if (left_x.density [li][lj][lk] < smallr) left_x.density [li][lj][lk] = local_subgrid.density [i][j][k];
            if (left_x.pressure[li][lj][lk] < smallp) left_x.pressure[li][lj][lk] = local_subgrid.pressure[i][j][k];
        }
        /* Right state at j-1/2 */
        if ((i > 1 && i < work_size) && j > 1 && (k > 1 && k < work_size)) {
            int ri=i-2, rj=j-2, rk=k-2;
            right_y.density   [ri][rj][rk] = cell.density    - sy.density;
            right_y.velocity_x[ri][rj][rk] = cell.velocity_x - sy.velocity_x;
            right_y.velocity_y[ri][rj][rk] = cell.velocity_y - sy.velocity_y;
            right_y.velocity_z[ri][rj][rk] = cell.velocity_z - sy.velocity_z;
            right_y.pressure  [ri][rj][rk] = cell.pressure   - sy.pressure;
            if (right_y.density [ri][rj][rk] < smallr) right_y.density [ri][rj][rk] = local_subgrid.density [i][j][k];
            if (right_y.pressure[ri][rj][rk] < smallp) right_y.pressure[ri][rj][rk] = local_subgrid.pressure[i][j][k];
        }
        /* Left state at j+1/2 */
        if ((i > 1 && i < work_size) && j < work_size && (k > 1 && k < work_size)) {
            int li=i-2, lj=j-1, lk=k-2;
            left_y.density   [li][lj][lk] = cell.density    + sy.density;
            left_y.velocity_x[li][lj][lk] = cell.velocity_x + sy.velocity_x;
            left_y.velocity_y[li][lj][lk] = cell.velocity_y + sy.velocity_y;
            left_y.velocity_z[li][lj][lk] = cell.velocity_z + sy.velocity_z;
            left_y.pressure  [li][lj][lk] = cell.pressure   + sy.pressure;
            if (left_y.density [li][lj][lk] < smallr) left_y.density [li][lj][lk] = local_subgrid.density [i][j][k];
            if (left_y.pressure[li][lj][lk] < smallp) left_y.pressure[li][lj][lk] = local_subgrid.pressure[i][j][k];
        }
        /* Right state at k-1/2 */
        if ((i > 1 && i < work_size) && (j > 1 && j < work_size) && k > 1) {
            int ri=i-2, rj=j-2, rk=k-2;
            right_z.density   [ri][rj][rk] = cell.density    - sz.density;
            right_z.velocity_x[ri][rj][rk] = cell.velocity_x - sz.velocity_x;
            right_z.velocity_y[ri][rj][rk] = cell.velocity_y - sz.velocity_y;
            right_z.velocity_z[ri][rj][rk] = cell.velocity_z - sz.velocity_z;
            right_z.pressure  [ri][rj][rk] = cell.pressure   - sz.pressure;
            if (right_z.density [ri][rj][rk] < smallr) right_z.density [ri][rj][rk] = local_subgrid.density [i][j][k];
            if (right_z.pressure[ri][rj][rk] < smallp) right_z.pressure[ri][rj][rk] = local_subgrid.pressure[i][j][k];
        }
        /* Left state at k+1/2 */
        if ((i > 1 && i < work_size) && (j > 1 && j < work_size) && k < work_size) {
            int li=i-2, lj=j-2, lk=k-1;
            left_z.density   [li][lj][lk] = cell.density    + sz.density;
            left_z.velocity_x[li][lj][lk] = cell.velocity_x + sz.velocity_x;
            left_z.velocity_y[li][lj][lk] = cell.velocity_y + sz.velocity_y;
            left_z.velocity_z[li][lj][lk] = cell.velocity_z + sz.velocity_z;
            left_z.pressure  [li][lj][lk] = cell.pressure   + sz.pressure;
            if (left_z.density [li][lj][lk] < smallr) left_z.density [li][lj][lk] = local_subgrid.density [i][j][k];
            if (left_z.pressure[li][lj][lk] < smallp) left_z.pressure[li][lj][lk] = local_subgrid.pressure[i][j][k];
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

/* ===========================================================================
 * Riemann solvers — mirror gpu_hydro.cuf (no MHD)
 * ========================================================================= */

/* hll_flux — single-field HLL flux */
float hll_flux(float sl, float sr, float fl, float fr, float cl, float cr) {
    return (sr*fl - sl*fr + sr*sl*(cr - cl)) / (sr - sl);
}

/* hll_fluxes — mirrors gpu_hydro.cuf; llf=true → LLF wave speeds */
conserved_t hll_fluxes(thread primitive_t &left, thread primitive_t &right,
                       float gamma, float smallr, float smallc2, bool llf) {
    float smallp = smallc2 / gamma;
    left.density   = max(left.density,   smallr);
    right.density  = max(right.density,  smallr);
    left.pressure  = max(left.pressure,  smallp * left.density);
    right.pressure = max(right.pressure, smallp * right.density);

    float cl = sound_speed(left,  gamma);
    float cr = sound_speed(right, gamma);
    float sl, sr;
    if (llf) {
        float spd = max(abs(left.velocity_x) + cl, abs(right.velocity_x) + cr);
        sl = -spd; sr = spd;
    } else {
        sl = min(min(left.velocity_x, right.velocity_x) - max(cl, cr), 0.0f);
        sr = max(max(left.velocity_x, right.velocity_x) + max(cl, cr), 0.0f);
    }

    conserved_t lc = primitive_2_conserved(left,  gamma);
    conserved_t rc = primitive_2_conserved(right, gamma);

    conserved_t lf, rf;
    lf.density    = lc.momentum_x;
    lf.momentum_x = left.pressure  + left.velocity_x  * lc.momentum_x;
    lf.momentum_y = left.velocity_x  * lc.momentum_y;
    lf.momentum_z = left.velocity_x  * lc.momentum_z;
    lf.energy     = left.velocity_x  * (left.pressure  + lc.energy);
    rf.density    = rc.momentum_x;
    rf.momentum_x = right.pressure + right.velocity_x * rc.momentum_x;
    rf.momentum_y = right.velocity_x * rc.momentum_y;
    rf.momentum_z = right.velocity_x * rc.momentum_z;
    rf.energy     = right.velocity_x * (right.pressure + rc.energy);

    conserved_t flux;
    flux.density    = hll_flux(sl, sr, lf.density,    rf.density,    lc.density,    rc.density);
    flux.momentum_x = hll_flux(sl, sr, lf.momentum_x, rf.momentum_x, lc.momentum_x, rc.momentum_x);
    flux.momentum_y = hll_flux(sl, sr, lf.momentum_y, rf.momentum_y, lc.momentum_y, rc.momentum_y);
    flux.momentum_z = hll_flux(sl, sr, lf.momentum_z, rf.momentum_z, lc.momentum_z, rc.momentum_z);
    flux.energy     = hll_flux(sl, sr, lf.energy,     rf.energy,     lc.energy,     rc.energy);
    return flux;
}

/* hllc_fluxes — mirrors gpu_hydro.cuf */
conserved_t hllc_fluxes(thread primitive_t &left, thread primitive_t &right,
                         float gamma, float smallr, float smallc2) {
    float smallp = smallc2 / gamma;
    left.density   = max(left.density,   smallr);
    right.density  = max(right.density,  smallr);
    left.pressure  = max(left.pressure,  smallp * left.density);
    right.pressure = max(right.pressure, smallp * right.density);

    float rl=left.density,  rr=right.density;
    float ul=left.velocity_x,  ur=right.velocity_x;
    float vl=left.velocity_y,  vr=right.velocity_y;
    float wl=left.velocity_z,  wr=right.velocity_z;
    float pl=left.pressure,    pr=right.pressure;

    float entho = 1.0f / (gamma - 1.0f);
    float etotl = pl*entho + 0.5f*rl*magnitude_squared(ul, vl, wl);
    float etotr = pr*entho + 0.5f*rr*magnitude_squared(ur, vr, wr);

    float cl = sound_speed(left,  gamma);
    float cr = sound_speed(right, gamma);
    float sl = min(ul, ur) - max(cl, cr);
    float sr = max(ul, ur) + max(cl, cr);

    float rcl   = rl*(ul - sl);
    float rcr   = rr*(sr - ur);
    float ustar = (rcr*ur + rcl*ul + (pl - pr)) / (rcr + rcl);
    float pstar = (rcr*pl + rcl*pr + rcl*rcr*(ul - ur)) / (rcr + rcl);
    float rstarl = rl*(sl - ul) / (sl - ustar);
    float rstarr = rr*(sr - ur) / (sr - ustar);
    float el_star = ((sl - ul)*etotl - pl*ul + pstar*ustar) / (sl - ustar);
    float er_star = ((sr - ur)*etotr - pr*ur + pstar*ustar) / (sr - ustar);

    float ro, uo, vo, wo, po, etoto;
    if      (sl    > 0.0f) { ro=rl;     uo=ul;    vo=vl; wo=wl; po=pl;    etoto=etotl;   }
    else if (ustar > 0.0f) { ro=rstarl; uo=ustar; vo=vl; wo=wl; po=pstar; etoto=el_star; }
    else if (sr    > 0.0f) { ro=rstarr; uo=ustar; vo=vr; wo=wr; po=pstar; etoto=er_star; }
    else                   { ro=rr;     uo=ur;    vo=vr; wo=wr; po=pr;    etoto=etotr;   }

    conserved_t flux;
    flux.density    = ro*uo;
    flux.momentum_x = ro*uo*uo + po;
    flux.momentum_y = ro*uo*vo;
    flux.momentum_z = ro*uo*wo;
    flux.energy     = (etoto + po)*uo;
    return flux;
}

/* riemann_fluxes — dispatcher, mirrors gpu_hydro.cuf (no MHD) */
conserved_t riemann_fluxes(thread primitive_t &left, thread primitive_t &right,
                            float gamma, float smallr, float smallc2, int riemann) {
    if      (riemann == SOLVER_HLLC) return hllc_fluxes(left, right, gamma, smallr, smallc2);
    else if (riemann == SOLVER_HLL)  return hll_fluxes (left, right, gamma, smallr, smallc2, false);
    else if (riemann == SOLVER_LLF)  return hll_fluxes (left, right, gamma, smallr, smallc2, true);
    else                             return hll_fluxes (left, right, gamma, smallr, smallc2, true); /* LLF */
}

/* ===========================================================================
 * riemann_driver — mirrors gpu_hydro.cuf (no MHD, no scalars)
 *
 * Each thread owns one or more interfaces (loop with stride = blockDim).
 * For y/z, velocity components are cyclically rotated into the x-normal frame
 * before the Riemann call and rotated back when storing fluxes.
 * The left_interfaces_* arrays are reused to hold the computed conserved fluxes.
 * ========================================================================= */
void riemann_driver(
    threadgroup interfaces_x_t &left_x,  threadgroup interfaces_x_t &right_x,
    threadgroup interfaces_y_t &left_y,  threadgroup interfaces_y_t &right_y,
    threadgroup interfaces_z_t &left_z,  threadgroup interfaces_z_t &right_z,
    int thread_idx, uint threads_per_tg,
    float gamma, float smallr, float smallc2, int riemann)
{
    /* interface_array_size = (2*NSUBGRID+1)*(2*NSUBGRID)^2 = 3*4 = 12 for nsubgrid=1 */
    const int ias = (2*NSUBGRID+1) * (2*NSUBGRID) * (2*NSUBGRID);

    for (int work_idx = thread_idx; work_idx < 3*ias; work_idx += int(threads_per_tg)) {
        int i, j, k;
        primitive_t L, R;
        conserved_t flux;

        if (work_idx < ias) {
            /* ---- x-interfaces ---- */
            index_1Dto3D(work_idx, 2*NSUBGRID+1, 2*NSUBGRID, i, j, k);
            L.density=left_x.density[i][j][k];  L.velocity_x=left_x.velocity_x[i][j][k];
            L.velocity_y=left_x.velocity_y[i][j][k]; L.velocity_z=left_x.velocity_z[i][j][k];
            L.pressure=left_x.pressure[i][j][k];
            R.density=right_x.density[i][j][k]; R.velocity_x=right_x.velocity_x[i][j][k];
            R.velocity_y=right_x.velocity_y[i][j][k]; R.velocity_z=right_x.velocity_z[i][j][k];
            R.pressure=right_x.pressure[i][j][k];

            flux = riemann_fluxes(L, R, gamma, smallr, smallc2, riemann);

            left_x.density   [i][j][k] = flux.density;
            left_x.velocity_x[i][j][k] = flux.momentum_x;
            left_x.velocity_y[i][j][k] = flux.momentum_y;
            left_x.velocity_z[i][j][k] = flux.momentum_z;
            left_x.pressure  [i][j][k] = flux.energy;

        } else if (work_idx < 2*ias) {
            /* ---- y-interfaces: rotate velocity (vy,vz,vx) → (vx,vy,vz) ---- */
            index_1Dto3D(work_idx - ias, 2*NSUBGRID, 2*NSUBGRID+1, i, j, k);
            L.density=left_y.density[i][j][k];  L.pressure=left_y.pressure[i][j][k];
            L.velocity_x=left_y.velocity_y[i][j][k]; /* vy becomes normal velocity */
            L.velocity_y=left_y.velocity_z[i][j][k];
            L.velocity_z=left_y.velocity_x[i][j][k];
            R.density=right_y.density[i][j][k]; R.pressure=right_y.pressure[i][j][k];
            R.velocity_x=right_y.velocity_y[i][j][k];
            R.velocity_y=right_y.velocity_z[i][j][k];
            R.velocity_z=right_y.velocity_x[i][j][k];

            flux = riemann_fluxes(L, R, gamma, smallr, smallc2, riemann);

            /* Rotate flux back: mx→vy slot, my→vz slot, mz→vx slot */
            left_y.density   [i][j][k] = flux.density;
            left_y.velocity_y[i][j][k] = flux.momentum_x;
            left_y.velocity_z[i][j][k] = flux.momentum_y;
            left_y.velocity_x[i][j][k] = flux.momentum_z;
            left_y.pressure  [i][j][k] = flux.energy;

        } else {
            /* ---- z-interfaces: rotate velocity (vz,vx,vy) → (vx,vy,vz) ---- */
            index_1Dto3D(work_idx - 2*ias, 2*NSUBGRID, 2*NSUBGRID, i, j, k);
            L.density=left_z.density[i][j][k];  L.pressure=left_z.pressure[i][j][k];
            L.velocity_x=left_z.velocity_z[i][j][k]; /* vz becomes normal velocity */
            L.velocity_y=left_z.velocity_x[i][j][k];
            L.velocity_z=left_z.velocity_y[i][j][k];
            R.density=right_z.density[i][j][k]; R.pressure=right_z.pressure[i][j][k];
            R.velocity_x=right_z.velocity_z[i][j][k];
            R.velocity_y=right_z.velocity_x[i][j][k];
            R.velocity_z=right_z.velocity_y[i][j][k];

            flux = riemann_fluxes(L, R, gamma, smallr, smallc2, riemann);

            /* Rotate flux back: mx→vz slot, my→vx slot, mz→vy slot */
            left_z.density   [i][j][k] = flux.density;
            left_z.velocity_z[i][j][k] = flux.momentum_x;
            left_z.velocity_x[i][j][k] = flux.momentum_y;
            left_z.velocity_y[i][j][k] = flux.momentum_z;
            left_z.pressure  [i][j][k] = flux.energy;
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

/* ===========================================================================
 * zero_fine_fluxes — mirrors gpu_hydro.cuf (no MHD, no scalars).
 *
 * For each interface between two cells where either is refined, zero the flux
 * so fine-fine fluxes do not contaminate the coarse update.
 * Called when ilevel < levelmax.
 * ========================================================================= */
void zero_fine_fluxes(
    threadgroup const local_subgrid_t &ls,
    int thread_idx, uint threads_per_tg,
    threadgroup interfaces_x_t &fluxes_x,
    threadgroup interfaces_y_t &fluxes_y,
    threadgroup interfaces_z_t &fluxes_z)
{
    const int ias = (2*NSUBGRID+1) * (2*NSUBGRID) * (2*NSUBGRID); /* =12 for nsubgrid=1 */

    for (int work_idx = thread_idx; work_idx < 3 * ias; work_idx += int(threads_per_tg)) {
        int i, j, k;
        if (work_idx < ias) {
            /* X interfaces */
            index_1Dto3D(work_idx, 2*NSUBGRID+1, 2*NSUBGRID, i, j, k);
            /* Fortran: refined(i+1,j+2,k+2) or refined(i+2,j+2,k+2)
             * Metal 0-based: refined[i][j+1][k+1] or refined[i+1][j+1][k+1] */
            if (ls.refined[i][j+1][k+1] || ls.refined[i+1][j+1][k+1]) {
                fluxes_x.density   [i][j][k] = 0.0f;
                fluxes_x.velocity_x[i][j][k] = 0.0f;
                fluxes_x.velocity_y[i][j][k] = 0.0f;
                fluxes_x.velocity_z[i][j][k] = 0.0f;
                fluxes_x.pressure  [i][j][k] = 0.0f;
            }
        } else if (work_idx < 2 * ias) {
            /* Y interfaces */
            index_1Dto3D(work_idx - ias, 2*NSUBGRID, 2*NSUBGRID+1, i, j, k);
            /* Fortran: refined(i+2,j+1,k+2) or refined(i+2,j+2,k+2)
             * Metal 0-based: refined[i+1][j][k+1] or refined[i+1][j+1][k+1] */
            if (ls.refined[i+1][j][k+1] || ls.refined[i+1][j+1][k+1]) {
                fluxes_y.density   [i][j][k] = 0.0f;
                fluxes_y.velocity_x[i][j][k] = 0.0f;
                fluxes_y.velocity_y[i][j][k] = 0.0f;
                fluxes_y.velocity_z[i][j][k] = 0.0f;
                fluxes_y.pressure  [i][j][k] = 0.0f;
            }
        } else {
            /* Z interfaces */
            index_1Dto3D(work_idx - 2*ias, 2*NSUBGRID, 2*NSUBGRID, i, j, k);
            /* Fortran: refined(i+2,j+2,k+1) or refined(i+2,j+2,k+2)
             * Metal 0-based: refined[i+1][j+1][k] or refined[i+1][j+1][k+1] */
            if (ls.refined[i+1][j+1][k] || ls.refined[i+1][j+1][k+1]) {
                fluxes_z.density   [i][j][k] = 0.0f;
                fluxes_z.velocity_x[i][j][k] = 0.0f;
                fluxes_z.velocity_y[i][j][k] = 0.0f;
                fluxes_z.velocity_z[i][j][k] = 0.0f;
                fluxes_z.pressure  [i][j][k] = 0.0f;
            }
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

/* ===========================================================================
 * coarse_cell_update — mirrors gpu_hydro.cuf (no MHD, no scalars).
 *
 * Applies boundary flux corrections to coarser-level cells in unew.
 * For nsubgrid=1: interface_array_size=1 → exactly 6 active threads (one per face).
 * The correction reaches ghost-cache octs (source_idx > ngridmax): their
 * father oct is the actual coarse grid oct to be corrected.
 * Called when ilevel > levelmin.
 * ========================================================================= */
void coarse_cell_update(
    device float              *unew,
    device const int          *nbor,
    device const oct_t        *grid,
    device const int          *father,
    threadgroup const interfaces_x_t &fluxes_x,
    threadgroup const interfaces_y_t &fluxes_y,
    threadgroup const interfaces_z_t &fluxes_z,
    int head_idx, int block_idx, int ngridmax,
    int thread_idx, float coarse_flux_scale)
{
    /* For nsubgrid=1: interface_array_size = nsubgrid^(ndim-1) = 1 */
    const int ias = NSUBGRID * NSUBGRID;   /* =1 for nsubgrid=1 */
    if (thread_idx >= 6 * ias) return;

    int subgrid_idx = head_idx + block_idx;   /* 1-based */
    int face    = thread_idx / ias;
    int work_idx = thread_idx % ias;

    /* Decode 2D position within the face (for nsubgrid=1, always 0) */
    int j_raw = work_idx % NSUBGRID;
    int k_raw = work_idx / NSUBGRID;

    float fd=0.0f, fmx=0.0f, fmy=0.0f, fmz=0.0f, fe=0.0f;
    int ind_nbor = 0;

    if (face == 0) {
        /* Left X: i_sg=0, j_sg=j_raw+1, k_sg=k_raw+1 */
        int j_sg = j_raw + 1, k_sg = k_raw + 1;
        for (int j = 2*j_sg-2; j <= 2*j_sg-1; j++)
            for (int k = 2*k_sg-2; k <= 2*k_sg-1; k++) {
                fd  -= fluxes_x.density   [0][j][k] * coarse_flux_scale;
                fmx -= fluxes_x.velocity_x[0][j][k] * coarse_flux_scale;
                fmy -= fluxes_x.velocity_y[0][j][k] * coarse_flux_scale;
                fmz -= fluxes_x.velocity_z[0][j][k] * coarse_flux_scale;
                fe  -= fluxes_x.pressure  [0][j][k] * coarse_flux_scale;
            }
        ind_nbor = 1 + 0 + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*k_sg;
    } else if (face == 1) {
        /* Right X: i_sg=NSUBGRID+1, j_sg=j_raw+1, k_sg=k_raw+1 */
        int j_sg = j_raw + 1, k_sg = k_raw + 1;
        for (int j = 2*j_sg-2; j <= 2*j_sg-1; j++)
            for (int k = 2*k_sg-2; k <= 2*k_sg-1; k++) {
                fd  += fluxes_x.density   [2*NSUBGRID][j][k] * coarse_flux_scale;
                fmx += fluxes_x.velocity_x[2*NSUBGRID][j][k] * coarse_flux_scale;
                fmy += fluxes_x.velocity_y[2*NSUBGRID][j][k] * coarse_flux_scale;
                fmz += fluxes_x.velocity_z[2*NSUBGRID][j][k] * coarse_flux_scale;
                fe  += fluxes_x.pressure  [2*NSUBGRID][j][k] * coarse_flux_scale;
            }
        ind_nbor = 1 + (NSUBGRID+1) + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*k_sg;
    } else if (face == 2) {
        /* Left Y: i_sg=j_raw+1, j_sg=0, k_sg=k_raw+1 */
        int i_sg = j_raw + 1, k_sg = k_raw + 1;
        for (int i = 2*i_sg-2; i <= 2*i_sg-1; i++)
            for (int k = 2*k_sg-2; k <= 2*k_sg-1; k++) {
                fd  -= fluxes_y.density   [i][0][k] * coarse_flux_scale;
                fmx -= fluxes_y.velocity_x[i][0][k] * coarse_flux_scale;
                fmy -= fluxes_y.velocity_y[i][0][k] * coarse_flux_scale;
                fmz -= fluxes_y.velocity_z[i][0][k] * coarse_flux_scale;
                fe  -= fluxes_y.pressure  [i][0][k] * coarse_flux_scale;
            }
        ind_nbor = 1 + i_sg + NSUBGRIDP2*0 + NSUBGRIDP2*NSUBGRIDP2*k_sg;
    } else if (face == 3) {
        /* Right Y: i_sg=j_raw+1, j_sg=NSUBGRID+1, k_sg=k_raw+1 */
        int i_sg = j_raw + 1, k_sg = k_raw + 1;
        for (int i = 2*i_sg-2; i <= 2*i_sg-1; i++)
            for (int k = 2*k_sg-2; k <= 2*k_sg-1; k++) {
                fd  += fluxes_y.density   [i][2*NSUBGRID][k] * coarse_flux_scale;
                fmx += fluxes_y.velocity_x[i][2*NSUBGRID][k] * coarse_flux_scale;
                fmy += fluxes_y.velocity_y[i][2*NSUBGRID][k] * coarse_flux_scale;
                fmz += fluxes_y.velocity_z[i][2*NSUBGRID][k] * coarse_flux_scale;
                fe  += fluxes_y.pressure  [i][2*NSUBGRID][k] * coarse_flux_scale;
            }
        ind_nbor = 1 + i_sg + NSUBGRIDP2*(NSUBGRID+1) + NSUBGRIDP2*NSUBGRIDP2*k_sg;
    } else if (face == 4) {
        /* Left Z: i_sg=j_raw+1, j_sg=k_raw+1, k_sg=0 */
        int i_sg = j_raw + 1, j_sg = k_raw + 1;
        for (int i = 2*i_sg-2; i <= 2*i_sg-1; i++)
            for (int j = 2*j_sg-2; j <= 2*j_sg-1; j++) {
                fd  -= fluxes_z.density   [i][j][0] * coarse_flux_scale;
                fmx -= fluxes_z.velocity_x[i][j][0] * coarse_flux_scale;
                fmy -= fluxes_z.velocity_y[i][j][0] * coarse_flux_scale;
                fmz -= fluxes_z.velocity_z[i][j][0] * coarse_flux_scale;
                fe  -= fluxes_z.pressure  [i][j][0] * coarse_flux_scale;
            }
        ind_nbor = 1 + i_sg + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*0;
    } else {
        /* Right Z: i_sg=j_raw+1, j_sg=k_raw+1, k_sg=NSUBGRID+1 */
        int i_sg = j_raw + 1, j_sg = k_raw + 1;
        for (int i = 2*i_sg-2; i <= 2*i_sg-1; i++)
            for (int j = 2*j_sg-2; j <= 2*j_sg-1; j++) {
                fd  += fluxes_z.density   [i][j][2*NSUBGRID] * coarse_flux_scale;
                fmx += fluxes_z.velocity_x[i][j][2*NSUBGRID] * coarse_flux_scale;
                fmy += fluxes_z.velocity_y[i][j][2*NSUBGRID] * coarse_flux_scale;
                fmz += fluxes_z.velocity_z[i][j][2*NSUBGRID] * coarse_flux_scale;
                fe  += fluxes_z.pressure  [i][j][2*NSUBGRID] * coarse_flux_scale;
            }
        ind_nbor = 1 + i_sg + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*(NSUBGRID+1);
    }

    /* Condition: source is ghost-cache oct → its father is the actual coarse oct */
    int source_idx = nbor_get(nbor, subgrid_idx, ind_nbor);
    if (source_idx > ngridmax) {
        int father_idx = father[source_idx - 1];
        int ic = grid[source_idx - 1].ckey[0] - 2 * grid[father_idx - 1].ckey[0];
        int jc = grid[source_idx - 1].ckey[1] - 2 * grid[father_idx - 1].ckey[1];
        int kc = grid[source_idx - 1].ckey[2] - 2 * grid[father_idx - 1].ckey[2];
        int cell_idx = 1 + ic + 2*jc + 4*kc;
        atomic_add_float((device atomic_uint*)&unew[u_flat(father_idx,1,cell_idx)], fd);
        atomic_add_float((device atomic_uint*)&unew[u_flat(father_idx,2,cell_idx)], fmx);
        atomic_add_float((device atomic_uint*)&unew[u_flat(father_idx,3,cell_idx)], fmy);
        atomic_add_float((device atomic_uint*)&unew[u_flat(father_idx,4,cell_idx)], fmz);
        atomic_add_float((device atomic_uint*)&unew[u_flat(father_idx,5,cell_idx)], fe);
    }
}

/* ===========================================================================
 * conservative_update — mirrors gpu_hydro.cuf (no MHD, no scalars)
 *
 * Applies flux divergence to unew for the 2×2×2 inner cells of the subgrid.
 * For nsubgrid=1 all 8 work items write to the central oct (nbor slot 14).
 * ========================================================================= */
void conservative_update(
    device float              *unew,
    device const int          *nbor,
    threadgroup const interfaces_x_t &fluxes_x,
    threadgroup const interfaces_y_t &fluxes_y,
    threadgroup const interfaces_z_t &fluxes_z,
    int head_idx, int block_idx,
    int thread_idx, uint threads_per_tg,
    float dtdx)
{
    const int work_size  = 2*NSUBGRID;         /* = 2 */
    const int total_work = work_size*work_size*work_size; /* = 8 */

    for (int work_idx = thread_idx; work_idx < total_work; work_idx += int(threads_per_tg)) {
        int subgrid_idx = head_idx + block_idx; /* 1-based */

        /* Oct-level index within the inner (nsubgrid×nsubgrid×nsubgrid) subgrid */
        int i_sg, j_sg, k_sg;
        index_1Dto3D(work_idx / 8, work_size/2, work_size/2, i_sg, j_sg, k_sg);
        i_sg += 1; j_sg += 1; k_sg += 1;

        /* For nsubgrid=1: ind_nbor=14 (central oct = the subgrid itself) */
        int ind_nbor = 1 + i_sg + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*k_sg;
        int oct_idx  = nbor_get(nbor, subgrid_idx, ind_nbor); /* 1-based */

        int cell_idx = work_idx % 8 + 1;       /* 1-based (1..8) */
        int i, j, k;
        index_1Dto3D(cell_idx - 1, 2, 2, i, j, k);
        i += 2*(i_sg - 1);
        j += 2*(j_sg - 1);
        k += 2*(k_sg - 1);

        /* Flux divergence: (F_i − F_{i+1})*dtdx — same formula as CUDA */
        float upd_rho = (fluxes_x.density   [i][j][k] - fluxes_x.density   [i+1][j  ][k  ]) * dtdx
                      + (fluxes_y.density   [i][j][k] - fluxes_y.density   [i  ][j+1][k  ]) * dtdx
                      + (fluxes_z.density   [i][j][k] - fluxes_z.density   [i  ][j  ][k+1]) * dtdx;
        float upd_mx  = (fluxes_x.velocity_x[i][j][k] - fluxes_x.velocity_x[i+1][j  ][k  ]) * dtdx
                      + (fluxes_y.velocity_x[i][j][k] - fluxes_y.velocity_x[i  ][j+1][k  ]) * dtdx
                      + (fluxes_z.velocity_x[i][j][k] - fluxes_z.velocity_x[i  ][j  ][k+1]) * dtdx;
        float upd_my  = (fluxes_x.velocity_y[i][j][k] - fluxes_x.velocity_y[i+1][j  ][k  ]) * dtdx
                      + (fluxes_y.velocity_y[i][j][k] - fluxes_y.velocity_y[i  ][j+1][k  ]) * dtdx
                      + (fluxes_z.velocity_y[i][j][k] - fluxes_z.velocity_y[i  ][j  ][k+1]) * dtdx;
        float upd_mz  = (fluxes_x.velocity_z[i][j][k] - fluxes_x.velocity_z[i+1][j  ][k  ]) * dtdx
                      + (fluxes_y.velocity_z[i][j][k] - fluxes_y.velocity_z[i  ][j+1][k  ]) * dtdx
                      + (fluxes_z.velocity_z[i][j][k] - fluxes_z.velocity_z[i  ][j  ][k+1]) * dtdx;
        float upd_e   = (fluxes_x.pressure  [i][j][k] - fluxes_x.pressure  [i+1][j  ][k  ]) * dtdx
                      + (fluxes_y.pressure  [i][j][k] - fluxes_y.pressure  [i  ][j+1][k  ]) * dtdx
                      + (fluxes_z.pressure  [i][j][k] - fluxes_z.pressure  [i  ][j  ][k+1]) * dtdx;

        u_set(unew, oct_idx, 1, cell_idx, u_get(unew, oct_idx, 1, cell_idx) + upd_rho);
        u_set(unew, oct_idx, 2, cell_idx, u_get(unew, oct_idx, 2, cell_idx) + upd_mx);
        u_set(unew, oct_idx, 3, cell_idx, u_get(unew, oct_idx, 3, cell_idx) + upd_my);
        u_set(unew, oct_idx, 4, cell_idx, u_get(unew, oct_idx, 4, cell_idx) + upd_mz);
        u_set(unew, oct_idx, 5, cell_idx, u_get(unew, oct_idx, 5, cell_idx) + upd_e);
    }
}

/* ===========================================================================
 * set_uold_kernel — mirrors attributes(global) set_uold_kernel in gpu_hydro.cuf
 *
 * Copies unew → uold for all octs at ilevel.
 * Thread layout mirrors CUDA dim3(8,16,1): threadIdx.x=cell_idx (0..7),
 * threadIdx.y=oct offset within block (0..15).
 * ========================================================================= */
kernel void set_uold_kernel(
    device float        *uold      [[buffer(0)]],
    device const float  *unew      [[buffer(1)]],
#ifdef MHD
    device float        *bold      [[buffer(2)]],
    device const float  *bnew      [[buffer(3)]],
    constant int        &head_idx  [[buffer(4)]],
    constant int        &num_octs  [[buffer(5)]],
#else
    constant int        &head_idx  [[buffer(2)]],
    constant int        &num_octs  [[buffer(3)]],
#endif
    uint2 tptg  [[thread_position_in_threadgroup]],
    uint2 blk   [[threadgroup_position_in_grid]])
{
    int oct_offset = int(blk.x) * NTHREADS_Y + int(tptg.y);
    if (oct_offset >= num_octs) return;
    int oct_idx  = head_idx + oct_offset; /* 1-based */
    int cell_idx = int(tptg.x) + 1;      /* 1-based (1..8) */

    for (int ivar = 1; ivar <= NVAR; ++ivar) u_set(uold, oct_idx, ivar, cell_idx, u_get(unew, oct_idx, ivar, cell_idx));
#ifdef MHD
    for (int ivar = 1; ivar <= 6; ++ivar) b_set(bold, oct_idx, ivar, cell_idx, b_get(bnew, oct_idx, ivar, cell_idx));
#endif
}

/* ===========================================================================
 * set_unew_kernel — mirrors attributes(global) set_unew_kernel in gpu_hydro.cuf
 *
 * Copies uold → unew for all octs at ilevel.
 * ========================================================================= */
kernel void set_unew_kernel(
    device const float  *uold      [[buffer(0)]],
    device float        *unew      [[buffer(1)]],
#ifdef MHD
    device const float  *bold      [[buffer(2)]],
    device float        *bnew      [[buffer(3)]],
    constant int        &head_idx  [[buffer(4)]],
    constant int        &num_octs  [[buffer(5)]],
#else
    constant int        &head_idx  [[buffer(2)]],
    constant int        &num_octs  [[buffer(3)]],
#endif
    uint2 tptg  [[thread_position_in_threadgroup]],
    uint2 blk   [[threadgroup_position_in_grid]])
{
    int oct_offset = int(blk.x) * NTHREADS_Y + int(tptg.y);
    if (oct_offset >= num_octs) return;
    int oct_idx  = head_idx + oct_offset;
    int cell_idx = int(tptg.x) + 1;

    for (int ivar = 1; ivar <= NVAR; ++ivar) u_set(unew, oct_idx, ivar, cell_idx, u_get(uold, oct_idx, ivar, cell_idx));
#ifdef MHD
    for (int ivar = 1; ivar <= 6; ++ivar) b_set(bnew, oct_idx, ivar, cell_idx, b_get(bold, oct_idx, ivar, cell_idx));
#endif
}

/* ===========================================================================
 * cmpdt_kernel — mirrors attributes(global) cmpdt_kernel in gpu_hydro.cuf
 *
 * 1024 threads/threadgroup, 32-wide SIMD groups → 32 SIMD groups per TG.
 * Each thread owns one cell (8 cells/oct × ceil(num_octs*8/1024) TGs).
 * Two-level reduction: simd_min/simd_sum within each SIMD group → tg_*[32],
 * then a second simd_min/simd_sum by SIMD-group 0 → one atomic per TG.
 * With 32 SIMD groups and 32 lanes in group 0, the second level is fully
 * loaded with no padding, and global atomics are reduced 4× vs 256-thread.
 *
 * data_buf layout: atomic_uint[5]
 *   [0..3] mass/ekin/eint/emag — fp32 via CAS atomic_add_float
 *   [4]    dt                  — fp32 min via uint bit-cast atomic_min_float_bits
 * ========================================================================= */
kernel void cmpdt_kernel(
    device const oct_t  *grid             [[buffer(0)]],
    device const float  *uold             [[buffer(1)]],
#ifdef MHD
    device const float  *bold             [[buffer(2)]],
    device atomic_uint  *data_buf         [[buffer(3)]],
    constant int        &head_idx         [[buffer(4)]],
    constant int        &num_octs         [[buffer(5)]],
    constant float      &dx               [[buffer(6)]],
    constant float      &gamma            [[buffer(7)]],
    constant float      &smallr           [[buffer(8)]],
    constant float      &smallc2          [[buffer(9)]],
    constant float      &courant_factor   [[buffer(10)]],
    constant int        &induction        [[buffer(11)]],
    constant float      *constant_gravity [[buffer(12)]],
    device const float  *f                [[buffer(13)]],
#else
    device atomic_uint  *data_buf         [[buffer(2)]],
    constant int        &head_idx         [[buffer(3)]],
    constant int        &num_octs         [[buffer(4)]],
    constant float      &dx               [[buffer(5)]],
    constant float      &gamma            [[buffer(6)]],
    constant float      &smallr           [[buffer(7)]],
    constant float      &smallc2          [[buffer(8)]],
    constant float      &courant_factor   [[buffer(9)]],
    constant float      *constant_gravity [[buffer(10)]],
    device const float  *f                [[buffer(11)]],
#endif
    uint tid    [[thread_position_in_threadgroup]],
    uint bid    [[threadgroup_position_in_grid]],
    uint lane   [[thread_index_in_simdgroup]],
    uint sg_idx [[simdgroup_index_in_threadgroup]],
    uint sg_size [[threads_per_simdgroup]])
{
    float mass_loc, ekin_loc, eint_loc, emag_loc, dt_loc;

    int oct_offset = int(bid * 1024 + tid) / TWOTONDIM;
    if (oct_offset >= num_octs) {
        dt_loc = HUGE_VALF; mass_loc = 0.0f; ekin_loc = 0.0f; eint_loc = 0.0f; emag_loc = 0.0f;
    } else {
        int oct_idx  = head_idx + oct_offset;          /* 1-based */
        int cell_idx = int(tid) % TWOTONDIM + 1;       /* 1-based (1..8) */
        if (grid[oct_idx - 1].refined[cell_idx - 1] != 0) {
            dt_loc = HUGE_VALF; mass_loc = 0.0f; ekin_loc = 0.0f; eint_loc = 0.0f; emag_loc = 0.0f;
        } else {
            conserved_t cv;
            cv.density    = u_get(uold, oct_idx, 1, cell_idx);
            cv.momentum_x = u_get(uold, oct_idx, 2, cell_idx);
            cv.momentum_y = u_get(uold, oct_idx, 3, cell_idx);
            cv.momentum_z = u_get(uold, oct_idx, 4, cell_idx);
            cv.energy     = u_get(uold, oct_idx, 5, cell_idx);
#ifdef MHD
            cv.Bx = 0.5f * (b_get(bold, oct_idx, 1, cell_idx) + b_get(bold, oct_idx, 4, cell_idx));
            cv.By = 0.5f * (b_get(bold, oct_idx, 2, cell_idx) + b_get(bold, oct_idx, 5, cell_idx));
            cv.Bz = 0.5f * (b_get(bold, oct_idx, 3, cell_idx) + b_get(bold, oct_idx, 6, cell_idx));
#endif
            float vol = dx*dx*dx;
            primitive_t pv = conserved_2_primitive(cv, gamma, smallr, smallc2);
#ifdef MHD
            float pressure = max(pv.pressure, pv.density * smallc2 / gamma);
            float b2 = magnitude_squared(pv.Bx, pv.By, pv.Bz);
            mass_loc = cv.density * vol;
#else
            mass_loc = pv.density * vol;
#endif
            ekin_loc = cv.energy  * vol;
#ifdef MHD
            eint_loc = pressure / (gamma - 1.0f) * vol;
            emag_loc = 0.5f * b2 * vol;
            float rinv = 1.0f / pv.density;
            float a2 = gamma * pressure * rinv;
            float ctot;
            if (induction != 0) {
                float smallc = sqrt(smallc2);
                ctot = abs(pv.velocity_x) + abs(pv.velocity_y) + abs(pv.velocity_z) + 3.0f * smallc;
            } else {
                float c2v = 0.5f * (b2 * rinv + a2);
                float cfx = sqrt(c2v + sqrt(max(c2v * c2v - a2 * pv.Bx * pv.Bx * rinv, 0.0f)));
                float cfy = sqrt(c2v + sqrt(max(c2v * c2v - a2 * pv.By * pv.By * rinv, 0.0f)));
                float cfz = sqrt(c2v + sqrt(max(c2v * c2v - a2 * pv.Bz * pv.Bz * rinv, 0.0f)));
                ctot = abs(pv.velocity_x) + cfx + abs(pv.velocity_y) + cfy + abs(pv.velocity_z) + cfz;
            }
#else
            eint_loc = pv.pressure / (gamma - 1.0f) * vol;
            emag_loc = 0.0f;
            float cs   = sqrt(gamma * pv.pressure / pv.density);
            float ctot = abs(pv.velocity_x) + abs(pv.velocity_y) + abs(pv.velocity_z) + 3.0f*cs;
#endif
            float grav;
#ifdef GRAV
            grav = abs(f[(oct_idx - 1)*3*8 + 0*8 + (cell_idx - 1)]) +
                   abs(f[(oct_idx - 1)*3*8 + 1*8 + (cell_idx - 1)]) +
                   abs(f[(oct_idx - 1)*3*8 + 2*8 + (cell_idx - 1)]);
#else
            grav = abs(constant_gravity[0]) + abs(constant_gravity[1]) + abs(constant_gravity[2]);
#endif
            grav = grav * dx / (ctot*ctot);
            grav = max(grav, 0.0001f);
            dt_loc = 2.0f*courant_factor*dx / (ctot * (1.0f + sqrt(1.0f + 2.0f*courant_factor*grav)));
        }
    }

    /* --- SIMD-group reduction (level 1): 32 groups × 32 lanes --- */
    threadgroup float tg_dt  [32];
    threadgroup float tg_mass[32];
    threadgroup float tg_ekin[32];
    threadgroup float tg_eint[32];
    threadgroup float tg_emag[32];

    float warp_dt   = simd_min(dt_loc);
    float warp_mass = simd_sum(mass_loc);
    float warp_ekin = simd_sum(ekin_loc);
    float warp_eint = simd_sum(eint_loc);
    float warp_emag = simd_sum(emag_loc);

    if (lane == 0) {
        tg_dt  [sg_idx] = warp_dt;
        tg_mass[sg_idx] = warp_mass;
        tg_ekin[sg_idx] = warp_ekin;
        tg_eint[sg_idx] = warp_eint;
        tg_emag[sg_idx] = warp_emag;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    /* --- Threadgroup reduction (level 2) — SIMD-group 0 only.
     * 32 SIMD groups fit exactly into the 32 lanes of group 0: no padding. --- */
    if (sg_idx == 0) {
        float blk_dt   = simd_min(tg_dt  [lane]);
        float blk_mass = simd_sum(tg_mass[lane]);
        float blk_ekin = simd_sum(tg_ekin[lane]);
        float blk_eint = simd_sum(tg_eint[lane]);
        float blk_emag = simd_sum(tg_emag[lane]);

        if (lane == 0) {
            atomic_add_float     (&data_buf[0], blk_mass);
            atomic_add_float     (&data_buf[1], blk_ekin);
            atomic_add_float     (&data_buf[2], blk_eint);
            atomic_add_float     (&data_buf[3], blk_emag);
            atomic_min_float_bits(&data_buf[4], blk_dt);
        }
    }
}

#ifndef MHD
/* ===========================================================================
 * hydro_integrator_kernel — mirrors attributes(global) hydro_integrator_kernel
 *                           in gpu_hydro.cuf (no MHD, no scalars, no GRAV).
 *
 * One threadgroup per oct (nsubgrid=1 → 1 subgrid = 1 oct).
 * 64 threads/threadgroup.
 * For levelmin==levelmax: zero_fine_fluxes and coarse_cell_update are skipped.
 * ========================================================================= */
kernel void hydro_integrator_kernel(
    device const oct_t  *grid           [[buffer(0)]],
    device const float  *uold           [[buffer(1)]],
    device float        *unew           [[buffer(2)]],
    device const int    *nbor           [[buffer(3)]],  /* int[ngridmax][SUBGRIDSIZE] */
    constant int        &head_idx       [[buffer(4)]],
    constant int        &num_subgrids   [[buffer(5)]],
    constant int        &ngridmax       [[buffer(6)]],
    constant int        &ilevel         [[buffer(7)]],
    constant int        &levelmin       [[buffer(8)]],
    constant int        &levelmax       [[buffer(9)]],
    constant float      &gamma          [[buffer(10)]],
    constant float      &smallr         [[buffer(11)]],
    constant float      &smallc2        [[buffer(12)]],
    constant float      &dt             [[buffer(13)]],
    constant float      &dx             [[buffer(14)]],
    constant int        &slope          [[buffer(15)]],
    constant int        &riemann        [[buffer(16)]],
    constant float      *constant_gravity [[buffer(17)]],
    device const int    *father           [[buffer(18)]],
    device const float  *f                [[buffer(19)]],
    uint block_idx      [[threadgroup_position_in_grid]],
    uint thread_idx     [[thread_position_in_threadgroup]],
    uint threads_per_tg [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_subgrids) return;

    float dtdx = dt / dx;

    /* Threadgroup-local data — mirrors the `shared` variables in gpu_hydro.cuf */
    threadgroup local_subgrid_t  local_subgrid;
    threadgroup interfaces_x_t   left_x, right_x;
    threadgroup interfaces_y_t   left_y, right_y;
    threadgroup interfaces_z_t   left_z, right_z;

    /* =========================================================================
     * Convert from conserved to primitive variables
     * ========================================================================= */
    subgrid_conserved_2_primitive(
        grid, uold, nbor, constant_gravity, f,
        head_idx, int(block_idx), int(thread_idx),
        gamma, smallr, smallc2, dt, threads_per_tg, local_subgrid);

    /* =========================================================================
     * MUSCL tracing: compute interface states
     * ========================================================================= */
    trace_3d(local_subgrid, int(thread_idx), threads_per_tg,
             left_x, right_x, left_y, right_y, left_z, right_z,
             gamma, smallr, smallc2, dtdx, slope);

    /* =========================================================================
     * Compute fluxes via Riemann solver (overwrites left_interfaces_*)
     * ========================================================================= */
    riemann_driver(left_x, right_x, left_y, right_y, left_z, right_z,
                   int(thread_idx), threads_per_tg, gamma, smallr, smallc2, riemann);

    /* Zero fine-level fluxes at refined faces so they do not corrupt coarse update. */
    if (ilevel < levelmax) {
        zero_fine_fluxes(local_subgrid, int(thread_idx), threads_per_tg,
                         left_x, left_y, left_z);
    }

    /* =========================================================================
     * Update conserved variables at current level
     * ========================================================================= */
    conservative_update(unew, nbor,
                        left_x, left_y, left_z,
                        head_idx, int(block_idx), int(thread_idx), threads_per_tg,
                        dtdx);

    /* Correct the coarser level's unew via boundary flux contributions. */
    if (ilevel > levelmin) {
        float cfs = dtdx / float(TWOTONDIM);   /* dtdx / 8 */
        coarse_cell_update(unew, nbor, grid, father,
                           left_x, left_y, left_z,
                           head_idx, int(block_idx), ngridmax,
                           int(thread_idx), cfs);
    }
}
#endif

/* ===========================================================================
 * upload_kernel — restriction (averaging down) for fine→coarse levels.
 * Mirrors attributes(global) upload_kernel in gpu_hydro.cuf (no MHD, no NENER).
 *
 * One thread per fine oct (ilevel+1).  Looks up the parent oct via father[],
 * then averages all 8 children uold values into the parent cell.
 * For internal_energy!=0: converts total→internal energy before averaging,
 * then converts back to total using the averaged parent momenta.
 * ========================================================================= */
kernel void upload_kernel(
    device const oct_t  *grid             [[buffer(0)]],
    device const int    *father           [[buffer(1)]],
    device float        *uold             [[buffer(2)]],
#ifdef MHD
    device float        *bold             [[buffer(3)]],
    constant int        &head_idx         [[buffer(4)]],
    constant int        &num_octs         [[buffer(5)]],
    constant int        &internal_energy  [[buffer(6)]],
    constant float      &gamma            [[buffer(7)]],
    constant float      &smallr           [[buffer(8)]],
    constant float      &smallc2          [[buffer(9)]],
#else
    constant int        &head_idx         [[buffer(3)]],
    constant int        &num_octs         [[buffer(4)]],
    constant int        &internal_energy  [[buffer(5)]],
    constant float      &gamma            [[buffer(6)]],
    constant float      &smallr           [[buffer(7)]],
    constant float      &smallc2          [[buffer(8)]],
#endif
    uint tid [[thread_position_in_grid]])
{
    if ((int)tid >= num_octs) return;
    int oct_idx    = head_idx + (int)tid;   /* 1-based fine oct */
    int father_idx = father[oct_idx - 1];   /* 1-based parent oct */

    /* Child cell position within parent: ckey difference (0 or 1 per axis) */
    int ic = grid[oct_idx - 1].ckey[0] - 2 * grid[father_idx - 1].ckey[0];
    int jc = grid[oct_idx - 1].ckey[1] - 2 * grid[father_idx - 1].ckey[1];
    int kc = grid[oct_idx - 1].ckey[2] - 2 * grid[father_idx - 1].ckey[2];
    int cell_idx = 1 + ic + 2*jc + 4*kc;   /* 1-based parent cell (1..8) */

    float inv8 = 1.0f / 8.0f;

    /* Average all NVAR conserved variables from 8 fine children → parent cell */
    for (int ivar = 1; ivar <= (NVAR); ivar++) {
        float avg = 0.0f;
        for (int ind = 1; ind <= 8; ind++)
            avg += u_get(uold, oct_idx, ivar, ind);
        u_set(uold, father_idx, ivar, cell_idx, avg * inv8);
    }
#ifdef MHD
    for (int idim = 0; idim < 3; idim++) {
        float low = 0.0f;
        float high = 0.0f;
        for (int ind = 0; ind < 8; ind++) {
            if (((ind >> idim) & 1) == 0)
                low += bold[(oct_idx - 1) * 48 + idim * 8 + ind];
            else
                high += bold[(oct_idx - 1) * 48 + (idim + 3) * 8 + ind];
        }
        bold[(father_idx - 1) * 48 + idim * 8 + cell_idx - 1] = 0.25f * low;
        bold[(father_idx - 1) * 48 + (idim + 3) * 8 + cell_idx - 1] = 0.25f * high;
    }
#endif

    /* Non-conservative upload: average internal energy, restore total */
    if (internal_energy != 0) {
        float smalle = smallc2 / gamma / (gamma - 1.0f);
        float eint_sum = 0.0f;
        for (int ind = 1; ind <= 8; ind++) {
            float dens = max(u_get(uold, oct_idx, 1, ind), smallr);
            float mx   = u_get(uold, oct_idx, 2, ind);
            float my   = u_get(uold, oct_idx, 3, ind);
            float mz   = u_get(uold, oct_idx, 4, ind);
            float etot = u_get(uold, oct_idx, 5, ind);
            float ekin = 0.5f * (mx*mx + my*my + mz*mz) / dens;
#ifdef MHD
            float eb = 0.0f;
            for (int idim = 0; idim < 3; idim++) {
                float b = 0.5f * (bold[(oct_idx - 1) * 48 + idim * 8 + ind - 1] + bold[(oct_idx - 1) * 48 + (idim + 3) * 8 + ind - 1]);
                eb += 0.5f * b * b;
            }
            eint_sum += max(etot - ekin - eb, smalle * dens);
#else
            eint_sum += max(etot - ekin, smalle * dens);
#endif
        }
        /* Recompute parent kinetic energy from averaged parent momenta */
        float dens = max(u_get(uold, father_idx, 1, cell_idx), smallr);
        float mx   = u_get(uold, father_idx, 2, cell_idx);
        float my   = u_get(uold, father_idx, 3, cell_idx);
        float mz   = u_get(uold, father_idx, 4, cell_idx);
        float ekin = 0.5f * (mx*mx + my*my + mz*mz) / dens;
#ifdef MHD
        float eb = 0.0f;
        for (int idim = 0; idim < 3; idim++) {
            float b = 0.5f * (bold[(father_idx - 1) * 48 + idim * 8 + cell_idx - 1] + bold[(father_idx - 1) * 48 + (idim + 3) * 8 + cell_idx - 1]);
            eb += 0.5f * b * b;
        }
        u_set(uold, father_idx, 5, cell_idx, eint_sum * inv8 + ekin + eb);
#else
        u_set(uold, father_idx, 5, cell_idx, eint_sum * inv8 + ekin);
#endif
    }
}

kernel void sync_hydro_kernel(
    device float       *uold             [[buffer(0)]],
    device const float *f                [[buffer(1)]],
    constant float     *constant_gravity [[buffer(2)]],
    constant int       &head_idx         [[buffer(3)]],
    constant int       &num_octs         [[buffer(4)]],
    constant float     &gamma            [[buffer(5)]],
    constant float     &smallr           [[buffer(6)]],
    constant float     &smallc2          [[buffer(7)]],
    constant float     &dt               [[buffer(8)]],
    uint2               thread_idx       [[thread_position_in_threadgroup]],
    uint2               group_idx        [[threadgroup_position_in_grid]],
    uint2               group_dim        [[threads_per_threadgroup]])
{
    int oct_idx = group_idx.x * group_dim.y + thread_idx.y;
    if (oct_idx >= num_octs) return;
    oct_idx = head_idx + oct_idx;  /* 1-based */
    int cell_idx = thread_idx.x + 1; /* 1-based */

    float rho = max(u_get(uold, oct_idx, 1, cell_idx), smallr);
    float mx  = u_get(uold, oct_idx, 2, cell_idx);
    float my  = u_get(uold, oct_idx, 3, cell_idx);
    float mz  = u_get(uold, oct_idx, 4, cell_idx);
    float ener = u_get(uold, oct_idx, 5, cell_idx) - 0.5f * (mx*mx + my*my + mz*mz) / rho;

#ifdef GRAV
    mx += rho * f[(oct_idx - 1)*3*8 + 0*8 + (cell_idx - 1)] * dt;
    my += rho * f[(oct_idx - 1)*3*8 + 1*8 + (cell_idx - 1)] * dt;
    mz += rho * f[(oct_idx - 1)*3*8 + 2*8 + (cell_idx - 1)] * dt;
#else
    mx += rho * constant_gravity[0] * dt;
    my += rho * constant_gravity[1] * dt;
    mz += rho * constant_gravity[2] * dt;
#endif

    ener += 0.5f * (mx*mx + my*my + mz*mz) / rho;
    u_set(uold, oct_idx, 2, cell_idx, mx);
    u_set(uold, oct_idx, 3, cell_idx, my);
    u_set(uold, oct_idx, 4, cell_idx, mz);
    u_set(uold, oct_idx, 5, cell_idx, ener);
}

kernel void grav_hydro_kernel(
    device const float *uold             [[buffer(0)]],
    device float       *unew             [[buffer(1)]],
    device const float *f                [[buffer(2)]],
    constant float     *constant_gravity [[buffer(3)]],
    constant int       &head_idx         [[buffer(4)]],
    constant int       &num_octs         [[buffer(5)]],
    constant float     &gamma            [[buffer(6)]],
    constant float     &smallr           [[buffer(7)]],
    constant float     &smallc2          [[buffer(8)]],
    constant float     &dt               [[buffer(9)]],
    uint2               thread_idx       [[thread_position_in_threadgroup]],
    uint2               group_idx        [[threadgroup_position_in_grid]],
    uint2               group_dim        [[threads_per_threadgroup]])
{
    int oct_idx = group_idx.x * group_dim.y + thread_idx.y;
    if (oct_idx >= num_octs) return;
    oct_idx = head_idx + oct_idx;  /* 1-based */
    int cell_idx = thread_idx.x + 1; /* 1-based */

    float d = max(u_get(unew, oct_idx, 1, cell_idx), smallr);
    float mx = u_get(unew, oct_idx, 2, cell_idx);
    float my = u_get(unew, oct_idx, 3, cell_idx);
    float mz = u_get(unew, oct_idx, 4, cell_idx);
    float e_prim = u_get(unew, oct_idx, 5, cell_idx) - 0.5f * (mx*mx + my*my + mz*mz) / d;
    float d_old = max(u_get(uold, oct_idx, 1, cell_idx), smallr);
    float fact = d_old / d * dt;

#ifdef GRAV
    mx += d * f[(oct_idx - 1)*3*8 + 0*8 + (cell_idx - 1)] * fact;
    my += d * f[(oct_idx - 1)*3*8 + 1*8 + (cell_idx - 1)] * fact;
    mz += d * f[(oct_idx - 1)*3*8 + 2*8 + (cell_idx - 1)] * fact;
#else
    mx += d * constant_gravity[0] * fact;
    my += d * constant_gravity[1] * fact;
    mz += d * constant_gravity[2] * fact;
#endif

    u_set(unew, oct_idx, 2, cell_idx, mx);
    u_set(unew, oct_idx, 3, cell_idx, my);
    u_set(unew, oct_idx, 4, cell_idx, mz);
    u_set(unew, oct_idx, 5, cell_idx, e_prim + 0.5f * (mx*mx + my*my + mz*mz) / d);
}

#ifdef MHD
#include "mhd_uct.metal"
#endif

/*
 * metal/kernels/refine.metal
 *
 * Metal port of gpu/gpu_refine.cuf — AMR neighbour-array construction.
 * Scope: NDIM=3; AMR supports NSUBGRID=1 or 2.
 *
 * Hash table layout (set by insert_hash_kernel in hash.metal):
 *   hash_key : device long[hash_size]   (64-bit Hilbert key; 0 == empty)
 *   hash_val : device int[hash_size]    (1-based oct index; 0 == not found)
 * Both are read non-atomically here because hash.metal's insert dispatch
 * has already completed (waitUntilCompleted in metal_bridge.mm).
 *
 * Kernel:
 *   build_nbor_kernel computes all neighbours for one subgrid.
 */

#include <metal_stdlib>
#include "../metal_types.h"
#include "../metal_config.h"
using namespace metal;

constant int NSUBGRID_RF       = NSUBGRID;
constant int NSUBGRIDP2_RF     = NSUBGRID_RF + 2;
constant int NSUBGRIDTONDIM_RF = NSUBGRID_RF * NSUBGRID_RF * NSUBGRID_RF;
constant int SUBGRIDSIZE_RF    = NSUBGRIDP2_RF * NSUBGRIDP2_RF * NSUBGRIDP2_RF;

/* =========================================================================
 * fnv64 / hash_bucket / hash_get — read-only hash helpers.
 * These mirror the insert-side versions in hash.metal but operate on plain
 * (non-atomic) device pointers, which is safe because the insert dispatch
 * completed before this kernel is enqueued.
 * ========================================================================= */
inline ulong fnv64_r(long key_signed)
{
    ulong key = as_type<ulong>(key_signed);
    ulong h   = 14695981039346656037UL;
    for (int j = 0; j < 8; j++) {
        ulong b = (key >> (8 * j)) & 0xFFUL;
        h ^= b;
        h *= 1099511628211UL;
    }
    return h;
}

inline int hash_bucket_r(long key, int hash_size)
{
    return int(fnv64_r(key) % ulong(hash_size)) + 1;
}

inline int hash_get(device const long* hash_key, device const int* hash_val,
                    int hash_size, long key)
{
    int b = hash_bucket_r(key, hash_size);
    for (;;) {
        long cur = hash_key[b - 1];
        if (cur == key) return hash_val[b - 1];
        if (cur == 0L)  return 0;
        b = (b % hash_size) + 1;
    }
}

/* =========================================================================
 * build_nbor_kernel — replaces the multi-launch update_nbor_array loop from
 *                     r_set_grid_device in gpu_manager.cuf.
 *
 * Each thread handles one subgrid and computes all neighbour indices.  The
 * inner loop body mirrors update_nbor_array:
 *
 *   ind = 1 + i_off + j_off*nsubgridp2 + k_off*nsubgridp2^2
 *   ckey_n[d] = (ckey[d]/nsubgrid)*nsubgrid + off - 1
 *   if periodic: wrap when ckey_n < box_min or ckey_n >= box_max
 *   key = key_off + ix + iy*nx + iz*nx*nx
 *   nbor[(subgrid_idx-1)*subgridsize + (ind-1)] = hash_get(key)
 *
 * Thread layout: 128 threads/threadgroup.
 * ========================================================================= */
kernel void build_nbor_kernel(
    device const oct_t   *grid         [[buffer(0)]],
    device int           *nbor         [[buffer(1)]],
    device const long    *hash_key     [[buffer(2)]],   /* non-atomic long */
    device const int     *hash_val     [[buffer(3)]],   /* non-atomic int  */
    constant int         &hash_size    [[buffer(4)]],
    constant int         &ckey_max_l   [[buffer(5)]],
    constant long        &key_off_l    [[buffer(6)]],
    constant int         *box_ckey_min [[buffer(7)]],   /* [3] */
    constant int         *box_ckey_max [[buffer(8)]],   /* [3] */
    constant int         *periodic     [[buffer(9)]],   /* [3]: non-zero = true */
    constant int         &head_idx     [[buffer(10)]],
    constant int         &num_subgrids [[buffer(11)]],
    uint thread_idx [[thread_position_in_grid]])
{
    if (int(thread_idx) >= num_subgrids) return;
    int subgrid_idx = head_idx + int(thread_idx);   /* 1-based */
    int oct_idx = (subgrid_idx - 1) * NSUBGRIDTONDIM_RF + 1;

    long nx = long(ckey_max_l);

    for (int k_off = 0; k_off < NSUBGRIDP2_RF; k_off++) {
        for (int j_off = 0; j_off < NSUBGRIDP2_RF; j_off++) {
            for (int i_off = 0; i_off < NSUBGRIDP2_RF; i_off++) {
                int ind = 1 + i_off
                            + NSUBGRIDP2_RF * j_off
                            + NSUBGRIDP2_RF * NSUBGRIDP2_RF * k_off;

                int ckey_n[3];
                ckey_n[0] = (grid[oct_idx - 1].ckey[0] / NSUBGRID_RF) * NSUBGRID_RF + i_off - 1;
                ckey_n[1] = (grid[oct_idx - 1].ckey[1] / NSUBGRID_RF) * NSUBGRID_RF + j_off - 1;
                ckey_n[2] = (grid[oct_idx - 1].ckey[2] / NSUBGRID_RF) * NSUBGRID_RF + k_off - 1;

                for (int d = 0; d < 3; d++) {
                    if (periodic[d]) {
                        if (ckey_n[d] <  box_ckey_min[d]) ckey_n[d] = box_ckey_max[d] - 1;
                        if (ckey_n[d] >= box_ckey_max[d]) ckey_n[d] = box_ckey_min[d];
                    }
                }

                long ix  = long(ckey_n[0]);
                long iy  = long(ckey_n[1]);
                long iz  = long(ckey_n[2]);
                long key = key_off_l + ix + iy * nx + iz * nx * nx;

                int nbor_idx = hash_get(hash_key, hash_val, hash_size, key);
                nbor[(subgrid_idx - 1) * SUBGRIDSIZE_RF + (ind - 1)] = nbor_idx;
            }
        }
    }
}

/* =========================================================================
 * update_father_kernel — populate father[oct_abs_0] with the 1-based parent
 * oct index by hashing the parent's Cartesian key.
 * Mirrors update_father_array<<<128>>> in gpu_refine.cuf.
 * Per-level Hilbert parameters are passed directly (not via s_ckey_max_dev)
 * so this kernel can be called once per level.
 * ========================================================================= */
kernel void update_father_kernel(
    device const oct_t  *grid     [[buffer(0)]],
    device int          *father   [[buffer(1)]],
    device const long   *hash_key [[buffer(2)]],
    device const int    *hash_val [[buffer(3)]],
    constant int        &hash_size  [[buffer(4)]],
    constant int        &ckey_max_l [[buffer(5)]],
    constant long       &key_off_l  [[buffer(6)]],
    constant int        &head_idx   [[buffer(7)]],
    constant int        &num_octs   [[buffer(8)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);

    /* Parent ckey: divide current ckey by 2 (right-shift) */
    long ix = long(grid[oct_abs_0].ckey[0] / 2);
    long iy = long(grid[oct_abs_0].ckey[1] / 2);
    long iz = long(grid[oct_abs_0].ckey[2] / 2);
    long nx = long(ckey_max_l);

    long parent_key = key_off_l + ix + iy * nx + iz * nx * nx;
    father[oct_abs_0] = hash_get(hash_key, hash_val, hash_size, parent_key);
}

/* =========================================================================
 * AMR refinement kernels — port of gpu/gpu_refine.cuf
 * Scope: HYDRO=1, NDIM=3, NSUBGRID=1 or 2, no CUB.
 *
 * Memory layout (0-based):
 *   uold    : float[8*(NVAR)*ntotal]  cell_0 + 8*ivar_0 + 8*(NVAR)*oct_abs_0
 *   flag1/2 : int[8*ntotal]           cell_0 + 8*oct_abs_0
 *   nbor    : int[subgridsize*ntotal] (subgrid_abs_0)*subgridsize + ind_0
 *   grid    : oct_t[ntotal]           0-based oct index
 *   father  : int[ntotal]             1-based parent oct index per oct
 *   swap_*  : int[ntotal]             1-based oct index
 *   prefix  : int[ntotal]
 *
 * Level params uploaded by mtl_upload_level_params (0-based in device buffers):
 *   ckey_max_dev[lev-1]            = ckey_max(lev)
 *   key_off_dev[lev-1]             = key_off(lev)   (long)
 *   box_ckey_min_dev[3*(lev-1)+d]  = box_ckey_min(d+1, lev)
 *   box_ckey_max_dev[3*(lev-1)+d]  = box_ckey_max(d+1, lev)
 *   periodic_dev[d]                = periodic(d+1)
 * ========================================================================= */

#ifndef NVAR
#define NVAR 5
#endif

inline float rf_u_get(device const float *u, int oct, int ivar, int cell)
{
    return u[(oct - 1) * NVAR * 8 + ivar * 8 + cell];
}

inline void rf_u_set(device float *u, int oct, int ivar, int cell, float value)
{
    u[(oct - 1) * NVAR * 8 + ivar * 8 + cell] = value;
}

inline float rf_b_get(device const float *b, int oct, int ivar, int cell)
{
    return b[(oct - 1) * 6 * 8 + ivar * 8 + cell];
}

inline void rf_b_set(device float *b, int oct, int ivar, int cell, float value)
{
    b[(oct - 1) * 6 * 8 + ivar * 8 + cell] = value;
}

inline void rf_limiter(thread const float *a, thread float *w, int interpol_type)
{
    if (interpol_type == 1) {
        for (int d = 0; d < 3; d++) {
            float dl = 0.5f * (a[2*d + 2] - a[0]);
            float dr = 0.5f * (a[0] - a[2*d + 1]);
            w[d] = dl * dr <= 0.0f ? 0.0f : copysign(min(abs(dl), abs(dr)), dl);
        }
        return;
    }

    for (int d = 0; d < 3; d++) w[d] = 0.25f * (a[2*d + 2] - a[2*d + 1]);
    if (interpol_type == 3) return;

    float corner_max = -INFINITY;
    float corner_min = INFINITY;
    for (int c = 0; c < 8; c++) {
        float value = a[0];
        for (int d = 0; d < 3; d++) value += 2.0f * w[d] * (float((c >> d) & 1) - 0.5f);
        corner_max = max(corner_max, value);
        corner_min = min(corner_min, value);
    }
    float kernel_max = a[1];
    float kernel_min = a[1];
    for (int n = 2; n <= 6; n++) {
        kernel_max = max(kernel_max, a[n]);
        kernel_min = min(kernel_min, a[n]);
    }
    float max_limiter = 0.0f;
    float min_limiter = 0.0f;
    float dk = a[0] - kernel_max;
    float dc = a[0] - corner_max;
    if (dk * dc > 0.0f) max_limiter = min(1.0f, dk / dc);
    dk = a[0] - kernel_min;
    dc = a[0] - corner_min;
    if (dk * dc > 0.0f) min_limiter = min(1.0f, dk / dc);
    float limiter = min(min_limiter, max_limiter);
    for (int d = 0; d < 3; d++) w[d] *= limiter;
}

inline void rf_tvd2(thread const float *b, thread float *s, int interpol_type)
{
    if (interpol_type == 3) {
        s[0] = 0.5f * (b[0] - b[1]) + 0.5f * (b[2] - b[0]);
        s[1] = 0.5f * (b[0] - b[3]) + 0.5f * (b[4] - b[0]);
        return;
    }
    float r = float(interpol_type);
    float dl = r * (b[0] - b[1]);
    float dr = r * (b[2] - b[0]);
    float dc = 0.5f * (dl + dr) / r;
    float limit = dl * dr <= 0.0f ? 0.0f : min(min(abs(dl), abs(dr)), abs(dc));
    s[0] = copysign(limit, dc);
    dl = r * (b[0] - b[3]);
    dr = r * (b[4] - b[0]);
    dc = 0.5f * (dl + dr) / r;
    limit = dl * dr <= 0.0f ? 0.0f : min(min(abs(dl), abs(dr)), abs(dc));
    s[1] = copysign(limit, dc);
}

inline int rf_uface(int x, int y, int z) { return ((x + 1) * 2 + y) * 2 + z; }
inline int rf_vface(int x, int y, int z) { return (x * 3 + y + 1) * 2 + z; }
inline int rf_wface(int x, int y, int z) { return (x * 2 + y) * 3 + z + 1; }

inline void rf_interpol_mhd(thread float *u1, thread float *u2,
                            thread float *b1, thread float *b2,
                            thread float *b3, thread bool *refined,
                            int interpol_var, int interpol_type, float smallr)
{
    if (interpol_var == 1) {
        for (int n = 0; n <= 6; n++) {
            float rho = max(u1[n*NVAR], smallr);
            float ekin = 0.0f;
            float emag = 0.0f;
            for (int d = 0; d < 3; d++) {
                float mom = u1[n*NVAR + d + 1];
                float bc = 0.5f * (b1[n*6 + d] + b1[n*6 + d + 3]);
                ekin += 0.5f * mom * mom / rho;
                emag += 0.5f * bc * bc;
            }
            u1[n*NVAR + 4] -= ekin + emag;
        }
    }

    for (int ivar = 0; ivar < NVAR; ivar++) {
        float a[7];
        float w[3] = {0.0f, 0.0f, 0.0f};
        for (int n = 0; n <= 6; n++) a[n] = u1[n*NVAR + ivar];
        if (interpol_type > 0) rf_limiter(a, w, interpol_type);
        for (int c = 0; c < 8; c++) {
            float value = a[0];
            for (int d = 0; d < 3; d++) value += w[d] * (float((c >> d) & 1) - 0.5f);
            u2[c*NVAR + ivar] = value;
        }
    }

    float uf[12];
    float vf[12];
    float wf[12];
    float b[5];
    float s[2];

    for (int side = 0; side < 2; side++) {
        int comp = side == 0 ? 0 : 3;
        b[0] = b1[comp];
        b[1] = b1[3*6 + comp]; b[2] = b1[4*6 + comp];
        b[3] = b1[5*6 + comp]; b[4] = b1[6*6 + comp];
        s[0] = s[1] = 0.0f;
        if (interpol_type > 0) rf_tvd2(b, s, interpol_type);
        int x = side == 0 ? -1 : 1;
        for (int j = 0; j < 2; j++) for (int k = 0; k < 2; k++)
            uf[rf_uface(x,j,k)] = b[0] + 0.5f*s[0]*(float(j)-0.5f) + 0.5f*s[1]*(float(k)-0.5f);
    }
    for (int side = 0; side < 2; side++) {
        int comp = side == 0 ? 1 : 4;
        b[0] = b1[comp];
        b[1] = b1[1*6 + comp]; b[2] = b1[2*6 + comp];
        b[3] = b1[5*6 + comp]; b[4] = b1[6*6 + comp];
        s[0] = s[1] = 0.0f;
        if (interpol_type > 0) rf_tvd2(b, s, interpol_type);
        int y = side == 0 ? -1 : 1;
        for (int i = 0; i < 2; i++) for (int k = 0; k < 2; k++)
            vf[rf_vface(i,y,k)] = b[0] + 0.5f*s[0]*(float(i)-0.5f) + 0.5f*s[1]*(float(k)-0.5f);
    }
    for (int side = 0; side < 2; side++) {
        int comp = side == 0 ? 2 : 5;
        b[0] = b1[comp];
        b[1] = b1[1*6 + comp]; b[2] = b1[2*6 + comp];
        b[3] = b1[3*6 + comp]; b[4] = b1[4*6 + comp];
        s[0] = s[1] = 0.0f;
        if (interpol_type > 0) rf_tvd2(b, s, interpol_type);
        int z = side == 0 ? -1 : 1;
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
            wf[rf_wface(i,j,z)] = b[0] + 0.5f*s[0]*(float(i)-0.5f) + 0.5f*s[1]*(float(j)-0.5f);
    }

    for (int j = 0; j < 2; j++) for (int k = 0; k < 2; k++) {
        if (refined[0]) uf[rf_uface(-1,j,k)] = b3[(0*8 + 1 + 2*j + 4*k)*6 + 3];
        if (refined[1]) uf[rf_uface( 1,j,k)] = b3[(1*8 +     2*j + 4*k)*6    ];
    }
    for (int i = 0; i < 2; i++) for (int k = 0; k < 2; k++) {
        if (refined[2]) vf[rf_vface(i,-1,k)] = b3[(2*8 + i + 2 + 4*k)*6 + 4];
        if (refined[3]) vf[rf_vface(i, 1,k)] = b3[(3*8 + i     + 4*k)*6 + 1];
    }
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
        if (refined[4]) wf[rf_wface(i,j,-1)] = b3[(4*8 + i + 2*j + 4)*6 + 5];
        if (refined[5]) wf[rf_wface(i,j, 1)] = b3[(5*8 + i + 2*j    )*6 + 2];
    }

    float uxx = 0.0f, vyy = 0.0f, wzz = 0.0f;
    float uxyz = 0.0f, vxyz = 0.0f, wxyz = 0.0f;
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) for (int k = 0; k < 2; k++) {
        float ii = float(2*i - 1), jj = float(2*j - 1), kk = float(2*k - 1);
        uxx += 0.125f * (ii*jj*vf[rf_vface(i,2*j-1,k)] + ii*kk*wf[rf_wface(i,j,2*k-1)]);
        vyy += 0.125f * (jj*kk*wf[rf_wface(i,j,2*k-1)] + ii*jj*uf[rf_uface(2*i-1,j,k)]);
        wzz += 0.125f * (ii*kk*uf[rf_uface(2*i-1,j,k)] + jj*kk*vf[rf_vface(i,2*j-1,k)]);
        uxyz += 0.125f * ii*jj*kk*uf[rf_uface(2*i-1,j,k)];
        vxyz += 0.125f * ii*jj*kk*vf[rf_vface(i,2*j-1,k)];
        wxyz += 0.125f * ii*jj*kk*wf[rf_wface(i,j,2*k-1)];
    }
    for (int j = 0; j < 2; j++) for (int k = 0; k < 2; k++)
        uf[rf_uface(0,j,k)] = 0.5f*(uf[rf_uface(-1,j,k)] + uf[rf_uface(1,j,k)]) + uxx
                              + (float(k)-0.5f)*vxyz + (float(j)-0.5f)*wxyz;
    for (int i = 0; i < 2; i++) for (int k = 0; k < 2; k++)
        vf[rf_vface(i,0,k)] = 0.5f*(vf[rf_vface(i,-1,k)] + vf[rf_vface(i,1,k)]) + vyy
                              + (float(i)-0.5f)*wxyz + (float(k)-0.5f)*uxyz;
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
        wf[rf_wface(i,j,0)] = 0.5f*(wf[rf_wface(i,j,-1)] + wf[rf_wface(i,j,1)]) + wzz
                              + (float(j)-0.5f)*uxyz + (float(i)-0.5f)*vxyz;

    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) for (int k = 0; k < 2; k++) {
        int c = i + 2*j + 4*k;
        b2[c*6    ] = uf[rf_uface(i-1,j,k)]; b2[c*6 + 3] = uf[rf_uface(i,j,k)];
        b2[c*6 + 1] = vf[rf_vface(i,j-1,k)]; b2[c*6 + 4] = vf[rf_vface(i,j,k)];
        b2[c*6 + 2] = wf[rf_wface(i,j,k-1)]; b2[c*6 + 5] = wf[rf_wface(i,j,k)];
    }

    if (interpol_var == 1) {
        for (int c = 0; c < 8; c++) {
            float rho = max(u2[c*NVAR], smallr);
            float ekin = 0.0f;
            float emag = 0.0f;
            for (int d = 0; d < 3; d++) {
                float mom = u2[c*NVAR + d + 1];
                float bc = 0.5f * (b2[c*6 + d] + b2[c*6 + d + 3]);
                ekin += 0.5f * mom * mom / rho;
                emag += 0.5f * bc * bc;
            }
            u2[c*NVAR + 4] += ekin + emag;
        }
    }
}

inline void rf_prolong_mhd(device const oct_t *grid, device const float *uold,
                           device const float *bold, device const long *hash_key,
                           device const int *hash_val, int hash_size,
                           device const int *ckey_max_dev, device const long *key_off_dev,
                           device const int *box_ckey_min, device const int *box_ckey_max,
                           device const int *periodic, int ilevel, thread const int *chkey,
                           int central_oct, int central_cell, int interpol_var,
                           int interpol_type, float smallr, thread float *u2, thread float *b2)
{
    float u1[7*NVAR];
    float b1[7*6];
    float b3[6*8*6];
    bool refined[6];
    for (int ivar = 0; ivar < NVAR; ivar++) u1[ivar] = rf_u_get(uold, central_oct, ivar, central_cell);
    for (int ib = 0; ib < 6; ib++) b1[ib] = rf_b_get(bold, central_oct, ib, central_cell);

    for (int n = 0; n < 6; n++) {
        refined[n] = false;
        int nkey[3] = {chkey[0], chkey[1], chkey[2]};
        int d = n / 2;
        nkey[d] += (n & 1) == 0 ? -1 : 1;
        for (int q = 0; q < 3; q++) if (periodic[q]) {
            int bmin = box_ckey_min[3*(ilevel-1) + q];
            int bmax = box_ckey_max[3*(ilevel-1) + q];
            if (nkey[q] < bmin) nkey[q] = bmax - 1;
            if (nkey[q] >= bmax) nkey[q] = bmin;
        }
        int fkey[3] = {nkey[0]/2, nkey[1]/2, nkey[2]/2};
        int cell = (nkey[0] - 2*fkey[0]) + 2*(nkey[1] - 2*fkey[1]) + 4*(nkey[2] - 2*fkey[2]);
        long nx = long(ckey_max_dev[ilevel-2]);
        long key = key_off_dev[ilevel-2] + long(fkey[0]) + long(fkey[1])*nx + long(fkey[2])*nx*nx;
        int parent = hash_get(hash_key, hash_val, hash_size, key);
        if (parent > 0) {
            for (int ivar = 0; ivar < NVAR; ivar++) u1[(n+1)*NVAR + ivar] = rf_u_get(uold, parent, ivar, cell);
            for (int ib = 0; ib < 6; ib++) b1[(n+1)*6 + ib] = rf_b_get(bold, parent, ib, cell);
            if (grid[parent-1].refined[cell]) {
                nx = long(ckey_max_dev[ilevel-1]);
                key = key_off_dev[ilevel-1] + long(nkey[0]) + long(nkey[1])*nx + long(nkey[2])*nx*nx;
                int child = hash_get(hash_key, hash_val, hash_size, key);
                if (child > 0) {
                    refined[n] = true;
                    for (int c = 0; c < 8; c++) for (int ib = 0; ib < 6; ib++)
                        b3[(n*8 + c)*6 + ib] = rf_b_get(bold, child, ib, c);
                }
            }
        } else {
            for (int ivar = 0; ivar < NVAR; ivar++) u1[(n+1)*NVAR + ivar] = u1[ivar];
            for (int ib = 0; ib < 6; ib++) b1[(n+1)*6 + ib] = b1[ib];
        }
    }
    rf_interpol_mhd(u1, u2, b1, b2, b3, refined, interpol_var, interpol_type, smallr);
}

/* ---- Hilbert state tables (NDIM=3), from gpu/gpu_hilbert.cuf ---------- */
/* left_shift=4, right_shift=-1 → ISHFT(hkey,4) then ISHFT(hkey,-1)       */
/* Net effect per loop iteration: hkey = (hkey << 4) >> 1  =  hkey << 3   */
constant long next_digits_rf[96] = {
    0, 1, 3, 2, 7, 6, 4, 5,
    0, 7, 1, 6, 3, 4, 2, 5,
    0, 3, 7, 4, 1, 2, 6, 5,
    2, 3, 1, 0, 5, 4, 6, 7,
    4, 3, 5, 2, 7, 0, 6, 1,
    6, 5, 1, 2, 7, 4, 0, 3,
    4, 7, 3, 0, 5, 6, 2, 1,
    6, 7, 5, 4, 1, 0, 2, 3,
    2, 5, 3, 4, 1, 6, 0, 7,
    2, 1, 5, 6, 3, 0, 4, 7,
    4, 5, 7, 6, 3, 2, 0, 1,
    6, 1, 7, 0, 5, 2, 4, 3
};

constant int next_state_rf[96] = {
     1,  2,  3,  2,  4,  5,  3,  5,
     2,  6,  0,  7,  8,  8,  0,  7,
     0,  9, 10,  9,  1,  1, 11, 11,
     6,  0,  6, 11,  9,  0,  9,  8,
    11, 11,  0,  7,  5,  9,  0,  7,
     4,  4,  8,  8,  0,  6, 10,  6,
     5,  7,  5,  3,  1,  1, 11, 11,
     6,  1,  6, 10,  9,  4,  9, 10,
    10,  3,  1,  1, 10,  3,  5,  9,
     4,  4,  8,  8,  2,  7,  2,  3,
     7,  2, 11,  2,  7,  5,  8,  5,
    10,  3,  2,  6, 10,  3,  4,  4
};

/* hilbert_key_rf — direct port of gpu_hilbert.cuf::hilbert_key (device fn).
 * Uses ulong internally (ISHFT is logical; C >> on signed is implementation-
 * defined, so we keep ulong throughout). Max level=21: 3*21=63 bits < 64. */
inline long hilbert_key_rf(int ckey[3], int level)
{
    ulong hkey   = 0;
    int   cstate = 0;
    for (int ibit = level - 1; ibit >= 0; ibit--) {
        hkey = (hkey << 4) >> 1;               /* ISHFT(hkey,4) then ISHFT(hkey,-1) */
        int sdigit = 0;
        if ((ckey[0] >> ibit) & 1) sdigit += 4; /* idim=1: add_digit=2^(3-1)=4 */
        if ((ckey[1] >> ibit) & 1) sdigit += 2; /* idim=2: add_digit=2^(3-2)=2 */
        if ((ckey[2] >> ibit) & 1) sdigit += 1; /* idim=3: add_digit=2^(3-3)=1 */
        int ind = cstate * 8 + sdigit;
        cstate  = next_state_rf[ind];
        hkey   += (ulong)next_digits_rf[ind];
    }
    return as_type<long>(hkey);
}

/* ---- Read-write hash helpers (device atomic_int* for val) ------------- */
/* hash_bucket_r is already defined above in the read-only section.        */

/* hash_set_r: insert or overwrite key→val using 32-bit CAS on hash_val.  */
inline void hash_set_r(device long* hash_key, device atomic_int* hash_val,
                       int hash_size, long key, int val)
{
    int b = hash_bucket_r(key, hash_size);
    for (;;) {
        int expected = 0;
        if (atomic_compare_exchange_weak_explicit(
                &hash_val[b - 1], &expected, val,
                memory_order_relaxed, memory_order_relaxed)) {
            hash_key[b - 1] = key;
            return;
        }
        if (expected == val) { hash_key[b - 1] = key; return; } /* same key re-insert */
        if (expected != 0)    b = (b % hash_size) + 1;          /* occupied → probe   */
        /* expected == 0: spurious weak-CAS failure → retry same slot */
    }
}

/* hash_free_r: set val=0 for key, leaving key in table (matches CUDA).
 * No atomics needed: each thread owns a unique key → unique bucket.
 * The prior kernel dispatch (a full barrier on Metal's serial queue) has
 * completed before this runs, so no concurrent insertion is possible.      */
inline void hash_free_r(device long* hash_key, device int* hash_val,
                        int hash_size, long key)
{
    int b = hash_bucket_r(key, hash_size);
    for (;;) {
        long cur = hash_key[b - 1];
        if (cur == key) { hash_val[b - 1] = 0; return; }
        if (cur == 0L) return;
        b = (b % hash_size) + 1;
    }
}

/* hash_update_r: find key, overwrite val (key must already be in table).
 * No atomics needed: same unique-key-per-thread guarantee as hash_free_r.  */
inline void hash_update_r(device long* hash_key, device int* hash_val,
                          int hash_size, long key, int val)
{
    int b = hash_bucket_r(key, hash_size);
    for (;;) {
        long cur = hash_key[b - 1];
        if (cur == key) { hash_val[b - 1] = val; return; }
        if (cur == 0L) return;
        b = (b % hash_size) + 1;
    }
}

/* =========================================================================
 * refine_kernel — create new child octs for flagged, unrefined cells.
 * Mirrors refine_kernel + make_new_oct (gpu_refine.cuf).
 * Uses MHD interpolation when MHD is enabled and straight injection otherwise.
 *
 * 1D thread layout: gid = oct_offset * 8 + cell_0 (< 8 * num_octs).
 * Each thread independently creates ONE child oct if flag1[cell_0, oct] is
 * set and grid[oct].refined[cell_0] is false.  atomicAdd on ifree_dev gives
 * each creating thread a unique 1-based child index.
 * ========================================================================= */
kernel void refine_kernel(
    device oct_t       *grid     [[buffer(0)]],
    device int         *flag1    [[buffer(1)]],
    device float       *uold     [[buffer(2)]],
#ifdef MHD
    device float       *bold     [[buffer(3)]],
    device atomic_int  *ifree_dev [[buffer(4)]],
    device const long  *hash_key [[buffer(5)]],
    device const int   *hash_val [[buffer(6)]],
    device const int   *ckey_max_dev [[buffer(7)]],
    device const long  *key_off_dev [[buffer(8)]],
    device const int   *box_ckey_min [[buffer(9)]],
    device const int   *box_ckey_max [[buffer(10)]],
    device const int   *periodic [[buffer(11)]],
    constant int       &hash_size [[buffer(12)]],
    constant int       &interpol_var [[buffer(13)]],
    constant int       &interpol_type [[buffer(14)]],
    constant float     &smallr [[buffer(15)]],
    constant int       &head_idx [[buffer(16)]],
    constant int       &num_octs [[buffer(17)]],
#else
    device atomic_int  *ifree_dev [[buffer(3)]],
    constant int       &head_idx [[buffer(4)]],
    constant int       &num_octs [[buffer(5)]],
#endif
#ifdef GRAV
#ifdef MHD
    device float       *f_grav   [[buffer(18)]],
    device float       *phi      [[buffer(19)]],
    device float       *phi_old  [[buffer(20)]],
#else
    device float       *f_grav   [[buffer(6)]],
    device float       *phi      [[buffer(7)]],
    device float       *phi_old  [[buffer(8)]],
#endif
#endif
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs * 8) return;
    int oct_offset = int(gid) / 8;
    int cell_0     = int(gid) % 8;    /* 0-based cell index */
    int oct_abs_0  = (head_idx - 1) + oct_offset;

    bool flagged = (flag1[cell_0 + 8 * oct_abs_0] == 1);
    bool refined = (grid[oct_abs_0].refined[cell_0] != 0);
    if (!flagged || refined) return;

    /* make_new_oct inline (HYDRO=1 only) */
    int i = cell_0 & 1;
    int j = (cell_0 >> 1) & 1;
    int k = (cell_0 >> 2) & 1;

    int ilevel  = grid[oct_abs_0].lev + 1;
    int cckey[3];
    cckey[0] = 2 * grid[oct_abs_0].ckey[0] + i;
    cckey[1] = 2 * grid[oct_abs_0].ckey[1] + j;
    cckey[2] = 2 * grid[oct_abs_0].ckey[2] + k;

    /* Claim the next free slot (1-based). */
    int child_1based = atomic_fetch_add_explicit(ifree_dev, 1, memory_order_relaxed);
    int child_abs_0  = child_1based - 1;

    /* Initialise child oct. */
    grid[child_abs_0].lev    = ilevel;
    grid[child_abs_0].ckey[0] = cckey[0];
    grid[child_abs_0].ckey[1] = cckey[1];
    grid[child_abs_0].ckey[2] = cckey[2];
    for (int c = 0; c < 8; c++) grid[child_abs_0].refined[c] = 0;
    for (int c = 0; c < 8; c++) flag1[c + 8 * child_abs_0] = 0;
    grid[child_abs_0].hkey = hilbert_key_rf(cckey, ilevel - 1);

#ifdef MHD
    float u2[8*NVAR];
    float b2[8*6];
    rf_prolong_mhd(grid, uold, bold, hash_key, hash_val, hash_size,
                   ckey_max_dev, key_off_dev, box_ckey_min, box_ckey_max,
                   periodic, ilevel, cckey, oct_abs_0 + 1, cell_0,
                   interpol_var, interpol_type, smallr, u2, b2);
    for (int c = 0; c < 8; c++) {
        for (int ivar_0 = 0; ivar_0 < NVAR; ivar_0++)
            rf_u_set(uold, child_1based, ivar_0, c, u2[c*NVAR + ivar_0]);
        for (int ib = 0; ib < 6; ib++)
            rf_b_set(bold, child_1based, ib, c, b2[c*6 + ib]);
    }
#else
    for (int ivar_0 = 0; ivar_0 < (NVAR); ivar_0++) {
        float val = uold[cell_0 + 8 * ivar_0 + 8 * (NVAR) * oct_abs_0];
        for (int c = 0; c < 8; c++)
            uold[c + 8 * ivar_0 + 8 * (NVAR) * child_abs_0] = val;
    }
#endif

#ifdef GRAV
    /* Straight injection for gravity variables if active */
    if (f_grav && phi && phi_old) {
        int parent_abs_0 = oct_abs_0;
        for (int c = 0; c < 8; c++) {
            f_grav[(child_abs_0) * 24 + 0 * 8 + c] = f_grav[(parent_abs_0) * 24 + 0 * 8 + cell_0];
            f_grav[(child_abs_0) * 24 + 1 * 8 + c] = f_grav[(parent_abs_0) * 24 + 1 * 8 + cell_0];
            f_grav[(child_abs_0) * 24 + 2 * 8 + c] = f_grav[(parent_abs_0) * 24 + 2 * 8 + cell_0];

            phi[(child_abs_0) * 8 + c] = phi[(parent_abs_0) * 8 + cell_0];
            phi_old[(child_abs_0) * 8 + c] = phi_old[(parent_abs_0) * 8 + cell_0];
        }
    }
#endif

    /* Mark parent cell as refined. */
    grid[oct_abs_0].refined[cell_0] = 1;
}

/* =========================================================================
 * derefine_kernel — remove child octs whose parent cells are no longer flagged.
 * Mirrors derefine_kernel in gpu_refine.cuf. 1 thread/oct.
 * ========================================================================= */
kernel void derefine_kernel(
    device oct_t       *grid        [[buffer(0)]],
    device const int   *flag1       [[buffer(1)]],
    device long        *hash_key    [[buffer(2)]],
    device int         *hash_val    [[buffer(3)]],
    device const int   *ckey_max_dev [[buffer(4)]],
    device const long  *key_off_dev  [[buffer(5)]],
    constant int       &hash_size   [[buffer(6)]],
    constant int       &head_idx    [[buffer(7)]],
    constant int       &num_octs    [[buffer(8)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);
    int ilevel    = grid[oct_abs_0].lev;
    if (ilevel <= 0) return;

    /* Compute oct's own Hilbert key and hash it. */
    long nx  = long(ckey_max_dev[ilevel - 1]);
    long ix  = long(grid[oct_abs_0].ckey[0]);
    long iy  = long(grid[oct_abs_0].ckey[1]);
    long iz  = long(grid[oct_abs_0].ckey[2]);
    long key = key_off_dev[ilevel - 1] + ix + iy * nx + iz * nx * nx;

    /* Compute parent level Hilbert key. */
    int parent_lev = ilevel - 1;
    if (parent_lev <= 0) return;
    long pnx = long(ckey_max_dev[parent_lev - 1]);
    long pix = ix / 2L;
    long piy = iy / 2L;
    long piz = iz / 2L;
    long parent_key = key_off_dev[parent_lev - 1] + pix + piy * pnx + piz * pnx * pnx;

    int parent_1based = hash_get(hash_key, hash_val, hash_size, parent_key);
    if (parent_1based <= 0) return;
    int parent_abs_0 = parent_1based - 1;

    int pi = int(ix - 2L * pix);
    int pj = int(iy - 2L * piy);
    int pk = int(iz - 2L * piz);
    int cell_0 = pi + 2 * pj + 4 * pk;

    bool flagged = (flag1[cell_0 + 8 * parent_abs_0] == 1);
    bool refined = (grid[parent_abs_0].refined[cell_0] != 0);

    if (!flagged && refined) {
        hash_free_r(hash_key, hash_val, hash_size, key);
        grid[oct_abs_0].lev = 0;
        grid[parent_abs_0].refined[cell_0] = 0;
    }
}

/* =========================================================================
 * free_hash_kernel — free hash entries for a range of octs (cache wipe).
 * Mirrors free_hash_kernel in gpu_refine.cuf. 1 thread/oct.
 * ========================================================================= */
kernel void free_hash_kernel(
    device const oct_t *grid        [[buffer(0)]],
    device long        *hash_key    [[buffer(1)]],
    device int         *hash_val    [[buffer(2)]],
    device const int   *ckey_max_dev [[buffer(3)]],
    device const long  *key_off_dev  [[buffer(4)]],
    constant int       &hash_size   [[buffer(5)]],
    constant int       &head_idx    [[buffer(6)]],
    constant int       &num_octs    [[buffer(7)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);
    int ilevel    = grid[oct_abs_0].lev;
    if (ilevel <= 0) return;

    long nx  = long(ckey_max_dev[ilevel - 1]);
    long ix  = long(grid[oct_abs_0].ckey[0]);
    long iy  = long(grid[oct_abs_0].ckey[1]);
    long iz  = long(grid[oct_abs_0].ckey[2]);
    long key = key_off_dev[ilevel - 1] + ix + iy * nx + iz * nx * nx;
    hash_free_r(hash_key, hash_val, hash_size, key);
}

/* =========================================================================
 * update_hash_kernel — update val for each oct's key after sort rearrangement.
 * Mirrors update_hash_kernel in gpu_refine.cuf. 1 thread/oct.
 * ========================================================================= */
kernel void update_hash_kernel(
    device const oct_t *grid        [[buffer(0)]],
    device long        *hash_key    [[buffer(1)]],
    device int         *hash_val    [[buffer(2)]],
    device const int   *ckey_max_dev [[buffer(3)]],
    device const long  *key_off_dev  [[buffer(4)]],
    constant int       &hash_size   [[buffer(5)]],
    constant int       &head_idx    [[buffer(6)]],
    constant int       &num_octs    [[buffer(7)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);
    int ilevel    = grid[oct_abs_0].lev;
    if (ilevel <= 0) return;

    long nx  = long(ckey_max_dev[ilevel - 1]);
    long ix  = long(grid[oct_abs_0].ckey[0]);
    long iy  = long(grid[oct_abs_0].ckey[1]);
    long iz  = long(grid[oct_abs_0].ckey[2]);
    long key = key_off_dev[ilevel - 1] + ix + iy * nx + iz * nx * nx;
    int  val = (head_idx - 1) + int(gid) + 1;   /* 1-based oct index */
    hash_update_r(hash_key, hash_val, hash_size, key, val);
}

/* =========================================================================
 * insert_hash_all_kernel — insert all octs in a range, reading level from grid.
 * Used for newly created octs (which may span multiple levels).
 * Mirrors insert_hash_kernel in gpu_refine.cuf (the multi-level version).
 * ========================================================================= */
kernel void insert_hash_all_kernel(
    device const oct_t *grid        [[buffer(0)]],
    device long        *hash_key    [[buffer(1)]],
    device atomic_int  *hash_val    [[buffer(2)]],
    device const int   *ckey_max_dev [[buffer(3)]],
    device const long  *key_off_dev  [[buffer(4)]],
    constant int       &hash_size   [[buffer(5)]],
    constant int       &head_idx    [[buffer(6)]],
    constant int       &num_octs    [[buffer(7)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);
    int ilevel    = grid[oct_abs_0].lev;
    if (ilevel <= 0) return;

    long nx  = long(ckey_max_dev[ilevel - 1]);
    long ix  = long(grid[oct_abs_0].ckey[0]);
    long iy  = long(grid[oct_abs_0].ckey[1]);
    long iz  = long(grid[oct_abs_0].ckey[2]);
    long key = key_off_dev[ilevel - 1] + ix + iy * nx + iz * nx * nx;
    int  val = oct_abs_0 + 1;   /* 1-based */
    hash_set_r(hash_key, hash_val, hash_size, key, val);
}

/* =========================================================================
 * enforce_rules_hash_kernel — same as ensure_ref_rules on the CPU, with hash
 * lookups instead of nbor, so that the octs created earlier in this refine
 * pass count. Checks that all octs at level ilevel are surrounded by 26
 * neighbouring octs, for entire groups of nskip=(nsubgrid/2)**3 octs.
 * Mirrors enforce_rules_hash in gpu_refine.cuf. 1 thread/group.
 * ========================================================================= */
kernel void enforce_rules_hash_kernel(
    device const oct_t *grid         [[buffer(0)]],
    device int         *flag1        [[buffer(1)]],
    device const long  *hash_key     [[buffer(2)]],
    device const int   *hash_val     [[buffer(3)]],
    device const int   *ckey_max_dev [[buffer(4)]],
    device const long  *key_off_dev  [[buffer(5)]],
    device const int   *box_ckey_min [[buffer(6)]],
    device const int   *box_ckey_max [[buffer(7)]],
    device const int   *periodic     [[buffer(8)]],
    constant int       &hash_size    [[buffer(9)]],
    constant int       &ngridmax     [[buffer(10)]],
    constant int       &nskip        [[buffer(11)]],
    constant int       &head_idx     [[buffer(12)]],
    constant int       &num_groups   [[buffer(13)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_groups) return;
    int oct_abs_0 = (head_idx - 1) + int(gid) * nskip;   /* leading oct of the group */

    bool ok = true;

    /* Loop over octs in the group */
    for (int iskip = 0; iskip < nskip; iskip++) {
        int ilevel = grid[oct_abs_0 + iskip].lev;
        long nx  = long(ckey_max_dev[ilevel - 1]);
        long off = key_off_dev[ilevel - 1];

        /* Loop over 3x3x3 neighbouring father cells */
        for (int kn = -1; kn <= 1; kn++) {
            for (int jn = -1; jn <= 1; jn++) {
                for (int in = -1; in <= 1; in++) {
                    int ckey_n[3];
                    ckey_n[0] = grid[oct_abs_0 + iskip].ckey[0] + in;
                    ckey_n[1] = grid[oct_abs_0 + iskip].ckey[1] + jn;
                    ckey_n[2] = grid[oct_abs_0 + iskip].ckey[2] + kn;

                    /* Periodic boundary conditions (level-specific box bounds). */
                    for (int d = 0; d < 3; d++) {
                        if (periodic[d]) {
                            int bmin = box_ckey_min[3 * (ilevel - 1) + d];
                            int bmax = box_ckey_max[3 * (ilevel - 1) + d];
                            if (ckey_n[d] <  bmin) ckey_n[d] = bmax - 1;
                            if (ckey_n[d] >= bmax) ckey_n[d] = bmin;
                        }
                    }

                    long key = off + long(ckey_n[0])
                                   + long(ckey_n[1]) * nx
                                   + long(ckey_n[2]) * nx * nx;
                    int nbor_idx = hash_get(hash_key, hash_val, hash_size, key);
                    ok = ok && (nbor_idx > 0) && (nbor_idx <= ngridmax);
                }
            }
        }
    }

    /* Only cancel new refinements: unflagging a refined cell would
       destroy its child and could orphan finer octs (2:1 violation) */
    if (!ok) {
        for (int iskip = 0; iskip < nskip; iskip++) {
            for (int ic = 0; ic < 8; ic++) {
                if (!(bool)grid[oct_abs_0 + iskip].refined[ic]) flag1[ic + 8 * (oct_abs_0 + iskip)] = 0;
            }
        }
    }
}

/* =========================================================================
 * init_global_swap_table_kernel — swap_global[oct] = oct (identity).
 * Mirrors init_global_swap_table_kernel in gpu_refine.cuf. 1 thread/oct.
 * ========================================================================= */
kernel void init_global_swap_table_kernel(
    device int   *swap_global [[buffer(0)]],
    constant int &head_idx    [[buffer(1)]],
    constant int &num_octs    [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);
    swap_global[oct_abs_0] = oct_abs_0 + 1;   /* 1-based */
}

/* =========================================================================
 * init_prefix_sum_level_kernel
 * prefix_sum[oct] = 0 if grid[swap_global[oct]-1].lev == ilevel, else 1.
 * Mirrors init_prefix_sum_level in gpu_refine.cuf. 1 thread/oct.
 * ========================================================================= */
kernel void init_prefix_sum_level_kernel(
    device const oct_t *grid        [[buffer(0)]],
    device const int   *swap_global [[buffer(1)]],
    device int         *prefix_sum  [[buffer(2)]],
    constant int       &head_idx    [[buffer(3)]],
    constant int       &num_octs    [[buffer(4)]],
    constant int       &ilevel      [[buffer(5)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0     = (head_idx - 1) + int(gid);
    int old_1based    = swap_global[oct_abs_0];
    int old_abs_0     = old_1based - 1;
    prefix_sum[oct_abs_0] = (grid[old_abs_0].lev == ilevel) ? 0 : 1;
}

/* =========================================================================
 * init_prefix_sum_bit_kernel
 * prefix_sum[oct] = bit `ibit` of grid[swap_global[oct]-1].hkey.
 * Mirrors init_prefix_sum_bit in gpu_refine.cuf. 1 thread/oct.
 * ========================================================================= */
kernel void init_prefix_sum_bit_kernel(
    device const oct_t *grid        [[buffer(0)]],
    device const int   *swap_global [[buffer(1)]],
    device int         *prefix_sum  [[buffer(2)]],
    constant int       &head_idx    [[buffer(3)]],
    constant int       &num_octs    [[buffer(4)]],
    constant int       &ibit        [[buffer(5)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0  = (head_idx - 1) + int(gid);
    int old_1based = swap_global[oct_abs_0];
    int old_abs_0  = old_1based - 1;
    ulong hkey     = as_type<ulong>(grid[old_abs_0].hkey);
    prefix_sum[oct_abs_0] = int((hkey >> ibit) & 1UL);
}

/* =========================================================================
 * compute_local_swap_table_kernel — LSD scatter step.
 * Mirrors compute_local_swap_table_impl in gpu_refine.cuf. 1 thread/oct.
 * ========================================================================= */
kernel void compute_local_swap_table_kernel(
    device int       *swap_local  [[buffer(0)]],
    device const int *swap_global [[buffer(1)]],
    device const int *prefix_sum  [[buffer(2)]],
    constant int     &head_idx    [[buffer(3)]],
    constant int     &num_octs    [[buffer(4)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);
    int tail_abs  = (head_idx - 1) + num_octs - 1;

    int nones  = prefix_sum[tail_abs];
    int nzeros = num_octs - nones;

    int prev_sum = (oct_abs_0 > (head_idx - 1)) ? prefix_sum[oct_abs_0 - 1] : 0;
    int bit      = prefix_sum[oct_abs_0] - prev_sum;
    int old_val  = swap_global[oct_abs_0];   /* 1-based oct index */

    int sorted_abs_0;
    if (bit == 0)
        sorted_abs_0 = oct_abs_0 - prev_sum;
    else
        sorted_abs_0 = (head_idx - 1) + nzeros + prev_sum;

    swap_local[sorted_abs_0] = old_val;
}

/* =========================================================================
 * update_global_swap_table_kernel — swap_global[oct] = swap_local[oct].
 * Mirrors update_global_swap_table_impl. 1 thread/oct.
 * ========================================================================= */
kernel void update_global_swap_table_kernel(
    device int       *swap_global [[buffer(0)]],
    device const int *swap_local  [[buffer(1)]],
    constant int     &head_idx    [[buffer(2)]],
    constant int     &num_octs    [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);
    swap_global[oct_abs_0] = swap_local[oct_abs_0];
}

/* =========================================================================
 * sort_gather_grid_kernel — pack grid properties into flag2 scratch.
 * flag2 slots (0-based, per oct):
 *   [0] = lev, [1]=ckey[0], [2]=ckey[1], [3]=ckey[2], [4]=refined_mask
 * Mirrors sort_gather_grid_impl. 1 thread/oct.
 * ========================================================================= */
kernel void sort_gather_grid_kernel(
    device int       *flag2       [[buffer(0)]],
    device const oct_t *grid      [[buffer(1)]],
    device const int *swap_global [[buffer(2)]],
    constant int     &head_idx    [[buffer(3)]],
    constant int     &num_octs    [[buffer(4)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0     = (head_idx - 1) + int(gid);
    int old_1based    = swap_global[oct_abs_0];
    int old_abs_0     = old_1based - 1;

    flag2[0 + 8 * oct_abs_0] = grid[old_abs_0].lev;
    flag2[1 + 8 * oct_abs_0] = grid[old_abs_0].ckey[0];
    flag2[2 + 8 * oct_abs_0] = grid[old_abs_0].ckey[1];
    flag2[3 + 8 * oct_abs_0] = grid[old_abs_0].ckey[2];
    int mask = 0;
    for (int c = 0; c < 8; c++)
        if (grid[old_abs_0].refined[c]) mask |= (1 << c);
    flag2[4 + 8 * oct_abs_0] = mask;
}

/* =========================================================================
 * sort_scatter_grid_kernel — unpack flag2 scratch back to grid, recompute hkey.
 * Mirrors sort_scatter_grid_impl. 1 thread/oct.
 * ========================================================================= */
kernel void sort_scatter_grid_kernel(
    device oct_t     *grid    [[buffer(0)]],
    device const int *flag2   [[buffer(1)]],
    constant int     &head_idx [[buffer(2)]],
    constant int     &num_octs [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);

    int ilevel    = flag2[0 + 8 * oct_abs_0];
    int ckey[3];
    ckey[0] = flag2[1 + 8 * oct_abs_0];
    ckey[1] = flag2[2 + 8 * oct_abs_0];
    ckey[2] = flag2[3 + 8 * oct_abs_0];
    int mask      = flag2[4 + 8 * oct_abs_0];

    grid[oct_abs_0].lev    = ilevel;
    grid[oct_abs_0].ckey[0] = ckey[0];
    grid[oct_abs_0].ckey[1] = ckey[1];
    grid[oct_abs_0].ckey[2] = ckey[2];
    for (int c = 0; c < 8; c++)
        grid[oct_abs_0].refined[c] = (mask >> c) & 1;
    grid[oct_abs_0].hkey = hilbert_key_rf(ckey, ilevel - 1);
}

/* =========================================================================
 * sort_gather_flag_kernel — flag2[:,oct] = flag1[:,swap_global[oct]-1].
 * Mirrors sort_gather_flag_impl. 1 thread/oct.
 * ========================================================================= */
kernel void sort_gather_flag_kernel(
    device int       *flag2       [[buffer(0)]],
    device const int *flag1       [[buffer(1)]],
    device const int *swap_global [[buffer(2)]],
    constant int     &head_idx    [[buffer(3)]],
    constant int     &num_octs    [[buffer(4)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0  = (head_idx - 1) + int(gid);
    int old_abs_0  = swap_global[oct_abs_0] - 1;
    for (int c = 0; c < 8; c++)
        flag2[c + 8 * oct_abs_0] = flag1[c + 8 * old_abs_0];
}

/* =========================================================================
 * sort_scatter_flag_kernel — flag1[:,oct] = flag2[:,oct].
 * Mirrors sort_scatter_flag_impl. 1 thread/oct.
 * ========================================================================= */
kernel void sort_scatter_flag_kernel(
    device int       *flag1   [[buffer(0)]],
    device const int *flag2   [[buffer(1)]],
    constant int     &head_idx [[buffer(2)]],
    constant int     &num_octs [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_abs_0 = (head_idx - 1) + int(gid);
    for (int c = 0; c < 8; c++)
        flag1[c + 8 * oct_abs_0] = flag2[c + 8 * oct_abs_0];
}

/* =========================================================================
 * sort_gather_hydro_kernel — gather hydro and optional magnetic state.
 * Mirrors sort_gather_hydro_impl.
 * 1D gid = oct_offset * 8 + cell_0; inner loop over NVAR variables.
 * ========================================================================= */
kernel void sort_gather_hydro_kernel(
    device float     *unew        [[buffer(0)]],
    device const float *uold      [[buffer(1)]],
#ifdef MHD
    device float     *bnew        [[buffer(2)]],
    device const float *bold      [[buffer(3)]],
    device const int *swap_global [[buffer(4)]],
    constant int     &head_idx    [[buffer(5)]],
    constant int     &num_octs    [[buffer(6)]],
#else
    device const int *swap_global [[buffer(2)]],
    constant int     &head_idx    [[buffer(3)]],
    constant int     &num_octs    [[buffer(4)]],
#endif
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs * 8) return;
    int oct_offset = int(gid) / 8;
    int cell_0     = int(gid) % 8;
    int oct_abs_0  = (head_idx - 1) + oct_offset;
    int old_abs_0  = swap_global[oct_abs_0] - 1;

    for (int ivar_0 = 0; ivar_0 < (NVAR); ivar_0++)
        unew[cell_0 + 8 * ivar_0 + 8 * (NVAR) * oct_abs_0] =
            uold[cell_0 + 8 * ivar_0 + 8 * (NVAR) * old_abs_0];
#ifdef MHD
    for (int ib = 0; ib < 6; ib++)
        bnew[cell_0 + 8 * ib + 8 * 6 * oct_abs_0] =
            bold[cell_0 + 8 * ib + 8 * 6 * old_abs_0];
#endif
}

/* =========================================================================
 * update_nbor_prefix_kernel — compute nbor[input_ind, subgrid] and set
 * prefix_sum[subgrid] = 1 if missing, 0 if found.
 * Mirrors update_nbor_array(..., prefix_sum=prefix_sum) in gpu_refine.cuf.
 * 1 thread/subgrid.
 * ========================================================================= */
kernel void update_nbor_prefix_kernel(
    device int         *nbor          [[buffer(0)]],
    device const oct_t *grid          [[buffer(1)]],
    device const long  *hash_key      [[buffer(2)]],
    device const int   *hash_val      [[buffer(3)]],
    device const int   *ckey_max_dev  [[buffer(4)]],
    device const long  *key_off_dev   [[buffer(5)]],
    device const int   *box_ckey_min  [[buffer(6)]],
    device const int   *box_ckey_max  [[buffer(7)]],
    device const int   *periodic      [[buffer(8)]],
    device int         *prefix_sum    [[buffer(9)]],
    constant int       &hash_size     [[buffer(10)]],
    constant int       &head_idx      [[buffer(11)]],
    constant int       &num_subgrids  [[buffer(12)]],
    constant int       &input_ind     [[buffer(13)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_subgrids) return;
    int subgrid_abs_0 = (head_idx - 1) + int(gid);
    int oct_abs_0     = subgrid_abs_0 * NSUBGRIDTONDIM_RF;

    int ind0 = input_ind - 1;
    int k_off = ind0 / (NSUBGRIDP2_RF * NSUBGRIDP2_RF);
    int j_off = (ind0 - k_off * NSUBGRIDP2_RF * NSUBGRIDP2_RF) / NSUBGRIDP2_RF;
    int i_off = ind0 - k_off * NSUBGRIDP2_RF * NSUBGRIDP2_RF - j_off * NSUBGRIDP2_RF;

    int ilevel = grid[oct_abs_0].lev;
    int ckey_n[3];
    ckey_n[0] = (grid[oct_abs_0].ckey[0] / NSUBGRID_RF) * NSUBGRID_RF + i_off - 1;
    ckey_n[1] = (grid[oct_abs_0].ckey[1] / NSUBGRID_RF) * NSUBGRID_RF + j_off - 1;
    ckey_n[2] = (grid[oct_abs_0].ckey[2] / NSUBGRID_RF) * NSUBGRID_RF + k_off - 1;

    /* Periodic boundary conditions (level-specific box bounds). */
    for (int d = 0; d < 3; d++) {
        if (periodic[d]) {
            int bmin = box_ckey_min[3 * (ilevel - 1) + d];
            int bmax = box_ckey_max[3 * (ilevel - 1) + d];
            if (ckey_n[d] <  bmin) ckey_n[d] = bmax - 1;
            if (ckey_n[d] >= bmax) ckey_n[d] = bmin;
        }
    }

    long nx  = long(ckey_max_dev[ilevel - 1]);
    long key = key_off_dev[ilevel - 1]
               + long(ckey_n[0])
               + long(ckey_n[1]) * nx
               + long(ckey_n[2]) * nx * nx;

    int nbor_val = hash_get(hash_key, hash_val, hash_size, key);
    nbor[subgrid_abs_0 * SUBGRIDSIZE_RF + (input_ind - 1)] = nbor_val;
    prefix_sum[subgrid_abs_0] = (nbor_val == 0) ? 1 : 0;
}

/* =========================================================================
 * compute_cache_swap_table_kernel — build swap_local for cache oct creation.
 * If prefix_sum bit is 1 (nbor missing), write subgrid_idx (1-based) to
 * swap_local[prev_sum] (0-based slot).
 * Mirrors compute_cache_swap_table_impl in gpu_refine.cuf. 1 thread/subgrid.
 * ========================================================================= */
kernel void compute_cache_swap_table_kernel(
    device int       *swap_local   [[buffer(0)]],
    device const int *prefix_sum   [[buffer(1)]],
    constant int     &head_idx     [[buffer(2)]],
    constant int     &num_subgrids [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_subgrids) return;
    int subgrid_abs_0 = (head_idx - 1) + int(gid);
    int prev_sum = (subgrid_abs_0 > (head_idx - 1)) ? prefix_sum[subgrid_abs_0 - 1] : 0;
    int bit      = prefix_sum[subgrid_abs_0] - prev_sum;
    if (bit == 1)
        swap_local[prev_sum] = subgrid_abs_0 + 1;   /* 1-based subgrid idx */
}

/* =========================================================================
 * make_cache_octs_kernel — create ghost octs in the cache region.
 * Mirrors make_cache_octs_impl.
 * 1 thread per new cache oct (gid < new_noct).
 * ========================================================================= */
kernel void make_cache_octs_kernel(
    device oct_t       *grid          [[buffer(0)]],
    device int         *flag1         [[buffer(1)]],
    device float       *uold          [[buffer(2)]],
#ifdef MHD
    device float       *bold          [[buffer(3)]],
    device const int   *swap_local    [[buffer(4)]],
    device int         *father        [[buffer(5)]],
    device int         *nbor          [[buffer(6)]],
    device const long  *hash_key      [[buffer(7)]],
    device const int   *hash_val      [[buffer(8)]],
    device const int   *ckey_max_dev  [[buffer(9)]],
    device const long  *key_off_dev   [[buffer(10)]],
    device const int   *box_ckey_min  [[buffer(11)]],
    device const int   *box_ckey_max  [[buffer(12)]],
    device const int   *periodic      [[buffer(13)]],
    constant int       &hash_size     [[buffer(14)]],
    constant int       &ngridmax      [[buffer(15)]],
    constant int       &ifree_cache   [[buffer(16)]],
    constant int       &new_noct      [[buffer(17)]],
    constant int       &input_ind     [[buffer(18)]],
    constant int       &interpol_var  [[buffer(19)]],
    constant int       &interpol_type [[buffer(20)]],
    constant float     &smallr        [[buffer(21)]],
#else
    device const int   *swap_local    [[buffer(3)]],
    device int         *father        [[buffer(4)]],
    device int         *nbor          [[buffer(5)]],
    device const long  *hash_key      [[buffer(6)]],
    device const int   *hash_val      [[buffer(7)]],
    device const int   *ckey_max_dev  [[buffer(8)]],
    device const long  *key_off_dev   [[buffer(9)]],
    device const int   *box_ckey_min  [[buffer(10)]],
    device const int   *box_ckey_max  [[buffer(11)]],
    device const int   *periodic      [[buffer(12)]],
    constant int       &hash_size     [[buffer(13)]],
    constant int       &ngridmax      [[buffer(14)]],
    constant int       &ifree_cache   [[buffer(15)]],   /* 1-based offset (host-read) */
    constant int       &new_noct      [[buffer(16)]],
    constant int       &input_ind     [[buffer(17)]],
#endif
#ifdef GRAV
#ifdef MHD
    device float       *f_grav        [[buffer(22)]],
    device float       *phi           [[buffer(23)]],
    device float       *phi_old       [[buffer(24)]],
#else
    device float       *f_grav        [[buffer(18)]],
    device float       *phi           [[buffer(19)]],
    device float       *phi_old       [[buffer(20)]],
#endif
#endif
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= new_noct) return;

    /* swap_local[gid] holds the 1-based subgrid idx (set by compute_cache_swap_table). */
    int subgrid_1based = swap_local[int(gid)];
    int subgrid_abs_0  = subgrid_1based - 1;
    int oct_abs_0      = subgrid_abs_0 * NSUBGRIDTONDIM_RF;

    /* Cache grid slot: 1-based = ngridmax + ifree_cache + gid */
    int cache_1based = ngridmax + ifree_cache + int(gid);
    int cache_abs_0  = cache_1based - 1;

    nbor[subgrid_abs_0 * SUBGRIDSIZE_RF + (input_ind - 1)] = cache_1based;

    /* Decode 3D offset for the neighbour direction. */
    int ind0  = input_ind - 1;
    int k_off = ind0 / (NSUBGRIDP2_RF * NSUBGRIDP2_RF);
    int j_off = (ind0 - k_off * NSUBGRIDP2_RF * NSUBGRIDP2_RF) / NSUBGRIDP2_RF;
    int i_off = ind0 - k_off * NSUBGRIDP2_RF * NSUBGRIDP2_RF - j_off * NSUBGRIDP2_RF;

    int ilevel = grid[oct_abs_0].lev;
    int ckey[3];
    ckey[0] = (grid[oct_abs_0].ckey[0] / NSUBGRID_RF) * NSUBGRID_RF + i_off - 1;
    ckey[1] = (grid[oct_abs_0].ckey[1] / NSUBGRID_RF) * NSUBGRID_RF + j_off - 1;
    ckey[2] = (grid[oct_abs_0].ckey[2] / NSUBGRID_RF) * NSUBGRID_RF + k_off - 1;

    for (int d = 0; d < 3; d++) {
        if (periodic[d]) {
            int bmin = box_ckey_min[3 * (ilevel - 1) + d];
            int bmax = box_ckey_max[3 * (ilevel - 1) + d];
            if (ckey[d] <  bmin) ckey[d] = bmax - 1;
            if (ckey[d] >= bmax) ckey[d] = bmin;
        }
    }

    /* Set cache oct grid properties. */
    grid[cache_abs_0].lev    = ilevel;
    grid[cache_abs_0].ckey[0] = ckey[0];
    grid[cache_abs_0].ckey[1] = ckey[1];
    grid[cache_abs_0].ckey[2] = ckey[2];
    for (int c = 0; c < 8; c++) grid[cache_abs_0].refined[c] = 0;
    for (int c = 0; c < 8; c++) flag1[c + 8 * cache_abs_0] = 0;
    grid[cache_abs_0].hkey = hilbert_key_rf(ckey, ilevel - 1);

    /* Compute parent hash key. */
    int parent_lev = ilevel - 1;
    long pnx = long(ckey_max_dev[parent_lev - 1]);
    long pix = long(ckey[0]) / 2L;
    long piy = long(ckey[1]) / 2L;
    long piz = long(ckey[2]) / 2L;
    long parent_key = key_off_dev[parent_lev - 1] + pix + piy * pnx + piz * pnx * pnx;

    int parent_1based = hash_get(hash_key, hash_val, hash_size, parent_key);
    father[cache_abs_0] = parent_1based;

    if (parent_1based > 0) {
        int pi = int(long(ckey[0]) - 2L * pix);
        int pj = int(long(ckey[1]) - 2L * piy);
        int pk = int(long(ckey[2]) - 2L * piz);
        int cell_0 = pi + 2 * pj + 4 * pk;
#ifdef MHD
        float u2[8*NVAR];
        float b2[8*6];
        rf_prolong_mhd(grid, uold, bold, hash_key, hash_val, hash_size,
                       ckey_max_dev, key_off_dev, box_ckey_min, box_ckey_max,
                       periodic, ilevel, ckey, parent_1based, cell_0,
                       interpol_var, interpol_type, smallr, u2, b2);
        for (int c = 0; c < 8; c++) {
            for (int ivar_0 = 0; ivar_0 < NVAR; ivar_0++)
                rf_u_set(uold, cache_1based, ivar_0, c, u2[c*NVAR + ivar_0]);
            for (int ib = 0; ib < 6; ib++)
                rf_b_set(bold, cache_1based, ib, c, b2[c*6 + ib]);
        }
#else
        int parent_abs_0 = parent_1based - 1;
        for (int ivar_0 = 0; ivar_0 < (NVAR); ivar_0++) {
            float val = uold[cell_0 + 8 * ivar_0 + 8 * (NVAR) * parent_abs_0];
            for (int c = 0; c < 8; c++)
                uold[c + 8 * ivar_0 + 8 * (NVAR) * cache_abs_0] = val;
        }
#ifdef GRAV
        if (f_grav && phi && phi_old) {
            int cell_idx_1based = cell_0 + 1;
            for (int c = 1; c <= 8; c++) {
                f_grav[(cache_1based - 1) * 24 + 0 * 8 + (c - 1)] = f_grav[(parent_1based - 1) * 24 + 0 * 8 + (cell_idx_1based - 1)];
                f_grav[(cache_1based - 1) * 24 + 1 * 8 + (c - 1)] = f_grav[(parent_1based - 1) * 24 + 1 * 8 + (cell_idx_1based - 1)];
                f_grav[(cache_1based - 1) * 24 + 2 * 8 + (c - 1)] = f_grav[(parent_1based - 1) * 24 + 2 * 8 + (cell_idx_1based - 1)];

                phi[(cache_1based - 1) * 8 + (c - 1)] = phi[(parent_1based - 1) * 8 + (cell_idx_1based - 1)];
                phi_old[(cache_1based - 1) * 8 + (c - 1)] = phi_old[(parent_1based - 1) * 8 + (cell_idx_1based - 1)];
            }
        }
#endif
#endif
    }
}

/* =========================================================================
 * insert_hash_cache_kernel — insert new cache octs into hash table.
 * Mirrors insert_hash_cache_impl in gpu_refine.cuf. 1 thread per cache oct.
 * ========================================================================= */
kernel void insert_hash_cache_kernel(
    device const oct_t *grid        [[buffer(0)]],
    device long        *hash_key    [[buffer(1)]],
    device atomic_int  *hash_val    [[buffer(2)]],
    device const int   *ckey_max_dev [[buffer(3)]],
    device const long  *key_off_dev  [[buffer(4)]],
    constant int       &hash_size   [[buffer(5)]],
    constant int       &ngridmax    [[buffer(6)]],
    constant int       &ifree_cache [[buffer(7)]],
    constant int       &new_noct    [[buffer(8)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= new_noct) return;
    int cache_1based = ngridmax + ifree_cache + int(gid);
    int cache_abs_0  = cache_1based - 1;
    int ilevel       = grid[cache_abs_0].lev;
    if (ilevel <= 0) return;

    long nx  = long(ckey_max_dev[ilevel - 1]);
    long ix  = long(grid[cache_abs_0].ckey[0]);
    long iy  = long(grid[cache_abs_0].ckey[1]);
    long iz  = long(grid[cache_abs_0].ckey[2]);
    long key = key_off_dev[ilevel - 1] + ix + iy * nx + iz * nx * nx;
    hash_set_r(hash_key, hash_val, hash_size, key, cache_1based);
}

/* =========================================================================
 * Gravity/Poisson variable sorting kernels (for refinement/reordering)
 * ========================================================================= */
kernel void sort_gather_force_kernel(
    device float       *scratch      [[buffer(0)]], // s_nref
    device const float *f_grav       [[buffer(1)]], // s_f_grav
    device const int   *swap_global  [[buffer(2)]], // s_swap_global
    constant int       &idim         [[buffer(3)]], // 1-based dimension
    constant int       &head_idx     [[buffer(4)]], // 1-based start
    constant int       &num_octs     [[buffer(5)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_idx = head_idx + int(gid);
    int old_oct_idx = swap_global[oct_idx - 1];

    for (int c = 0; c < 8; c++) {
        scratch[(oct_idx - 1) * 8 + c] = f_grav[(old_oct_idx - 1) * 24 + (idim - 1) * 8 + c];
    }
}

kernel void sort_scatter_force_kernel(
    device float       *f_grav       [[buffer(0)]], // s_f_grav
    device const float *scratch      [[buffer(1)]], // s_nref
    constant int       &idim         [[buffer(2)]], // 1-based dimension
    constant int       &head_idx     [[buffer(3)]], // 1-based start
    constant int       &num_octs     [[buffer(4)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_idx = head_idx + int(gid);

    for (int c = 0; c < 8; c++) {
        f_grav[(oct_idx - 1) * 24 + (idim - 1) * 8 + c] = scratch[(oct_idx - 1) * 8 + c];
    }
}

kernel void sort_gather_phi_kernel(
    device float       *scratch      [[buffer(0)]], // s_nref
    device const float *phi          [[buffer(1)]], // s_phi or s_phi_old
    device const int   *swap_global  [[buffer(2)]], // s_swap_global
    constant int       &head_idx     [[buffer(3)]], // 1-based start
    constant int       &num_octs     [[buffer(4)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_idx = head_idx + int(gid);
    int old_oct_idx = swap_global[oct_idx - 1];

    for (int c = 0; c < 8; c++) {
        scratch[(oct_idx - 1) * 8 + c] = phi[(old_oct_idx - 1) * 8 + c];
    }
}

kernel void sort_scatter_phi_kernel(
    device float       *phi          [[buffer(0)]], // s_phi or s_phi_old
    device const float *scratch      [[buffer(1)]], // s_nref
    constant int       &head_idx     [[buffer(2)]], // 1-based start
    constant int       &num_octs     [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if (int(gid) >= num_octs) return;
    int oct_idx = head_idx + int(gid);

    for (int c = 0; c < 8; c++) {
        phi[(oct_idx - 1) * 8 + c] = scratch[(oct_idx - 1) * 8 + c];
    }
}

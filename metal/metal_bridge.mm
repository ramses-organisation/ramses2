#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mach/mach.h>

/* Returns resident set size in pages (4096 bytes each), matching Linux /proc/self/stat field 24. */
extern "C" long getmem_mac(void)
{
    struct task_basic_info info;
    mach_msg_type_number_t count = TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&info, &count) != KERN_SUCCESS)
        return 0L;
    return (long)(info.resident_size / 4096);
}

#include "metal_types.h"
#include "metal_config.h"

/* -----------------------------------------------------------------------
 * Dispatch helpers — thread layout macros.
 * ----------------------------------------------------------------------- */
#define DISPATCH_2D_8_16(enc_, n_) \
    [enc_ dispatchThreadgroups:{((NSUInteger)(n_)+15)/16,1,1} threadsPerThreadgroup:{8,16,1}]

#define DISPATCH_1D_256_OCT(enc_, n_) \
    [enc_ dispatchThreadgroups:{((NSUInteger)(n_)*8+255)/256,1,1} threadsPerThreadgroup:{256,1,1}]

#define DISPATCH_1D_1024_OCT(enc_, n_) \
    [enc_ dispatchThreadgroups:{((NSUInteger)(n_)*8+1023)/1024,1,1} threadsPerThreadgroup:{1024,1,1}]

#define DISPATCH_1D_128(enc_, n_) \
    [enc_ dispatchThreadgroups:{((NSUInteger)(n_)+127)/128,1,1} threadsPerThreadgroup:{128,1,1}]

#define DISPATCH_2D_4_32(enc_, n_) \
    [enc_ dispatchThreadgroups:{((NSUInteger)(n_)+31)/32,1,1} threadsPerThreadgroup:{4,32,1}]

/* -----------------------------------------------------------------------
 * File-scope Metal objects retained for the program lifetime.
 * Mirrors the module-level device arrays in gpu_manager.cuf.
 * ----------------------------------------------------------------------- */
static id<MTLDevice>               s_device       = nil;
static id<MTLCommandQueue>         s_queue        = nil;
static id<MTLLibrary>              s_library      = nil;
static id<MTLBuffer>               s_uold         = nil;
static id<MTLBuffer>               s_unew         = nil;
static id<MTLBuffer>               s_bold         = nil;
static id<MTLBuffer>               s_bnew         = nil;
static id<MTLBuffer>               s_uct_velocity = nil;
static id<MTLBuffer>               s_uct_face_product = nil;
static bool                        s_uct_face_reuse_disabled = false;
static id<MTLBuffer>               s_grid         = nil;
static id<MTLBuffer>               s_nbor         = nil;
static id<MTLBuffer>               s_hash_key     = nil;
static id<MTLBuffer>               s_hash_val     = nil;
/* AMR refinement buffers (allocated by mtl_alloc_refine) */
static id<MTLBuffer>               s_flag1            = nil;
static id<MTLBuffer>               s_flag2            = nil;
static id<MTLBuffer>               s_father           = nil;
static id<MTLBuffer>               s_swap_local       = nil;
static id<MTLBuffer>               s_swap_global      = nil;
static id<MTLBuffer>               s_prefix_sum       = nil;
static id<MTLBuffer>               s_partial_sums     = nil;  /* level-1 scratch (ps0): ceil(n/256)      */
static id<MTLBuffer>               s_partial_sums_2   = nil;  /* level-2 scratch (ps1): ceil(n/256^2)    */
static id<MTLBuffer>               s_partial_sums_3   = nil;  /* level-3 scratch (ps2): ceil(n/256^3)    */
static id<MTLBuffer>               s_partial_sums_4   = nil;  /* dummy sink for deepest single-block pass */
static id<MTLBuffer>               s_ifree_dev        = nil;
static id<MTLBuffer>               s_ifree_cache_dev  = nil;
static id<MTLBuffer>               s_ckey_max_dev     = nil;
static id<MTLBuffer>               s_key_off_dev      = nil;
static id<MTLBuffer>               s_box_ckey_min_dev = nil;
static id<MTLBuffer>               s_box_ckey_max_dev = nil;
static id<MTLBuffer>               s_periodic_dev     = nil;

/* Fence for inter-encoder ordering within a single command buffer.
 * Used by mtl_hilbert_sort_level to synchronize between 6 encoders per pass
 * without needing 6 separate commit/waitUntilCompleted cycles.            */
static id<MTLFence>                s_sort_fence       = nil;
static id<MTLBuffer>               s_count_buf        = nil;  /* 1×int, persistent for batched flag reductions */
static int                         s_ngridmax         = 0;
static int                         s_hash_size        = 0;

static id<MTLComputePipelineState> s_pso_set_unew = nil;
static id<MTLComputePipelineState> s_pso_set_uold = nil;
static id<MTLComputePipelineState> s_pso_cmpdt    = nil;
static id<MTLComputePipelineState> s_pso_sync_hydro = nil;
static id<MTLComputePipelineState> s_pso_grav_hydro = nil;
static id<MTLComputePipelineState> s_pso_godunov      = nil;
static id<MTLComputePipelineState> s_pso_uct_velocity = nil;
static id<MTLComputePipelineState> s_pso_uct_face_product = nil;
static id<MTLComputePipelineState> s_pso_uct_reuse = nil;
static id<MTLComputePipelineState> s_pso_build_nbor   = nil;
static id<MTLComputePipelineState> s_pso_scan_block        = nil;
static id<MTLComputePipelineState> s_pso_scan_fixup        = nil;
static id<MTLComputePipelineState> s_pso_reset_flag1       = nil;
static id<MTLComputePipelineState> s_pso_reset_flag2       = nil;
static id<MTLComputePipelineState> s_pso_init_flag         = nil;
static id<MTLComputePipelineState> s_pso_count_flag1       = nil;
static id<MTLComputePipelineState> s_pso_hydro_flag        = nil;
static id<MTLComputePipelineState> s_pso_poisson_flag      = nil;
static id<MTLComputePipelineState> s_pso_count_neighbors   = nil;
static id<MTLComputePipelineState> s_pso_flag_count        = nil;
static id<MTLComputePipelineState> s_pso_enforce_subgrid  = nil;
static id<MTLComputePipelineState> s_pso_enforce_rules     = nil;
static id<MTLComputePipelineState> s_pso_update_father     = nil;

/* PSOs for AMR refine/sort/cache kernels (refine.metal) */
static id<MTLComputePipelineState> s_pso_refine            = nil;
static id<MTLComputePipelineState> s_pso_gather_force      = nil;
static id<MTLComputePipelineState> s_pso_scatter_force     = nil;
static id<MTLComputePipelineState> s_pso_gather_phi        = nil;
static id<MTLComputePipelineState> s_pso_scatter_phi       = nil;
static id<MTLComputePipelineState> s_pso_derefine          = nil;
static id<MTLComputePipelineState> s_pso_free_hash         = nil;
static id<MTLComputePipelineState> s_pso_update_hash       = nil;
static id<MTLComputePipelineState> s_pso_insert_hash_all   = nil;
static id<MTLComputePipelineState> s_pso_init_swap_table   = nil;
static id<MTLComputePipelineState> s_pso_prefix_level      = nil;
static id<MTLComputePipelineState> s_pso_prefix_bit        = nil;
static id<MTLComputePipelineState> s_pso_local_swap        = nil;
static id<MTLComputePipelineState> s_pso_global_swap       = nil;
static id<MTLComputePipelineState> s_pso_gather_grid       = nil;
static id<MTLComputePipelineState> s_pso_scatter_grid      = nil;
static id<MTLComputePipelineState> s_pso_gather_flag       = nil;
static id<MTLComputePipelineState> s_pso_scatter_flag      = nil;
static id<MTLComputePipelineState> s_pso_gather_hydro      = nil;
static id<MTLComputePipelineState> s_pso_nbor_prefix       = nil;
static id<MTLComputePipelineState> s_pso_cache_swap        = nil;
static id<MTLComputePipelineState> s_pso_make_cache        = nil;
static id<MTLComputePipelineState> s_pso_insert_hash_cache = nil;
static id<MTLComputePipelineState> s_pso_upload            = nil;

/* PSOs for gravity / rho kernels (rho.metal) */
static id<MTLComputePipelineState> s_pso_reset_rho              = nil;
static id<MTLComputePipelineState> s_pso_multipole_leaf         = nil;
static id<MTLComputePipelineState> s_pso_multipole_upload       = nil;
static id<MTLComputePipelineState> s_pso_multipole_tot          = nil;
static id<MTLComputePipelineState> s_pso_deposit_rho            = nil;

/* PSOs for multigrid Poisson kernels (mg.metal) */
static id<MTLComputePipelineState> s_pso_save_phi_old           = nil;
static id<MTLComputePipelineState> s_pso_reset_phi_mg           = nil;
static id<MTLComputePipelineState> s_pso_make_initial_phi       = nil;
static id<MTLComputePipelineState> s_pso_reset_mask_mg          = nil;
static id<MTLComputePipelineState> s_pso_reset_rhs_mg           = nil;
static id<MTLComputePipelineState> s_pso_update_father_array    = nil;
static id<MTLComputePipelineState> s_pso_init_prefix_sum_mg     = nil;
static id<MTLComputePipelineState> s_pso_compute_father_swap    = nil;
static id<MTLComputePipelineState> s_pso_make_father_octs       = nil;
static id<MTLComputePipelineState> s_pso_restrict_mask_mg       = nil;
static id<MTLComputePipelineState> s_pso_volume_to_mask         = nil;
static id<MTLComputePipelineState> s_pso_cmp_residual           = nil;
static id<MTLComputePipelineState> s_pso_gauss_seidel           = nil;
static id<MTLComputePipelineState> s_pso_reset_phi_val          = nil;
static id<MTLComputePipelineState> s_pso_restrict_residual      = nil;
static id<MTLComputePipelineState> s_pso_interpolate_correct    = nil;
static id<MTLComputePipelineState> s_pso_residual_norm          = nil;
static id<MTLComputePipelineState> s_pso_rhs_norm               = nil;
static id<MTLComputePipelineState> s_pso_cmp_epot               = nil;
static id<MTLComputePipelineState> s_pso_cmp_rhomax             = nil;
static id<MTLComputePipelineState> s_pso_gradient_phi           = nil;
static id<MTLComputePipelineState> s_pso_update_nbor_array_mg   = nil;

/* PSOs for cooling */
static id<MTLComputePipelineState> s_pso_cooling = nil;

/* Cooling buffers */
static id<MTLBuffer> s_nH_tbl_d = nil;
static id<MTLBuffer> s_T2_tbl_d = nil;
static id<MTLBuffer> s_cool_d = nil;
static id<MTLBuffer> s_heat_d = nil;
static id<MTLBuffer> s_cool_com_d = nil;
static id<MTLBuffer> s_heat_com_d = nil;
static id<MTLBuffer> s_metal_d = nil;
static id<MTLBuffer> s_cool_prime_d = nil;
static id<MTLBuffer> s_heat_prime_d = nil;
static id<MTLBuffer> s_cool_com_prime_d = nil;
static id<MTLBuffer> s_heat_com_prime_d = nil;
static id<MTLBuffer> s_metal_prime_d = nil;
static BOOL s_table_uploaded = NO;
static int s_table_n1 = 0;
static int s_table_n2 = 0;
static float s_table_dlog_nH = 0.0f;
static float s_table_dlog_T2 = 0.0f;

/* Gravity AMR buffers */
static id<MTLBuffer> s_rho          = nil;   /* rho(8, ngridmax+ncachemax) float */
static id<MTLBuffer> s_nref         = nil;   /* nref(8, ngridmax+ncachemax) int  */
static id<MTLBuffer> s_phi          = nil;   /* phi(8, ngridmax+ncachemax) float */
static id<MTLBuffer> s_phi_old      = nil;   /* phi_old same shape               */
static id<MTLBuffer> s_f_grav       = nil;   /* f(8,3, ngridmax+ncachemax) float */
static id<MTLBuffer> s_multipole_buf = nil;  /* multipole(4,ngridmax) float       */
static id<MTLBuffer> s_scalar_buf   = nil;   /* single float for scalar reductions */

/* MG buffers */
static id<MTLBuffer> s_grid_mg      = nil;   /* oct_t[ngridmax_mg]           */
static id<MTLBuffer> s_phi_mg       = nil;   /* phi(8, ncachemax_mg) float    */
static id<MTLBuffer> s_f_mg         = nil;   /* f(8,3, ncachemax_mg) float    */
static id<MTLBuffer> s_nbor_mg      = nil;   /* nbor(27, ncachemax_mg) int    */
static id<MTLBuffer> s_father_mg    = nil;   /* father(ngridmax) int          */
static id<MTLBuffer> s_hash_key_mg  = nil;   /* hash_key_mg(hash_size_mg) long */
static id<MTLBuffer> s_hash_val_mg  = nil;   /* hash_val_mg(hash_size_mg) int  */

/* MG sizing kept for runtime use */
static int s_ngridmax_mg  = 0;
static int s_ncachemax_mg = 0;
static int s_hash_size_mg = 0;

static int s_nvar = 5;   /* set in mtl_alloc_amr; used in mtl_blit_unew_to_uold */
static int s_twotondim = 8;

/* PSOs for particle kernels (part.metal) */
static id<MTLComputePipelineState> s_pso_kick_drift_part      = nil;
static id<MTLComputePipelineState> s_pso_newdt_part           = nil;
static id<MTLComputePipelineState> s_pso_bucket_part          = nil;
static id<MTLComputePipelineState> s_pso_init_ps_hilbert_part = nil;
static id<MTLComputePipelineState> s_pso_write_swap_partition = nil;
static id<MTLComputePipelineState> s_pso_write_sortp_part     = nil;
static id<MTLComputePipelineState> s_pso_hkey_part            = nil;
static id<MTLComputePipelineState> s_pso_prefix_part_bit      = nil;
static id<MTLComputePipelineState> s_pso_gather_real_col      = nil;
static id<MTLComputePipelineState> s_pso_scatter_real_col     = nil;
static id<MTLComputePipelineState> s_pso_gather_real_1d       = nil;
static id<MTLComputePipelineState> s_pso_scatter_real_1d      = nil;
static id<MTLComputePipelineState> s_pso_gather_i8_1d         = nil;
static id<MTLComputePipelineState> s_pso_scatter_i8_1d        = nil;
static id<MTLComputePipelineState> s_pso_gather_i4_1d         = nil;
static id<MTLComputePipelineState> s_pso_scatter_i4_1d        = nil;
static id<MTLComputePipelineState> s_pso_cic_part_medium      = nil;
static id<MTLComputePipelineState> s_pso_multipole_q_part     = nil;

/* Particle buffers */
static id<MTLBuffer> s_xp          = nil;
static id<MTLBuffer> s_vp          = nil;
static id<MTLBuffer> s_mp          = nil;
static id<MTLBuffer> s_idp         = nil;
static id<MTLBuffer> s_levelp      = nil;
static id<MTLBuffer> s_sortp       = nil;
static id<MTLBuffer> s_xp_swap     = nil;
static id<MTLBuffer> s_isp_swap    = nil;
static id<MTLBuffer> s_idp_swap    = nil;
static id<MTLBuffer> s_prefix_sum_part = nil;
static id<MTLBuffer> s_multipole_q_part_buf = nil; /* float[4]: q[0..3] for particle monopole/dipole */
static id<MTLBuffer> s_part_scalar_buf = nil;
static id<MTLBuffer> s_partial_sums_part  = nil;
static id<MTLBuffer> s_partial_sums_part2 = nil;
static id<MTLBuffer> s_partial_sums_part3 = nil;
static id<MTLBuffer> s_partial_sums_part4 = nil;
static int s_npartmax = 0;

/* ----------------------------------------------------------------------- */
static id<MTLComputePipelineState> make_pso(NSString *name)
{
    NSError *error = nil;
    id<MTLFunction> fn = [s_library newFunctionWithName:name];
    if (!fn) {
        fprintf(stderr, "[metal] kernel not found: %s\n", [name UTF8String]);
        exit(1);
    }
    id<MTLComputePipelineState> pso =
        [s_device newComputePipelineStateWithFunction:fn error:&error];
    if (!pso) {
        fprintf(stderr, "[metal] pipeline error for %s: %s\n",
                [name UTF8String], [[error localizedDescription] UTF8String]);
        exit(1);
    }
    return pso;
}

/* -----------------------------------------------------------------------
 * mtl_init — create device, queue, load metallib, build pipeline states.
 * Mirrors mdl_initialize / cudaSetDevice in the CUDA path.
 * The metallib is expected next to the running executable (bin/).
 * ----------------------------------------------------------------------- */
extern "C" void mtl_init(void)
{
    s_device = MTLCreateSystemDefaultDevice();
    if (!s_device) {
        fprintf(stderr, "[metal] no Metal device found\n");
        exit(1);
    }
    s_queue = [s_device newCommandQueue];

    NSString *exec_path = [[[NSProcessInfo processInfo] arguments] firstObject];
    NSString *exec_dir  = [exec_path stringByDeletingLastPathComponent];
    NSString *lib_path  = [exec_dir stringByAppendingPathComponent:
                           @"ramses_kernels.metallib"];
    NSURL    *lib_url   = [NSURL fileURLWithPath:lib_path];

    NSError *error = nil;
    s_library = [s_device newLibraryWithURL:lib_url error:&error];
    if (!s_library) {
        fprintf(stderr, "[metal] cannot load %s: %s\n",
                [lib_path UTF8String],
                [[error localizedDescription] UTF8String]);
        exit(1);
    }

    s_pso_set_unew     = make_pso(@"set_unew_kernel");
    s_pso_set_uold     = make_pso(@"set_uold_kernel");
    s_pso_cmpdt        = make_pso(@"cmpdt_kernel");
    s_pso_godunov      = make_pso(@"hydro_integrator_kernel");
#ifdef MHD
    s_pso_uct_velocity = make_pso(@"uct_velocity_kernel");
    s_pso_uct_face_product = make_pso(@"mhd_uct_face_product_kernel");
    s_pso_uct_reuse = make_pso(@"hydro_integrator_uct_reuse_kernel");
#endif
    s_pso_sync_hydro   = make_pso(@"sync_hydro_kernel");
    s_pso_grav_hydro   = make_pso(@"grav_hydro_kernel");
    s_pso_build_nbor   = make_pso(@"build_nbor_kernel");
    s_pso_scan_block        = make_pso(@"scan_block_kernel");
    s_pso_scan_fixup        = make_pso(@"scan_fixup_kernel");
    s_pso_reset_flag1       = make_pso(@"reset_flag1_kernel");
    s_pso_reset_flag2       = make_pso(@"reset_flag2_kernel");
    s_pso_init_flag         = make_pso(@"init_flag_kernel");
    s_pso_count_flag1       = make_pso(@"count_flag1_kernel");
    s_pso_hydro_flag        = make_pso(@"hydro_flag_kernel");
    s_pso_poisson_flag      = make_pso(@"poisson_flag_kernel");
    s_pso_count_neighbors   = make_pso(@"count_neighbors_kernel");
    s_pso_flag_count        = make_pso(@"flag_count_kernel");
    s_pso_enforce_subgrid   = make_pso(@"enforce_subgrid_kernel");
    s_pso_enforce_rules     = make_pso(@"enforce_rules_kernel");
    s_pso_update_father     = make_pso(@"update_father_kernel");

    s_pso_refine            = make_pso(@"refine_kernel");
    s_pso_derefine          = make_pso(@"derefine_kernel");
    s_pso_free_hash         = make_pso(@"free_hash_kernel");
    s_pso_update_hash       = make_pso(@"update_hash_kernel");
    s_pso_insert_hash_all   = make_pso(@"insert_hash_all_kernel");
    s_pso_init_swap_table   = make_pso(@"init_global_swap_table_kernel");
    s_pso_prefix_level      = make_pso(@"init_prefix_sum_level_kernel");
    s_pso_prefix_bit        = make_pso(@"init_prefix_sum_bit_kernel");
    s_pso_local_swap        = make_pso(@"compute_local_swap_table_kernel");
    s_pso_global_swap       = make_pso(@"update_global_swap_table_kernel");
    s_pso_gather_grid       = make_pso(@"sort_gather_grid_kernel");
    s_pso_gather_force      = make_pso(@"sort_gather_force_kernel");
    s_pso_scatter_force     = make_pso(@"sort_scatter_force_kernel");
    s_pso_gather_phi        = make_pso(@"sort_gather_phi_kernel");
    s_pso_scatter_phi       = make_pso(@"sort_scatter_phi_kernel");
    s_pso_scatter_grid      = make_pso(@"sort_scatter_grid_kernel");
    s_pso_gather_flag       = make_pso(@"sort_gather_flag_kernel");
    s_pso_scatter_flag      = make_pso(@"sort_scatter_flag_kernel");
    s_pso_gather_hydro      = make_pso(@"sort_gather_hydro_kernel");
    s_pso_nbor_prefix       = make_pso(@"update_nbor_prefix_kernel");
    s_pso_cache_swap        = make_pso(@"compute_cache_swap_table_kernel");
    s_pso_make_cache        = make_pso(@"make_cache_octs_kernel");
    s_pso_insert_hash_cache = make_pso(@"insert_hash_cache_kernel");
    s_pso_upload            = make_pso(@"upload_kernel");

    /* gravity / rho kernels */
    s_pso_reset_rho              = make_pso(@"reset_rho_kernel");
    s_pso_multipole_leaf         = make_pso(@"multipole_leaf_kernel");
    s_pso_multipole_upload       = make_pso(@"multipole_upload_kernel");
    s_pso_multipole_tot          = make_pso(@"multipole_tot_kernel");
    s_pso_deposit_rho            = make_pso(@"deposit_rho_kernel");

    /* multigrid Poisson kernels */
    s_pso_save_phi_old           = make_pso(@"save_phi_old_kernel");
    s_pso_reset_phi_mg           = make_pso(@"reset_phi_kernel_mg");
    s_pso_make_initial_phi       = make_pso(@"make_initial_phi_kernel");
    s_pso_reset_mask_mg          = make_pso(@"reset_mask_kernel_mg");
    s_pso_reset_rhs_mg           = make_pso(@"reset_rhs_kernel_mg");
    s_pso_update_father_array    = make_pso(@"update_father_array_kernel");
    s_pso_init_prefix_sum_mg     = make_pso(@"init_prefix_sum_mg_kernel");
    s_pso_compute_father_swap    = make_pso(@"compute_father_swap_kernel");
    s_pso_make_father_octs       = make_pso(@"make_father_octs_kernel");
    s_pso_restrict_mask_mg       = make_pso(@"restrict_mask_kernel_mg");
    s_pso_volume_to_mask         = make_pso(@"volume_to_mask_kernel");
    s_pso_cmp_residual           = make_pso(@"cmp_residual_kernel");
    s_pso_gauss_seidel           = make_pso(@"gauss_seidel_kernel");
    s_pso_reset_phi_val          = make_pso(@"reset_phi_val_kernel");
    s_pso_restrict_residual      = make_pso(@"restrict_residual_kernel");
    s_pso_interpolate_correct    = make_pso(@"interpolate_correct_kernel");
    s_pso_residual_norm          = make_pso(@"residual_norm_kernel");
    s_pso_rhs_norm               = make_pso(@"rhs_norm_kernel");
    s_pso_cmp_epot               = make_pso(@"cmp_epot_kernel");
    s_pso_cmp_rhomax             = make_pso(@"cmp_rhomax_kernel");
    s_pso_gradient_phi           = make_pso(@"gradient_phi_kernel");
    s_pso_update_nbor_array_mg   = make_pso(@"update_nbor_array_mg_kernel");
    s_pso_cooling                = make_pso(@"cooling_kernel");

    /* particle kernels */
    s_pso_kick_drift_part      = make_pso(@"kick_drift_part_kernel");
    s_pso_newdt_part           = make_pso(@"newdt_part_kernel");
    s_pso_bucket_part          = make_pso(@"bucket_part_kernel");
    s_pso_init_ps_hilbert_part = make_pso(@"init_prefix_sum_part_hilbert_fine");
    s_pso_write_swap_partition = make_pso(@"write_swap_global_hilbert_partition");
    s_pso_write_sortp_part     = make_pso(@"write_sortp_part");
    s_pso_hkey_part            = make_pso(@"compute_hkey_part_kernel");
    s_pso_prefix_part_bit      = make_pso(@"init_prefix_sum_part_bit");
    s_pso_gather_real_col      = make_pso(@"sort_gather_part_real_col");
    s_pso_scatter_real_col     = make_pso(@"sort_scatter_part_real_col");
    s_pso_gather_real_1d       = make_pso(@"sort_gather_part_real_1d");
    s_pso_scatter_real_1d      = make_pso(@"sort_scatter_part_real_1d");
    s_pso_gather_i8_1d         = make_pso(@"sort_gather_part_i8_1d");
    s_pso_scatter_i8_1d        = make_pso(@"sort_scatter_part_i8_1d");
    s_pso_gather_i4_1d         = make_pso(@"sort_gather_part_i4_1d");
    s_pso_scatter_i4_1d        = make_pso(@"sort_scatter_part_i4_1d");
    s_pso_cic_part_medium      = make_pso(@"cic_part_medium_kernel");
    s_pso_multipole_q_part     = make_pso(@"multipole_q_part_kernel");

    fprintf(stdout, " Launching METAL.\n");
    s_sort_fence = [s_device newFence];
    s_count_buf  = [s_device newBufferWithLength:sizeof(int)
                                         options:MTLResourceStorageModeShared];
    *(int *)s_count_buf.contents = 0;

    fprintf(stdout, " Device name:       %s\n", [[s_device name] UTF8String]);
    fprintf(stdout, " Unified memory:    %s\n",
            s_device.hasUnifiedMemory ? "yes" : "no");
    fprintf(stdout, " Max working set:   %llu MB\n",
            (unsigned long long)s_device.recommendedMaxWorkingSetSize / (1024*1024));
    fprintf(stdout, " Max buffer length: %llu MB\n",
            (unsigned long long)s_device.maxBufferLength / (1024*1024));
    fprintf(stdout, " Max threads/tg:    %lu\n",
            (unsigned long)s_device.maxThreadsPerThreadgroup.width);
}

/* -----------------------------------------------------------------------
 * mtl_alloc_amr — allocate Metal-owned buffers for uold, unew, grid.
 * MTLResourceStorageModeShared: buffer lives in CPU/GPU shared DRAM.
 * Data is copied from the Fortran arrays in mtl_set_grid_device.
 * ----------------------------------------------------------------------- */
/* ncachemax is included so grid, nbor, and the hydro arrays cover the full
 * [1..ngridmax+ncachemax] range required by AMR ghost-zone caching.
 * For the unigrid PoC (ncachemax>0 but cache never populated) the extra
 * allocation is harmless — behaviour identical to the old 4-argument form. */
extern "C" void mtl_alloc_amr(int ngridmax, int ncachemax,
                               int nvar, int twotondim, int hash_size)
{
    s_nvar = nvar;
    s_twotondim = twotondim;
    s_ngridmax = ngridmax;
    s_hash_size = hash_size;
    int ntotal = ngridmax + ncachemax;
    NSUInteger u_bytes        = (NSUInteger)ntotal    * nvar * twotondim * sizeof(float);
    NSUInteger grid_bytes     = (NSUInteger)ntotal    * sizeof(oct_t);
#ifdef MHD
    NSUInteger nbor_bytes     = (NSUInteger)ntotal * (NSUBGRID + 2) * (NSUBGRID + 2) * (NSUBGRID + 2) * sizeof(int);
#else
    NSUInteger nbor_bytes     = (NSUInteger)ntotal * 27 * sizeof(int);
#endif
    NSUInteger hash_key_bytes = (NSUInteger)hash_size * sizeof(long);
    NSUInteger hash_val_bytes = (NSUInteger)hash_size * sizeof(int);

    s_uold     = [s_device newBufferWithLength:u_bytes
                                       options:MTLResourceStorageModeShared];
    s_unew     = [s_device newBufferWithLength:u_bytes
                                       options:MTLResourceStorageModeShared];
#ifdef MHD
    NSUInteger b_bytes = (NSUInteger)ntotal * 6 * twotondim * sizeof(float);
    s_bold = [s_device newBufferWithLength:b_bytes options:MTLResourceStorageModeShared];
    s_bnew = [s_device newBufferWithLength:b_bytes options:MTLResourceStorageModeShared];
#endif
    s_grid     = [s_device newBufferWithLength:grid_bytes
                                       options:MTLResourceStorageModeShared];
    s_nbor     = [s_device newBufferWithLength:nbor_bytes
                                       options:MTLResourceStorageModeShared];
    s_hash_key = [s_device newBufferWithLength:hash_key_bytes
                                       options:MTLResourceStorageModeShared];
    s_hash_val = [s_device newBufferWithLength:hash_val_bytes
                                       options:MTLResourceStorageModeShared];
}

/* -----------------------------------------------------------------------
 * mtl_set_grid_device — copy host arrays into Metal buffers (H->D).
 * Mirrors the cudaMemcpy calls in r_set_grid_device (gpu_manager.cuf).
 * On Apple Silicon the memcpy stays within DRAM (no PCIe), but is still
 * needed because the Metal buffer and the Fortran array are distinct
 * allocations at different addresses.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_set_grid_device(void *uold_ptr, void *unew_ptr,
                                    void *bold_ptr,
                                    void *grid_ptr,
                                    int ngridmax, int nvar, int twotondim)
{
    size_t u_bytes    = (size_t)ngridmax * nvar * twotondim * sizeof(float);
    size_t grid_bytes = (size_t)ngridmax * sizeof(oct_t);
    if (uold_ptr && s_uold) {
        memcpy(s_uold.contents, uold_ptr, u_bytes);
    }
    if (unew_ptr && s_unew) {
        memcpy(s_unew.contents, unew_ptr, u_bytes);
    }
    size_t b_bytes = (size_t)ngridmax * 6 * twotondim * sizeof(float);
    if (bold_ptr && s_bold) memcpy(s_bold.contents, bold_ptr, b_bytes);
    if (bold_ptr && s_bnew) memcpy(s_bnew.contents, bold_ptr, b_bytes);
    memcpy(s_grid.contents, grid_ptr, grid_bytes);
}

/* -----------------------------------------------------------------------
 * mtl_upload_flag1 — copy host flag1(8,ngridmax) to device s_flag1.
 * Called from r_set_grid_device when nlevelmax > levelmin so that
 * derefine_kernel reads the correct refinement flags, not stale zeros.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_upload_flag1(void *flag1_host, int ngridmax)
{
    if (flag1_host && s_flag1) {
        size_t nbytes = (size_t)8 * ngridmax * sizeof(int);
        memcpy(s_flag1.contents, flag1_host, nbytes);
    }
}

/* -----------------------------------------------------------------------
 * mtl_transfer_grid_host — copy Metal uold buffer back to host (D->H).
 * Called before each output dump so m%uold reflects the GPU result.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_transfer_grid_host(void *uold_ptr,
                                       void *bold_ptr,
                                       int ngridmax, int nvar, int twotondim)
{
    if (uold_ptr && s_uold) {
        size_t u_bytes = (size_t)ngridmax * nvar * twotondim * sizeof(float);
        memcpy(uold_ptr, s_uold.contents, u_bytes);
    }
    if (bold_ptr && s_bold) {
        size_t b_bytes = (size_t)ngridmax * 6 * twotondim * sizeof(float);
        memcpy(bold_ptr, s_bold.contents, b_bytes);
    }
}

/* -----------------------------------------------------------------------
 * mtl_transfer_grid_struct_host — copy Metal s_grid buffer back to host.
 * Required for AMR runs (levelmin < levelmax): metal_refine reorders octs
 * via Hilbert sort + scatter, updating s_grid on the device.  Without this
 * readback, output_amr reads stale host ckey/refined values and amr2map
 * produces a garbled level/density map.
 * Only the first ngridmax slots are written (cache octs start at ngridmax+1
 * and are never referenced by the output routines).
 * ----------------------------------------------------------------------- */
extern "C" void mtl_transfer_grid_struct_host(void *grid_ptr, int ngridmax)
{
    size_t grid_bytes = (size_t)ngridmax * sizeof(oct_t);
    memcpy(grid_ptr, s_grid.contents, grid_bytes);
}

/* -----------------------------------------------------------------------
 * mtl_device_sync — block until all previously submitted Metal work completes.
 * An empty command buffer committed to the queue is sufficient: Metal
 * serialises command buffers in submission order, so waiting on this empty
 * one guarantees all prior dispatches have finished.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_device_sync(void)
{
    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];
    [cmd commit];
    [cmd waitUntilCompleted];
}

extern "C" int mtl_nsubgrid(void)
{
    return NSUBGRID;
}


/* -----------------------------------------------------------------------
 * mtl_build_nbor — build the device nbor array from the already-populated
 * hash table by dispatching build_nbor_kernel.
 * Thread layout: 128 threads/threadgroup.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_build_nbor(int head_idx, int num_subgrids,
                                int hash_size,
                                int ckey_max_l,    long key_off_l,
                                int *box_ckey_min, int *box_ckey_max,
                                int *periodic)
{
    NSUInteger tg128  = 128;
    MTLSize tg_size   = {tg128, 1, 1};
    MTLSize grid_size = {((NSUInteger)num_subgrids + tg128 - 1) / tg128, 1, 1};

    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_build_nbor];
    [enc setBuffer:s_grid     offset:0 atIndex:0];
    [enc setBuffer:s_nbor     offset:0 atIndex:1];
    [enc setBuffer:s_hash_key offset:0 atIndex:2];
    [enc setBuffer:s_hash_val offset:0 atIndex:3];
    [enc setBytes:&hash_size      length:sizeof(int)      atIndex:4];
    [enc setBytes:&ckey_max_l     length:sizeof(int)      atIndex:5];
    [enc setBytes:&key_off_l      length:sizeof(long)     atIndex:6];
    [enc setBytes:box_ckey_min    length:3 * sizeof(int)  atIndex:7];
    [enc setBytes:box_ckey_max    length:3 * sizeof(int)  atIndex:8];
    [enc setBytes:periodic        length:3 * sizeof(int)  atIndex:9];
    [enc setBytes:&head_idx       length:sizeof(int)      atIndex:10];
    [enc setBytes:&num_subgrids   length:sizeof(int)      atIndex:11];
    [enc dispatchThreadgroups:grid_size threadsPerThreadgroup:tg_size];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_set_unew — dispatch set_unew_kernel: unew = uold for octs at ilevel.
 * Thread layout mirrors CUDA: dim3(8, 16, 1) per threadgroup,
 * ceil(num_octs / 16) threadgroups.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_set_unew(int head_idx, int num_octs)
{
    MTLSize tg_size   = {8, 16, 1};
    MTLSize grid_size = {((NSUInteger)num_octs + 15) / 16, 1, 1};

    id<MTLCommandBuffer>        cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_set_unew];
    [enc setBuffer:s_uold   offset:0 atIndex:0];
    [enc setBuffer:s_unew   offset:0 atIndex:1];
#ifdef MHD
    [enc setBuffer:s_bold   offset:0 atIndex:2];
    [enc setBuffer:s_bnew   offset:0 atIndex:3];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:4];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:5];
#else
    [enc setBytes:&head_idx length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:3];
#endif
    [enc dispatchThreadgroups:grid_size threadsPerThreadgroup:tg_size];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_set_uold — dispatch set_uold_kernel: uold = unew for octs at ilevel.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_set_uold(int head_idx, int num_octs)
{
    MTLSize tg_size   = {8, 16, 1};
    MTLSize grid_size = {((NSUInteger)num_octs + 15) / 16, 1, 1};

    id<MTLCommandBuffer>        cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_set_uold];
    [enc setBuffer:s_uold   offset:0 atIndex:0];
    [enc setBuffer:s_unew   offset:0 atIndex:1];
#ifdef MHD
    [enc setBuffer:s_bold   offset:0 atIndex:2];
    [enc setBuffer:s_bnew   offset:0 atIndex:3];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:4];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:5];
#else
    [enc setBytes:&head_idx length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:3];
#endif
    [enc dispatchThreadgroups:grid_size threadsPerThreadgroup:tg_size];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_cmpdt — dispatch cmpdt_kernel and read back results.
 *
 * data_buf layout: atomic_uint[5] reinterpreted as float[5] on readback.
 *   [0..3] fp32 accumulated via CAS atomic_add_float
 *   [4]    fp32 min via uint bit-cast atomic_min_float_bits
 *
 * 1024 threads/threadgroup; SIMD reduction inside the kernel collapses to
 * one atomic write per threadgroup.  dispatchThreadgroups with
 * ceil(num_octs*8 / 1024) threadgroups.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_cmpdt(int head_idx, int num_octs,
                           float dx, float gamma, float smallr, float smallc2,
                           float courant_factor,
#ifdef MHD
                           int induction,
#endif
                           float *constant_gravity,
                           float *mass, float *ekin, float *eint, float *emag,
                           float *dt)
{
    float dt_init = courant_factor * dx / sqrtf(smallc2);
    /* data_buf layout: atomic_uint[5] = {mass, ekin, eint, emag, dt}
     * [0..3] initialised to 0 (float bit-pattern); [4] initialised to dt_init. */
    uint32_t h_data[5] = {0, 0, 0, 0, 0};
    memcpy(&h_data[4], &dt_init, sizeof(float));
    id<MTLBuffer> data_buf =
        [s_device newBufferWithBytes:h_data
                              length:5 * sizeof(uint32_t)
                             options:MTLResourceStorageModeShared];

    float cg[3] = {constant_gravity[0], constant_gravity[1], constant_gravity[2]};

    NSUInteger total_cells = (NSUInteger)num_octs * 8;   /* 8 = twotondim for NDIM=3 */
    MTLSize tg_size   = {1024, 1, 1};
    MTLSize grid_size = {(total_cells + 1023) / 1024, 1, 1};

    id<MTLCommandBuffer>        cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_cmpdt];
    [enc setBuffer:s_grid    offset:0 atIndex:0];
    [enc setBuffer:s_uold    offset:0 atIndex:1];
#ifdef MHD
    [enc setBuffer:s_bold    offset:0 atIndex:2];
    [enc setBuffer:data_buf  offset:0 atIndex:3];
    [enc setBytes:&head_idx       length:sizeof(int)       atIndex:4];
    [enc setBytes:&num_octs       length:sizeof(int)       atIndex:5];
    [enc setBytes:&dx             length:sizeof(float)     atIndex:6];
    [enc setBytes:&gamma          length:sizeof(float)     atIndex:7];
    [enc setBytes:&smallr         length:sizeof(float)     atIndex:8];
    [enc setBytes:&smallc2        length:sizeof(float)     atIndex:9];
    [enc setBytes:&courant_factor length:sizeof(float)     atIndex:10];
    [enc setBytes:&induction      length:sizeof(int)       atIndex:11];
    [enc setBytes:cg              length:3 * sizeof(float) atIndex:12];
    [enc setBuffer:s_f_grav       offset:0                 atIndex:13];
#else
    [enc setBuffer:data_buf  offset:0 atIndex:2];
    [enc setBytes:&head_idx       length:sizeof(int)       atIndex:3];
    [enc setBytes:&num_octs       length:sizeof(int)       atIndex:4];
    [enc setBytes:&dx             length:sizeof(float)     atIndex:5];
    [enc setBytes:&gamma          length:sizeof(float)     atIndex:6];
    [enc setBytes:&smallr         length:sizeof(float)     atIndex:7];
    [enc setBytes:&smallc2        length:sizeof(float)     atIndex:8];
    [enc setBytes:&courant_factor length:sizeof(float)     atIndex:9];
    [enc setBytes:cg              length:3 * sizeof(float) atIndex:10];
    [enc setBuffer:s_f_grav       offset:0                 atIndex:11];
#endif

    [enc dispatchThreadgroups:grid_size threadsPerThreadgroup:tg_size];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    float *result = (float *)data_buf.contents;
    *mass = result[0];
    *ekin = result[1];
    *eint = result[2];
    *emag = result[3];
    *dt   = result[4];
}

#ifdef MHD
static bool mtl_uct_face_reuse_enabled(int ilevel, int levelmin, int levelmax)
{
    // Relative product indexing is valid only when this dispatch covers the complete periodic level.
    if (!s_periodic_dev) return false;
    int *periodic = (int *)s_periodic_dev.contents;
    return ilevel == levelmin && ilevel == levelmax &&
           periodic[0] != 0 && periodic[1] != 0 && periodic[2] != 0;
}
#endif

/* -----------------------------------------------------------------------
 * mtl_godunov — dispatch the hydro or MHD integrator.
 * One threadgroup handles each subgrid.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_godunov(int head_idx, int num_subgrids, int ngridmax,
                             int ilevel, int levelmin, int levelmax,
                             float gamma, float smallr, float smallc2,
                             float dt, float dx, int slope,
#ifdef MHD
                             int slope_mag, int riemann, int riemann2d,
                             float switch_llf_dmin, float switch_llf_pmin,
                             int induction, float etamag,
#else
                             int riemann,
#endif
                             float *constant_gravity)
{
    float cg[3] = {constant_gravity[0], constant_gravity[1], constant_gravity[2]};

#ifdef MHD
    MTLSize tg_size = {256, 1, 1};
#else
    MTLSize tg_size   = {64, 1, 1};
#endif
    MTLSize grid_size = {(NSUInteger)num_subgrids, 1, 1};

    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc;
#ifdef MHD
    bool use_uct_face_reuse = !s_uct_face_reuse_disabled && riemann == MR_SOLVER_UCT_HLLD &&
                              mtl_uct_face_reuse_enabled(ilevel, levelmin, levelmax);
    if (riemann == MR_SOLVER_UCT_HLLD) {
        if (use_uct_face_reuse) {
            NSUInteger product_bytes = (NSUInteger)num_subgrids * MR_UCT_PRODUCT_FIELDS *
                                       MR_UCT_PRODUCT_FACES * sizeof(float);
            if (product_bytes > s_device.maxBufferLength) {
                fprintf(stderr, "[metal] UCT face-product buffer exceeds maxBufferLength (%llu bytes); using legacy path\n",
                        (unsigned long long)product_bytes);
                s_uct_face_reuse_disabled = true;
                use_uct_face_reuse = false;
            }
            if (use_uct_face_reuse && (!s_uct_face_product || s_uct_face_product.length < product_bytes)) {
                s_uct_face_product = nil;
                s_uct_face_product = [s_device newBufferWithLength:product_bytes
                                                           options:MTLResourceStorageModePrivate];
                if (!s_uct_face_product) {
                    fprintf(stderr, "[metal] cannot allocate UCT face-product buffer (%llu bytes); using legacy path\n",
                            (unsigned long long)product_bytes);
                    s_uct_face_reuse_disabled = true;
                    use_uct_face_reuse = false;
                }
            }
        }
        if (!use_uct_face_reuse) {
            NSUInteger velocity_bytes = (NSUInteger)ngridmax * 3 * s_twotondim * 2 * sizeof(float);
            if (velocity_bytes > s_device.maxBufferLength) {
                fprintf(stderr, "[metal] UCT velocity buffer exceeds maxBufferLength (%llu bytes)\n",
                        (unsigned long long)velocity_bytes);
                exit(1);
            }
            if (!s_uct_velocity || s_uct_velocity.length < velocity_bytes) {
                s_uct_velocity = [s_device newBufferWithLength:velocity_bytes
                                                       options:MTLResourceStorageModePrivate];
                if (!s_uct_velocity) {
                    fprintf(stderr, "[metal] cannot allocate UCT velocity buffer (%llu bytes)\n",
                            (unsigned long long)velocity_bytes);
                    exit(1);
                }
            }
        }
        enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:use_uct_face_reuse ? s_pso_uct_face_product : s_pso_uct_velocity];
        [enc setBuffer:s_grid         offset:0 atIndex:0];
        [enc setBuffer:s_uold         offset:0 atIndex:1];
        [enc setBuffer:s_bold         offset:0 atIndex:2];
        [enc setBuffer:s_nbor         offset:0 atIndex:3];
        [enc setBytes:&head_idx       length:sizeof(int)       atIndex:4];
        [enc setBytes:&num_subgrids   length:sizeof(int)       atIndex:5];
        [enc setBytes:&gamma          length:sizeof(float)     atIndex:6];
        [enc setBytes:&smallr         length:sizeof(float)     atIndex:7];
        [enc setBytes:&smallc2        length:sizeof(float)     atIndex:8];
        [enc setBytes:&dt             length:sizeof(float)     atIndex:9];
        [enc setBytes:&dx             length:sizeof(float)     atIndex:10];
        [enc setBytes:&slope          length:sizeof(int)       atIndex:11];
        [enc setBytes:&slope_mag      length:sizeof(int)       atIndex:12];
        [enc setBytes:&switch_llf_dmin length:sizeof(float)    atIndex:13];
        [enc setBytes:&switch_llf_pmin length:sizeof(float)    atIndex:14];
        [enc setBytes:&induction      length:sizeof(int)       atIndex:15];
        [enc setBytes:cg              length:3 * sizeof(float) atIndex:16];
        [enc setBuffer:s_f_grav       offset:0                 atIndex:17];
        [enc setBuffer:use_uct_face_reuse ? s_uct_face_product : s_uct_velocity offset:0 atIndex:18];
        [enc dispatchThreadgroups:grid_size threadsPerThreadgroup:tg_size];
        [enc endEncoding];
        if (use_uct_face_reuse) {
            enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:s_pso_uct_reuse];
            [enc setBuffer:s_unew             offset:0 atIndex:0];
            [enc setBuffer:s_bold             offset:0 atIndex:1];
            [enc setBuffer:s_bnew             offset:0 atIndex:2];
            [enc setBuffer:s_nbor             offset:0 atIndex:3];
            [enc setBytes:&head_idx           length:sizeof(int)   atIndex:4];
            [enc setBytes:&num_subgrids       length:sizeof(int)   atIndex:5];
            [enc setBytes:&dt                 length:sizeof(float) atIndex:6];
            [enc setBytes:&dx                 length:sizeof(float) atIndex:7];
            [enc setBytes:&slope_mag          length:sizeof(int)   atIndex:8];
            [enc setBytes:&etamag             length:sizeof(float) atIndex:9];
            [enc setBuffer:s_uct_face_product offset:0             atIndex:10];
            [enc dispatchThreadgroups:grid_size threadsPerThreadgroup:tg_size];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
            return;
        }
    }
#endif
    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_godunov];
    [enc setBuffer:s_grid  offset:0 atIndex:0];
    [enc setBuffer:s_uold  offset:0 atIndex:1];
    [enc setBuffer:s_unew  offset:0 atIndex:2];
#ifdef MHD
    [enc setBuffer:s_bold  offset:0 atIndex:3];
    [enc setBuffer:s_bnew  offset:0 atIndex:4];
    [enc setBuffer:s_father offset:0 atIndex:5];
    [enc setBuffer:s_nbor  offset:0 atIndex:6];
    [enc setBytes:&head_idx     length:sizeof(int)       atIndex:7];
    [enc setBytes:&num_subgrids length:sizeof(int)       atIndex:8];
    [enc setBytes:&ngridmax    length:sizeof(int)       atIndex:9];
    [enc setBytes:&ilevel       length:sizeof(int)       atIndex:10];
    [enc setBytes:&levelmin     length:sizeof(int)       atIndex:11];
    [enc setBytes:&levelmax     length:sizeof(int)       atIndex:12];
    [enc setBytes:&gamma        length:sizeof(float)     atIndex:13];
    [enc setBytes:&smallr       length:sizeof(float)     atIndex:14];
    [enc setBytes:&smallc2      length:sizeof(float)     atIndex:15];
    [enc setBytes:&dt           length:sizeof(float)     atIndex:16];
    [enc setBytes:&dx           length:sizeof(float)     atIndex:17];
    [enc setBytes:&slope        length:sizeof(int)       atIndex:18];
    [enc setBytes:&slope_mag    length:sizeof(int)       atIndex:19];
    [enc setBytes:&riemann      length:sizeof(int)       atIndex:20];
    [enc setBytes:&riemann2d    length:sizeof(int)       atIndex:21];
    [enc setBytes:&switch_llf_dmin length:sizeof(float)  atIndex:22];
    [enc setBytes:&switch_llf_pmin length:sizeof(float)  atIndex:23];
    [enc setBytes:&induction    length:sizeof(int)       atIndex:24];
    [enc setBytes:&etamag       length:sizeof(float)     atIndex:25];
    [enc setBytes:cg            length:3 * sizeof(float) atIndex:26];
    [enc setBuffer:s_f_grav     offset:0                 atIndex:27];
    [enc setBuffer:s_uct_velocity ? s_uct_velocity : s_uold offset:0 atIndex:28];
#else
    [enc setBuffer:s_nbor  offset:0 atIndex:3];
    [enc setBytes:&head_idx     length:sizeof(int)       atIndex:4];
    [enc setBytes:&num_subgrids length:sizeof(int)       atIndex:5];
    [enc setBytes:&ngridmax     length:sizeof(int)       atIndex:6];
    [enc setBytes:&ilevel       length:sizeof(int)       atIndex:7];
    [enc setBytes:&levelmin     length:sizeof(int)       atIndex:8];
    [enc setBytes:&levelmax     length:sizeof(int)       atIndex:9];
    [enc setBytes:&gamma        length:sizeof(float)     atIndex:10];
    [enc setBytes:&smallr       length:sizeof(float)     atIndex:11];
    [enc setBytes:&smallc2      length:sizeof(float)     atIndex:12];
    [enc setBytes:&dt           length:sizeof(float)     atIndex:13];
    [enc setBytes:&dx           length:sizeof(float)     atIndex:14];
    [enc setBytes:&slope        length:sizeof(int)       atIndex:15];
    [enc setBytes:&riemann      length:sizeof(int)       atIndex:16];
    [enc setBytes:cg            length:3 * sizeof(float) atIndex:17];
    [enc setBuffer:s_father     offset:0               atIndex:18];
    [enc setBuffer:s_f_grav     offset:0               atIndex:19];
#endif
    [enc dispatchThreadgroups:grid_size threadsPerThreadgroup:tg_size];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

extern "C" void mtl_sync_hydro(int head_idx, int num_octs,
                               float gamma, float smallr, float smallc2,
                               float dt, float *constant_gravity)
{
    if (num_octs <= 0) return;
    float cg[3] = {constant_gravity[0], constant_gravity[1], constant_gravity[2]};
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_sync_hydro];
    [enc setBuffer:s_uold     offset:0 atIndex:0];
    [enc setBuffer:s_f_grav   offset:0 atIndex:1];
    [enc setBytes:cg          length:3 * sizeof(float) atIndex:2];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:4];
    [enc setBytes:&gamma      length:sizeof(float) atIndex:5];
    [enc setBytes:&smallr     length:sizeof(float) atIndex:6];
    [enc setBytes:&smallc2    length:sizeof(float) atIndex:7];
    [enc setBytes:&dt         length:sizeof(float) atIndex:8];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}


extern "C" void mtl_grav_hydro(int head_idx, int num_octs,
                               float gamma, float smallr, float smallc2,
                               float dt, float *constant_gravity)
{
    if (num_octs <= 0) return;
    float cg[3] = {constant_gravity[0], constant_gravity[1], constant_gravity[2]};
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_grav_hydro];
    [enc setBuffer:s_uold     offset:0 atIndex:0];
    [enc setBuffer:s_unew     offset:0 atIndex:1];
    [enc setBuffer:s_f_grav   offset:0 atIndex:2];
    [enc setBytes:cg          length:3 * sizeof(float) atIndex:3];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:4];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:5];
    [enc setBytes:&gamma      length:sizeof(float) atIndex:6];
    [enc setBytes:&smallr     length:sizeof(float) atIndex:7];
    [enc setBytes:&smallc2    length:sizeof(float) atIndex:8];
    [enc setBytes:&dt         length:sizeof(float) atIndex:9];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_upload — dispatch upload_kernel (restriction: fine → coarse level).
 * One thread per fine oct (128 threads/threadgroup).
 * ----------------------------------------------------------------------- */
extern "C" void mtl_upload(int head_idx, int num_octs,
                            int internal_energy,
                            float gamma, float smallr, float smallc2)
{
    if (num_octs <= 0) return;
    NSUInteger tg_size = 128;
    NSUInteger num_tg  = ((NSUInteger)num_octs + tg_size - 1) / tg_size;

    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_upload];
    [enc setBuffer:s_grid          offset:0 atIndex:0];
    [enc setBuffer:s_father        offset:0 atIndex:1];
    [enc setBuffer:s_uold          offset:0 atIndex:2];
#ifdef MHD
    [enc setBuffer:s_bold          offset:0 atIndex:3];
    [enc setBytes:&head_idx        length:sizeof(int)   atIndex:4];
    [enc setBytes:&num_octs        length:sizeof(int)   atIndex:5];
    [enc setBytes:&internal_energy length:sizeof(int)   atIndex:6];
    [enc setBytes:&gamma           length:sizeof(float) atIndex:7];
    [enc setBytes:&smallr          length:sizeof(float) atIndex:8];
    [enc setBytes:&smallc2         length:sizeof(float) atIndex:9];
#else
    [enc setBytes:&head_idx        length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs        length:sizeof(int)   atIndex:4];
    [enc setBytes:&internal_energy length:sizeof(int)   atIndex:5];
    [enc setBytes:&gamma           length:sizeof(float) atIndex:6];
    [enc setBytes:&smallr          length:sizeof(float) atIndex:7];
    [enc setBytes:&smallc2         length:sizeof(float) atIndex:8];
#endif
    [enc dispatchThreadgroups:MTLSizeMake(num_tg, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_alloc_refine — allocate device buffers for AMR refinement, sorting,
 * ghost-zone cache, and per-level Hilbert parameters.
 * All buffers are MTLResourceStorageModeShared (unified memory) and zeroed.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_alloc_refine(int ngridmax, int ncachemax, int nlevelmax)
{
    int ntotal    = ngridmax + ncachemax;
    int npartial  = (ntotal    + 255) / 256;  /* ps0: ceil(n/256)   */
    int npartial2 = (npartial  + 255) / 256;  /* ps1: ceil(n/256^2) */
    int npartial3 = (npartial2 + 255) / 256;  /* ps2: ceil(n/256^3) */
    int nlevels   = nlevelmax + 1;

    NSUInteger flag_bytes     = (NSUInteger)8 * ntotal   * sizeof(int);  /* 8 = twotondim NDIM=3 */
    NSUInteger oct_bytes      = (NSUInteger)ntotal        * sizeof(int);  /* father/swap/prefix */
    NSUInteger partial_bytes  = (NSUInteger)npartial      * sizeof(int);
    NSUInteger partial2_bytes = (NSUInteger)npartial2     * sizeof(int);
    NSUInteger partial3_bytes = (NSUInteger)npartial3     * sizeof(int);
    NSUInteger partial4_bytes = sizeof(int);                              /* 1-element dummy sink */
    NSUInteger scalar_bytes  = sizeof(int);
    NSUInteger ckey_bytes    = (NSUInteger)nlevels      * sizeof(int);
    NSUInteger koff_bytes    = (NSUInteger)nlevels      * sizeof(long);
    NSUInteger box_bytes     = (NSUInteger)3 * nlevels  * sizeof(int);
    NSUInteger per_bytes     = 3 * sizeof(int);

    s_flag1        = [s_device newBufferWithLength:flag_bytes    options:MTLResourceStorageModeShared];
    s_flag2        = [s_device newBufferWithLength:flag_bytes    options:MTLResourceStorageModeShared];
    s_father       = [s_device newBufferWithLength:oct_bytes     options:MTLResourceStorageModeShared];
    s_swap_local   = [s_device newBufferWithLength:oct_bytes     options:MTLResourceStorageModeShared];
    s_swap_global  = [s_device newBufferWithLength:oct_bytes     options:MTLResourceStorageModeShared];
    s_prefix_sum   = [s_device newBufferWithLength:oct_bytes     options:MTLResourceStorageModeShared];
    s_partial_sums   = [s_device newBufferWithLength:partial_bytes  options:MTLResourceStorageModeShared];
    s_partial_sums_2 = [s_device newBufferWithLength:partial2_bytes options:MTLResourceStorageModeShared];
    s_partial_sums_3 = [s_device newBufferWithLength:partial3_bytes options:MTLResourceStorageModeShared];
    s_partial_sums_4 = [s_device newBufferWithLength:partial4_bytes options:MTLResourceStorageModeShared];
    s_ifree_dev        = [s_device newBufferWithLength:scalar_bytes options:MTLResourceStorageModeShared];
    s_ifree_cache_dev  = [s_device newBufferWithLength:scalar_bytes options:MTLResourceStorageModeShared];
    s_ckey_max_dev     = [s_device newBufferWithLength:ckey_bytes   options:MTLResourceStorageModeShared];
    s_key_off_dev      = [s_device newBufferWithLength:koff_bytes   options:MTLResourceStorageModeShared];
    s_box_ckey_min_dev = [s_device newBufferWithLength:box_bytes    options:MTLResourceStorageModeShared];
    s_box_ckey_max_dev = [s_device newBufferWithLength:box_bytes    options:MTLResourceStorageModeShared];
    s_periodic_dev     = [s_device newBufferWithLength:per_bytes    options:MTLResourceStorageModeShared];

    memset(s_flag1.contents,           0, flag_bytes);
    memset(s_flag2.contents,           0, flag_bytes);
    memset(s_father.contents,          0, oct_bytes);
    memset(s_swap_local.contents,      0, oct_bytes);
    memset(s_swap_global.contents,     0, oct_bytes);
    memset(s_prefix_sum.contents,      0, oct_bytes);
    memset(s_partial_sums.contents,    0, partial_bytes);
    memset(s_partial_sums_2.contents,  0, partial2_bytes);
    memset(s_partial_sums_3.contents,  0, partial3_bytes);
    memset(s_partial_sums_4.contents,  0, partial4_bytes);
    memset(s_ifree_dev.contents,       0, scalar_bytes);
    memset(s_ifree_cache_dev.contents, 0, scalar_bytes);
}

/* -----------------------------------------------------------------------
 * mtl_upload_level_params — copy per-level Hilbert parameters from host to
 * device.  Called once from metal_allocate_amr after init_amr populates them.
 *
 * box_ckey_min/max are Fortran (ndim, nlevelmax+1) column-major arrays;
 * the C pointer receives them in the same byte order, so kernels index as
 * box_ckey_min_dev[3*(lev-1) + d] (0-based) to get box_ckey_min(d+1, lev).
 * ----------------------------------------------------------------------- */
extern "C" void mtl_upload_level_params(void *ckey_max, void *key_off,
                                         void *box_ckey_min, void *box_ckey_max,
                                         int  *periodic,    int nlevelmax)
{
    int nlevels = nlevelmax + 1;
    memcpy(s_ckey_max_dev.contents,     ckey_max,     nlevels * sizeof(int));
    memcpy(s_key_off_dev.contents,      key_off,      nlevels * sizeof(long));
    memcpy(s_box_ckey_min_dev.contents, box_ckey_min, 3 * nlevels * sizeof(int));
    memcpy(s_box_ckey_max_dev.contents, box_ckey_max, 3 * nlevels * sizeof(int));
    memcpy(s_periodic_dev.contents,     periodic,     3 * sizeof(int));
}

/* -----------------------------------------------------------------------
 * scan_phase — helper that dispatches one pass of scan_block_kernel or
 * scan_fixup_kernel over the given buffer range.
 * offset and n are 0-based (converted from Fortran 1-based by the caller).
 * ----------------------------------------------------------------------- */
static void scan_phase(id<MTLComputePipelineState> pso,
                       id<MTLBuffer> data, id<MTLBuffer> psums,
                       int offset, int n,
                       id<MTLCommandBuffer> cmd)
{
    NSUInteger tg   = 256;
    NSUInteger nblk = ((NSUInteger)n + tg - 1) / tg;
    MTLSize tg_size   = {tg, 1, 1};
    MTLSize grid_size = {nblk, 1, 1};

    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:pso];
    [enc setBuffer:data  offset:0 atIndex:0];
    [enc setBuffer:psums offset:0 atIndex:1];
    [enc setBytes:&offset length:sizeof(int) atIndex:2];
    [enc setBytes:&n      length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:grid_size threadsPerThreadgroup:tg_size];
    [enc endEncoding];
}

/* -----------------------------------------------------------------------
 * mtl_prefix_scan — inclusive prefix scan of s_prefix_sum[offset..offset+n-1].
 * offset is 0-based (Fortran head_idx - 1).
 *
 * Three-phase approach (mirrors block_scan + scan(partial_sums) + uniform_add):
 *   Phase 1: scan_block_kernel  — within-block scans, block totals → s_partial_sums
 *   Phase 2: scan_block_kernel  — scan the partial_sums array (if n > 256)
 *   Phase 3: scan_fixup_kernel  — add partial_sums[bid-1] to block bid (if n > 256)
 *
 * All phases are submitted as separate command buffers and serialised via
 * waitUntilCompleted between phases (phase 2 needs phase 1 done; phase 3 needs
 * phase 2 done).
 * ----------------------------------------------------------------------- */
/* -----------------------------------------------------------------------
 * mtl_prefix_scan — inclusive prefix scan of s_prefix_sum[offset..offset+n-1].
 *
 * Three separate scratch buffers
 * (s_partial_sums / _2 / _3) hold block totals at successive levels so that
 * no scan pass ever aliases its data buffer with its partial-sums output.
 * s_partial_sums_4 is a 1-element dummy sink for the deepest single-block
 * pass (which must write its block total somewhere but the value is unused).
 *
 * Three cases, identical to the CUDA port:
 *   n ≤ 256^2 = 65,536          : 2-level (ps0 only)
 *   n ≤ 256^3 = 16,777,216      : 3-level (ps0, ps1)
 *   n ≤ INT_MAX                 : 4-level (ps0, ps1, ps2)
 * ----------------------------------------------------------------------- */
static void mtl_prefix_scan_buffer_impl(id<MTLBuffer> buf, id<MTLBuffer> ps1,
					id<MTLBuffer> ps2, id<MTLBuffer> ps3,
					id<MTLBuffer> ps4, int offset, int n, bool wait)
{
    if (n <= 0) return;

    const int BS  = 256;
    int nb0 = (n   + BS - 1) / BS;
    int nb1 = (nb0 + BS - 1) / BS;
    int nb2 = (nb1 + BS - 1) / BS;

#define CMD_WAIT(body) do { \
    id<MTLCommandBuffer> _c = [s_queue commandBuffer]; \
    body; \
    [_c commit]; if (wait) [_c waitUntilCompleted]; \
} while(0)

    if (n <= BS * BS) {
        /* ---- 2-level: n ≤ 65,536 ---------------------------------------- */
        CMD_WAIT(scan_phase(s_pso_scan_block, buf, ps1, offset, n, _c));
        if (nb0 == 1) goto done;
        /* single-block scan of ps0; total → ps1[0] (unused dummy) */
        CMD_WAIT(scan_phase(s_pso_scan_block, ps1, ps2, 0, nb0, _c));
        CMD_WAIT(scan_phase(s_pso_scan_fixup, buf, ps1, offset, n, _c));

    } else if (n <= BS * BS * BS) {
        /* ---- 3-level: n ≤ 16,777,216 ------------------------------------- */
        CMD_WAIT(scan_phase(s_pso_scan_block, buf, ps1, offset, n, _c));
        /* scan ps0 with nb1 blocks; block totals → ps1 */
        CMD_WAIT(scan_phase(s_pso_scan_block, ps1, ps2, 0, nb0, _c));
        if (nb1 > 1) {
            /* single-block scan of ps1; total → ps2[0] (unused dummy) */
            CMD_WAIT(scan_phase(s_pso_scan_block, ps2, ps3, 0, nb1, _c));
            /* fixup ps0 using ps1 */
            CMD_WAIT(scan_phase(s_pso_scan_fixup, ps1, ps2, 0, nb0, _c));
        }
        CMD_WAIT(scan_phase(s_pso_scan_fixup, buf, ps1, offset, n, _c));

    } else {
        /* ---- 4-level: n ≤ INT_MAX ---------------------------------------- */
        CMD_WAIT(scan_phase(s_pso_scan_block, buf, ps1, offset, n, _c));
        /* scan ps0 with nb1 blocks; block totals → ps1 */
        CMD_WAIT(scan_phase(s_pso_scan_block, ps1, ps2, 0, nb0, _c));
        /* scan ps1 with nb2 blocks; block totals → ps2 */
        CMD_WAIT(scan_phase(s_pso_scan_block, ps2, ps3, 0, nb1, _c));
        if (nb2 > 1) {
            /* single-block scan of ps2; total → ps3[0] (unused dummy) */
            CMD_WAIT(scan_phase(s_pso_scan_block, ps3, ps4, 0, nb2, _c));
            /* fixup ps1 using ps2 */
            CMD_WAIT(scan_phase(s_pso_scan_fixup, ps2, ps3, 0, nb1, _c));
        }
        /* fixup ps0 using ps1 */
        CMD_WAIT(scan_phase(s_pso_scan_fixup, ps1, ps2, 0, nb0, _c));
        CMD_WAIT(scan_phase(s_pso_scan_fixup, buf, ps1, offset, n, _c));
    }

done:;
#undef CMD_WAIT
}

static void mtl_prefix_scan_buffer(id<MTLBuffer> buf, id<MTLBuffer> ps1,
				   id<MTLBuffer> ps2, id<MTLBuffer> ps3,
				   id<MTLBuffer> ps4, int offset, int n)
{
    mtl_prefix_scan_buffer_impl(buf, ps1, ps2, ps3, ps4, offset, n, true);
}

static void mtl_prefix_scan_buffer_async(id<MTLBuffer> buf,
					 id<MTLBuffer> ps1,
					 id<MTLBuffer> ps2,
					 id<MTLBuffer> ps3,
					 id<MTLBuffer> ps4,
					 int offset, int n)
{
    mtl_prefix_scan_buffer_impl(buf, ps1, ps2, ps3, ps4, offset, n, false);
}

static void mtl_prefix_scan_buffer_cb(id<MTLBuffer> buf,
				      id<MTLBuffer> ps1,
				      id<MTLBuffer> ps2,
				      id<MTLBuffer> ps3,
				      id<MTLBuffer> ps4,
				      int offset, int n,
				      id<MTLCommandBuffer> cmd)
{
    if (n <= 0) return;

    const int BS  = 256;
    int nb0 = (n   + BS - 1) / BS;
    int nb1 = (nb0 + BS - 1) / BS;
    int nb2 = (nb1 + BS - 1) / BS;

    if (n <= BS * BS) {
        /* ---- 2-level: n ≤ 65,536 ---------------------------------------- */
        scan_phase(s_pso_scan_block, buf, ps1, offset, n, cmd);
        if (nb0 == 1) return;
        /* single-block scan of ps0; total → ps1[0] (unused dummy) */
        scan_phase(s_pso_scan_block, ps1, ps2, 0, nb0, cmd);
        scan_phase(s_pso_scan_fixup, buf, ps1, offset, n, cmd);

    } else if (n <= BS * BS * BS) {
        /* ---- 3-level: n ≤ 16,777,216 ------------------------------------- */
        scan_phase(s_pso_scan_block, buf, ps1, offset, n, cmd);
        /* scan ps0 with nb1 blocks; block totals → ps1 */
        scan_phase(s_pso_scan_block, ps1, ps2, 0, nb0, cmd);
        if (nb1 > 1) {
            /* single-block scan of ps1; total → ps2[0] (unused dummy) */
            scan_phase(s_pso_scan_block, ps2, ps3, 0, nb1, cmd);
            /* fixup ps0 using ps1 */
            scan_phase(s_pso_scan_fixup, ps1, ps2, 0, nb0, cmd);
        }
        scan_phase(s_pso_scan_fixup, buf, ps1, offset, n, cmd);

    } else {
        /* ---- 4-level: n ≤ INT_MAX ---------------------------------------- */
        scan_phase(s_pso_scan_block, buf, ps1, offset, n, cmd);
        /* scan ps0 with nb1 blocks; block totals → ps1 */
        scan_phase(s_pso_scan_block, ps1, ps2, 0, nb0, cmd);
        /* scan ps1 with nb2 blocks; block totals → ps2 */
        scan_phase(s_pso_scan_block, ps2, ps3, 0, nb1, cmd);
        if (nb2 > 1) {
            /* single-block scan of ps2; total → ps3[0] (unused dummy) */
            scan_phase(s_pso_scan_block, ps3, ps4, 0, nb2, cmd);
            /* fixup ps1 using ps2 */
            scan_phase(s_pso_scan_fixup, ps2, ps3, 0, nb1, cmd);
        }
        /* fixup ps0 using ps1 */
        scan_phase(s_pso_scan_fixup, ps1, ps2, 0, nb0, cmd);
        scan_phase(s_pso_scan_fixup, buf, ps1, offset, n, cmd);
    }
}

extern "C" void mtl_prefix_scan(int offset, int n)
{
    mtl_prefix_scan_buffer(s_prefix_sum, s_partial_sums, s_partial_sums_2,
			   s_partial_sums_3, s_partial_sums_4, offset, n);
}

/* -----------------------------------------------------------------------
 * mtl_get_prefix_total — read the inclusive sum of prefix_sum[offset..offset+n-1]
 * after mtl_prefix_scan.  Equivalent to get_total_sum in gpu_scan.cuf.
 * Returns prefix_sum[offset + n - 1] directly from shared memory.
 * ----------------------------------------------------------------------- */
extern "C" int mtl_get_prefix_total(int offset, int n)
{
    if (n <= 0) return 0;
    return ((int *)s_prefix_sum.contents)[offset + n - 1];
}

/* -----------------------------------------------------------------------
 * dispatch_2d_flag — helper for the 2D flag kernels (8 cells × 16 octs).
 * Mirrors CUDA <<<dim3(N,1,1), dim3(8,16,1)>>>.
 * ----------------------------------------------------------------------- */
static void dispatch_2d_flag(id<MTLComputePipelineState> pso,
                              id<MTLComputeCommandEncoder> enc,
                              int num_octs)
{
    NSUInteger tg128  = 128;   /* 8 × 16 = 128 */
    NSUInteger nblk   = ((NSUInteger)num_octs + 15) / 16;
    MTLSize tg_size   = {tg128, 1, 1};
    MTLSize grid_size = {nblk,  1, 1};
    [enc setComputePipelineState:pso];
    [enc dispatchThreadgroups:grid_size threadsPerThreadgroup:tg_size];
}

/* -----------------------------------------------------------------------
 * mtl_reset_flag1 — zero flag1 for [head_idx .. head_idx+num_octs-1].
 * ----------------------------------------------------------------------- */
extern "C" void mtl_reset_flag1(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_flag1];
    [enc setBuffer:s_flag1    offset:0 atIndex:0];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:1];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:2];
    NSUInteger nblk = ((NSUInteger)num_octs + 15) / 16;
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_reset_flag2 — zero flag2 for [head_idx .. head_idx+num_octs-1].
 * ----------------------------------------------------------------------- */
extern "C" void mtl_reset_flag2(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_flag2];
    [enc setBuffer:s_flag2    offset:0 atIndex:0];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:1];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:2];
    NSUInteger nblk = ((NSUInteger)num_octs + 15) / 16;
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_init_flag — flag parent cell for each fine oct at ilevel+1.
 * head_idx / num_octs refer to the FINE level (ilevel+1).
 * ----------------------------------------------------------------------- */
extern "C" void mtl_init_flag(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg128  = 128;
    NSUInteger nblk   = ((NSUInteger)num_octs + tg128 - 1) / tg128;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_init_flag];
    [enc setBuffer:s_flag1    offset:0 atIndex:0];
    [enc setBuffer:s_grid     offset:0 atIndex:1];
    [enc setBuffer:s_father   offset:0 atIndex:2];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg128,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_count_flag1 — reduce sum of flag1 cells; returns the count.
 * head_idx / num_octs refer to the coarse level (ilevel).
 * ----------------------------------------------------------------------- */
extern "C" int mtl_count_flag1(int head_idx, int num_octs)
{
    if (num_octs <= 0) return 0;
    /* Allocate a zeroed one-int buffer for the atomic result */
    int zero = 0;
    id<MTLBuffer> result_buf =
        [s_device newBufferWithBytes:&zero
                              length:sizeof(int)
                             options:MTLResourceStorageModeShared];

    NSUInteger tg1024 = 1024;
    NSUInteger total_cells = (NSUInteger)num_octs * 8;
    NSUInteger nblk  = (total_cells + tg1024 - 1) / tg1024;

    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_count_flag1];
    [enc setBuffer:s_flag1    offset:0 atIndex:0];
    [enc setBuffer:result_buf offset:0 atIndex:1];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg1024,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    int count = *(int *)result_buf.contents;
    return count;
}

/* -----------------------------------------------------------------------
 * mtl_hydro_flag — gradient-based hydro and MHD refinement criteria.
 * head_idx / num_octs: octs at ilevel.  No GRAV.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_hydro_flag(int head_idx, int num_octs,
                                float gamma, float smallr, float smallc2,
                                float err_grad_d, float err_grad_p,
                                float floor_d,   float floor_p
#ifdef MHD
                                , float err_grad_b2, float floor_b2,
                                float err_grad_A, float floor_A,
                                float err_grad_B, float floor_B,
                                float err_grad_C, float floor_C
#endif
                                )
{
    if (num_octs <= 0) return;
    NSUInteger nblk = ((NSUInteger)num_octs + 15) / 16;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_hydro_flag];
    [enc setBuffer:s_flag1      offset:0 atIndex:0];
    [enc setBuffer:s_nbor       offset:0 atIndex:1];
    [enc setBuffer:s_uold       offset:0 atIndex:2];
    [enc setBytes:&head_idx     length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs     length:sizeof(int)   atIndex:4];
    [enc setBytes:&gamma        length:sizeof(float) atIndex:5];
    [enc setBytes:&smallr       length:sizeof(float) atIndex:6];
    [enc setBytes:&smallc2      length:sizeof(float) atIndex:7];
    [enc setBytes:&err_grad_d   length:sizeof(float) atIndex:8];
    [enc setBytes:&err_grad_p   length:sizeof(float) atIndex:9];
    [enc setBytes:&floor_d      length:sizeof(float) atIndex:10];
    [enc setBytes:&floor_p      length:sizeof(float) atIndex:11];
#ifdef MHD
    [enc setBuffer:s_bold       offset:0 atIndex:12];
    [enc setBytes:&err_grad_b2  length:sizeof(float) atIndex:13];
    [enc setBytes:&floor_b2     length:sizeof(float) atIndex:14];
    [enc setBytes:&err_grad_A   length:sizeof(float) atIndex:15];
    [enc setBytes:&floor_A      length:sizeof(float) atIndex:16];
    [enc setBytes:&err_grad_B   length:sizeof(float) atIndex:17];
    [enc setBytes:&floor_B      length:sizeof(float) atIndex:18];
    [enc setBytes:&err_grad_C   length:sizeof(float) atIndex:19];
    [enc setBytes:&floor_C      length:sizeof(float) atIndex:20];
#endif
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_count_neighbors — write flag2[cell,oct] = # flagged face neighbours.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_count_neighbors(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger nblk = ((NSUInteger)num_octs + 15) / 16;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_count_neighbors];
    [enc setBuffer:s_flag2    offset:0 atIndex:0];
    [enc setBuffer:s_flag1    offset:0 atIndex:1];
    [enc setBuffer:s_nbor     offset:0 atIndex:2];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_flag_count — promote flag1 if flag2 >= num_nbors; clear flag2 if
 * flag1 is already set.  num_nbors = n_nbor(idim) = {1,2,2} for idim=1..3.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_flag_count(int head_idx, int num_octs, int num_nbors)
{
    if (num_octs <= 0) return;
    NSUInteger nblk = ((NSUInteger)num_octs + 15) / 16;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_flag_count];
    [enc setBuffer:s_flag1      offset:0 atIndex:0];
    [enc setBuffer:s_flag2      offset:0 atIndex:1];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:3];
    [enc setBytes:&num_nbors    length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_enforce_rules — clear flag1 if any nbor slot is 0 or > ngridmax.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_enforce_rules(int head_idx, int num_octs, int ngridmax)
{
    if (num_octs <= 0) return;
    NSUInteger tg128 = 128;
    NSUInteger nblk  = ((NSUInteger)num_octs + tg128 - 1) / tg128;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_enforce_rules];
    [enc setBuffer:s_flag1    offset:0 atIndex:0];
    [enc setBuffer:s_nbor     offset:0 atIndex:1];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:3];
    [enc setBytes:&ngridmax   length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg128,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

extern "C" void mtl_enforce_subgrid(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_enforce_subgrid];
    [enc setBuffer:s_flag1 offset:0 atIndex:0];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:1];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:2];
    NSUInteger nblk = ((NSUInteger)num_octs + 15) / 16;
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * Batched flag helpers — one command buffer per operation, using
 * s_sort_fence for inter-encoder ordering and s_count_buf for the
 * atomic reduction result (zeroed via blit at the start of each cmd buf).
 *
 * mtl_init_flag_batch  : reset_flag1 [→ init_flag] → count_flag1  (1 sync)
 * mtl_user_flag_batch  : hydro_flag  → count_flag1                (1 sync)
 * mtl_smooth_flag_batch: 3×(count_neighbors → flag_count) → count (1 sync)
 * ----------------------------------------------------------------------- */

static int run_count_enc(id<MTLCommandBuffer> cmd,
                         int head_idx, int num_octs)
{
    NSUInteger tg1024 = 1024;
    NSUInteger total_cells = (NSUInteger)num_octs * 8;
    NSUInteger nblk  = (total_cells + tg1024 - 1) / tg1024;
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc waitForFence:s_sort_fence];
    [enc setComputePipelineState:s_pso_count_flag1];
    [enc setBuffer:s_flag1    offset:0 atIndex:0];
    [enc setBuffer:s_count_buf offset:0 atIndex:1];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg1024,1,1}];
    [enc endEncoding];   /* last encoder — no outgoing fence */
    return 0;            /* caller reads s_count_buf after waitUntilCompleted */
}

extern "C" int mtl_init_flag_batch(int head_coarse, int noct_coarse,
                                    int head_fine,   int noct_fine)
{
    if (noct_coarse <= 0) return 0;

    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];

    /* Zero the persistent count buffer */
    { id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
      [blit fillBuffer:s_count_buf range:NSMakeRange(0, sizeof(int)) value:0];
      [blit updateFence:s_sort_fence];
      [blit endEncoding]; }

    /* reset_flag1 for coarse level */
    { NSUInteger nblk = ((NSUInteger)noct_coarse + 15) / 16;
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc waitForFence:s_sort_fence];
      [enc setComputePipelineState:s_pso_reset_flag1];
      [enc setBuffer:s_flag1    offset:0 atIndex:0];
      [enc setBytes:&head_coarse length:sizeof(int) atIndex:1];
      [enc setBytes:&noct_coarse length:sizeof(int) atIndex:2];
      [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
      [enc updateFence:s_sort_fence];
      [enc endEncoding]; }

    /* init_flag: propagate fine-level flags to coarse parent cells */
    if (noct_fine > 0) {
        NSUInteger nblk = ((NSUInteger)noct_fine + 127) / 128;
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc waitForFence:s_sort_fence];
        [enc setComputePipelineState:s_pso_init_flag];
        [enc setBuffer:s_flag1    offset:0 atIndex:0];
        [enc setBuffer:s_grid     offset:0 atIndex:1];
        [enc setBuffer:s_father   offset:0 atIndex:2];
        [enc setBytes:&head_fine  length:sizeof(int) atIndex:3];
        [enc setBytes:&noct_fine  length:sizeof(int) atIndex:4];
        [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
        [enc updateFence:s_sort_fence];
        [enc endEncoding]; }

    run_count_enc(cmd, head_coarse, noct_coarse);
    [cmd commit]; [cmd waitUntilCompleted];
    return *(int *)s_count_buf.contents;
}

extern "C" int mtl_user_flag_batch(int head_idx, int num_octs,
                                    float gamma,      float smallr,  float smallc2,
                                    float err_grad_d, float err_grad_p,
                                    float floor_d,    float floor_p,
                                    float mass_sph,   float m_refine, float jeans_refine,
                                    float factG,      float dx_loc
#ifdef MHD
                                    , float err_grad_b2, float floor_b2,
                                    float err_grad_A, float floor_A,
                                    float err_grad_B, float floor_B,
                                    float err_grad_C, float floor_C
#endif
                                    )
{
    if (num_octs <= 0) return 0;

    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];

    /* Zero count buffer */
    { id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
      [blit fillBuffer:s_count_buf range:NSMakeRange(0, sizeof(int)) value:0];
      [blit updateFence:s_sort_fence];
      [blit endEncoding]; }

    /* poisson_flag and hydro_flag */
    { NSUInteger nblk = ((NSUInteger)num_octs + 15) / 16;
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc waitForFence:s_sort_fence];

      /* 1. poisson_flag if gravity is active (s_pso_poisson_flag) */
      if (s_pso_poisson_flag) {
          [enc setComputePipelineState:s_pso_poisson_flag];
          [enc setBuffer:s_flag1      offset:0 atIndex:0];
          [enc setBuffer:s_nref       offset:0 atIndex:1];
          [enc setBuffer:s_uold       offset:0 atIndex:2];
#ifdef MHD
          [enc setBuffer:s_bold       offset:0 atIndex:3];
#else
          [enc setBuffer:nil          offset:0 atIndex:3]; /* bold: nil for MHD=0 */
#endif
          [enc setBytes:&head_idx     length:sizeof(int)   atIndex:4];
          [enc setBytes:&num_octs     length:sizeof(int)   atIndex:5];
          [enc setBytes:&gamma        length:sizeof(float) atIndex:6];
          [enc setBytes:&smallr       length:sizeof(float) atIndex:7];
          [enc setBytes:&smallc2      length:sizeof(float) atIndex:8];
          [enc setBytes:&mass_sph     length:sizeof(float) atIndex:9];
          [enc setBytes:&m_refine     length:sizeof(float) atIndex:10];
          [enc setBytes:&jeans_refine length:sizeof(float) atIndex:11];
          [enc setBytes:&factG        length:sizeof(float) atIndex:12];
          [enc setBytes:&dx_loc       length:sizeof(float) atIndex:13];
          [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
      }

      /* 2. hydro_flag */
      [enc setComputePipelineState:s_pso_hydro_flag];
      [enc setBuffer:s_flag1      offset:0 atIndex:0];
      [enc setBuffer:s_nbor       offset:0 atIndex:1];
      [enc setBuffer:s_uold       offset:0 atIndex:2];
      [enc setBytes:&head_idx     length:sizeof(int)   atIndex:3];
      [enc setBytes:&num_octs     length:sizeof(int)   atIndex:4];
      [enc setBytes:&gamma        length:sizeof(float) atIndex:5];
      [enc setBytes:&smallr       length:sizeof(float) atIndex:6];
      [enc setBytes:&smallc2      length:sizeof(float) atIndex:7];
      [enc setBytes:&err_grad_d   length:sizeof(float) atIndex:8];
      [enc setBytes:&err_grad_p   length:sizeof(float) atIndex:9];
      [enc setBytes:&floor_d      length:sizeof(float) atIndex:10];
      [enc setBytes:&floor_p      length:sizeof(float) atIndex:11];
#ifdef MHD
      [enc setBuffer:s_bold       offset:0 atIndex:12];
      [enc setBytes:&err_grad_b2  length:sizeof(float) atIndex:13];
      [enc setBytes:&floor_b2     length:sizeof(float) atIndex:14];
      [enc setBytes:&err_grad_A   length:sizeof(float) atIndex:15];
      [enc setBytes:&floor_A      length:sizeof(float) atIndex:16];
      [enc setBytes:&err_grad_B   length:sizeof(float) atIndex:17];
      [enc setBytes:&floor_B      length:sizeof(float) atIndex:18];
      [enc setBytes:&err_grad_C   length:sizeof(float) atIndex:19];
      [enc setBytes:&floor_C      length:sizeof(float) atIndex:20];
#endif
      [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];

      [enc updateFence:s_sort_fence];
      [enc endEncoding]; }

    run_count_enc(cmd, head_idx, num_octs);
    [cmd commit]; [cmd waitUntilCompleted];
    return *(int *)s_count_buf.contents;
}

extern "C" int mtl_smooth_flag_batch(int head_idx, int num_octs)
{
    if (num_octs <= 0) return 0;

    /* n_nbor = {1, 2, 2} for NDIM=3 — mirrors metal_runner.f90 metal_smooth_flag */
    static const int n_nbor[3] = {1, 2, 2};
    NSUInteger nblk = ((NSUInteger)num_octs + 15) / 16;

    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];

    /* Zero count buffer */
    { id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
      [blit fillBuffer:s_count_buf range:NSMakeRange(0, sizeof(int)) value:0];
      [blit updateFence:s_sort_fence];
      [blit endEncoding]; }

    /* 3 dilatation passes: count_neighbors → flag_count */
    for (int idim = 0; idim < 3; idim++) {
        int nn = n_nbor[idim];

        /* count_neighbors: flag2[cell,oct] = # flagged face-adjacent nbors */
        { id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
          [enc waitForFence:s_sort_fence];
          [enc setComputePipelineState:s_pso_count_neighbors];
          [enc setBuffer:s_flag2    offset:0 atIndex:0];
          [enc setBuffer:s_flag1    offset:0 atIndex:1];
          [enc setBuffer:s_nbor     offset:0 atIndex:2];
          [enc setBytes:&head_idx   length:sizeof(int) atIndex:3];
          [enc setBytes:&num_octs   length:sizeof(int) atIndex:4];
          [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
          [enc updateFence:s_sort_fence];
          [enc endEncoding]; }

        /* flag_count: promote flag1 if flag2 >= nn; clear flag2 if flag1 set */
        { id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
          [enc waitForFence:s_sort_fence];
          [enc setComputePipelineState:s_pso_flag_count];
          [enc setBuffer:s_flag1    offset:0 atIndex:0];
          [enc setBuffer:s_flag2    offset:0 atIndex:1];
          [enc setBytes:&head_idx   length:sizeof(int) atIndex:2];
          [enc setBytes:&num_octs   length:sizeof(int) atIndex:3];
          [enc setBytes:&nn         length:sizeof(int) atIndex:4];
          [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{128,1,1}];
          [enc updateFence:s_sort_fence];
          [enc endEncoding]; }
    }

    run_count_enc(cmd, head_idx, num_octs);
    [cmd commit]; [cmd waitUntilCompleted];
    return *(int *)s_count_buf.contents;
}

/* -----------------------------------------------------------------------
 * mtl_build_father — populate father[oct] = 1-based parent oct for each
 * oct at the given level.  Uses the hash table already populated by
 * mtl_insert_hash.  Needed before init_flag_kernel can run.
 * ckey_max_l and key_off_l are the per-level Hilbert parameters for
 * the PARENT level (ilevel - 1).
 * ----------------------------------------------------------------------- */
extern "C" void mtl_build_father(int head_idx, int num_octs,
                                  int hash_size,
                                  int ckey_max_l, long key_off_l)
{
    if (num_octs <= 0) return;
    NSUInteger tg128 = 128;
    NSUInteger nblk  = ((NSUInteger)num_octs + tg128 - 1) / tg128;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_update_father];
    [enc setBuffer:s_grid       offset:0 atIndex:0];
    [enc setBuffer:s_father     offset:0 atIndex:1];
    [enc setBuffer:s_hash_key   offset:0 atIndex:2];
    [enc setBuffer:s_hash_val   offset:0 atIndex:3];
    [enc setBytes:&hash_size    length:sizeof(int)  atIndex:4];
    [enc setBytes:&ckey_max_l   length:sizeof(int)  atIndex:5];
    [enc setBytes:&key_off_l    length:sizeof(long) atIndex:6];
    [enc setBytes:&head_idx     length:sizeof(int)  atIndex:7];
    [enc setBytes:&num_octs     length:sizeof(int)  atIndex:8];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg128,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* =======================================================================
 * AMR refine / sort / cache bridge functions
 * All mirrors of gpu_refine subroutines from gpu_runner.cuf.
 * ======================================================================= */

/* --- ifree_dev / ifree_cache_dev direct access (unified memory) -------- */

extern "C" void mtl_set_ifree(int val)
{
    ((int *)s_ifree_dev.contents)[0] = val;
}

extern "C" int mtl_get_ifree(void)
{
    return ((int *)s_ifree_dev.contents)[0];
}

extern "C" void mtl_set_ifree_cache(int val)
{
    ((int *)s_ifree_cache_dev.contents)[0] = val;
}

extern "C" int mtl_get_ifree_cache(void)
{
    return ((int *)s_ifree_cache_dev.contents)[0];
}

/* CPU-side increment — no kernel dispatch needed for unified memory. */
extern "C" void mtl_advance_ifree_cache(int new_noct)
{
    ((int *)s_ifree_cache_dev.contents)[0] += new_noct;
}

/* --- refine_kernel: create child octs for flagged cells ---------------- */

extern "C" void mtl_refine_cells(int head_idx, int num_octs, int hash_size,
                                  int interpol_var, int interpol_type, float smallr)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, total = (NSUInteger)num_octs * 8;
    NSUInteger nblk = (total + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_refine];
    [enc setBuffer:s_grid      offset:0 atIndex:0];
    [enc setBuffer:s_flag1     offset:0 atIndex:1];
    [enc setBuffer:s_uold      offset:0 atIndex:2];
#ifdef MHD
    [enc setBuffer:s_bold      offset:0 atIndex:3];
    [enc setBuffer:s_ifree_dev offset:0 atIndex:4];
    [enc setBuffer:s_hash_key  offset:0 atIndex:5];
    [enc setBuffer:s_hash_val  offset:0 atIndex:6];
    [enc setBuffer:s_ckey_max_dev offset:0 atIndex:7];
    [enc setBuffer:s_key_off_dev offset:0 atIndex:8];
    [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:9];
    [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:10];
    [enc setBuffer:s_periodic_dev offset:0 atIndex:11];
    [enc setBytes:&hash_size length:sizeof(int) atIndex:12];
    [enc setBytes:&interpol_var length:sizeof(int) atIndex:13];
    [enc setBytes:&interpol_type length:sizeof(int) atIndex:14];
    [enc setBytes:&smallr length:sizeof(float) atIndex:15];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:16];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:17];
#else
    [enc setBuffer:s_ifree_dev offset:0 atIndex:3];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:4];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:5];
#ifdef GRAV
    [enc setBuffer:s_f_grav    offset:0 atIndex:6];
    [enc setBuffer:s_phi       offset:0 atIndex:7];
    [enc setBuffer:s_phi_old   offset:0 atIndex:8];
#endif
#endif
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- derefine_kernel: free child octs whose parent is no longer flagged */

extern "C" void mtl_derefine_cells(int head_idx, int num_octs, int hash_size)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_derefine];
    [enc setBuffer:s_grid         offset:0 atIndex:0];
    [enc setBuffer:s_flag1        offset:0 atIndex:1];
    [enc setBuffer:s_hash_key     offset:0 atIndex:2];
    [enc setBuffer:s_hash_val     offset:0 atIndex:3];
    [enc setBuffer:s_ckey_max_dev offset:0 atIndex:4];
    [enc setBuffer:s_key_off_dev  offset:0 atIndex:5];
    [enc setBytes:&hash_size length:sizeof(int) atIndex:6];
    [enc setBytes:&head_idx  length:sizeof(int) atIndex:7];
    [enc setBytes:&num_octs  length:sizeof(int) atIndex:8];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- free_hash_kernel: wipe hash entries for a range of octs ----------- */

extern "C" void mtl_free_hash_range(int head_idx, int num_octs, int hash_size)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_free_hash];
    [enc setBuffer:s_grid         offset:0 atIndex:0];
    [enc setBuffer:s_hash_key     offset:0 atIndex:1];
    [enc setBuffer:s_hash_val     offset:0 atIndex:2];
    [enc setBuffer:s_ckey_max_dev offset:0 atIndex:3];
    [enc setBuffer:s_key_off_dev  offset:0 atIndex:4];
    [enc setBytes:&hash_size length:sizeof(int) atIndex:5];
    [enc setBytes:&head_idx  length:sizeof(int) atIndex:6];
    [enc setBytes:&num_octs  length:sizeof(int) atIndex:7];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- update_hash_kernel: update hash entries after sort rearrangement --- */

extern "C" void mtl_update_hash_range(int head_idx, int num_octs, int hash_size)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_update_hash];
    [enc setBuffer:s_grid         offset:0 atIndex:0];
    [enc setBuffer:s_hash_key     offset:0 atIndex:1];
    [enc setBuffer:s_hash_val     offset:0 atIndex:2];
    [enc setBuffer:s_ckey_max_dev offset:0 atIndex:3];
    [enc setBuffer:s_key_off_dev  offset:0 atIndex:4];
    [enc setBytes:&hash_size length:sizeof(int) atIndex:5];
    [enc setBytes:&head_idx  length:sizeof(int) atIndex:6];
    [enc setBytes:&num_octs  length:sizeof(int) atIndex:7];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- insert_hash_all_kernel: insert newly created octs into hash ------- */

extern "C" void mtl_insert_hash_all(int head_idx, int num_octs, int hash_size)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_insert_hash_all];
    [enc setBuffer:s_grid         offset:0 atIndex:0];
    [enc setBuffer:s_hash_key     offset:0 atIndex:1];
    [enc setBuffer:s_hash_val     offset:0 atIndex:2];
    [enc setBuffer:s_ckey_max_dev offset:0 atIndex:3];
    [enc setBuffer:s_key_off_dev  offset:0 atIndex:4];
    [enc setBytes:&hash_size length:sizeof(int) atIndex:5];
    [enc setBytes:&head_idx  length:sizeof(int) atIndex:6];
    [enc setBytes:&num_octs  length:sizeof(int) atIndex:7];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- mtl_reset_hash: mirror gpu_reset_hash in gpu_runner.cuf ----------- *
 * Called once per coarse step (ilevel==levelmin).                         *
 * 1. Blit-zero s_hash_key / s_hash_val.                                   *
 * 2. Re-insert all real octs  [1 .. ifree-1].                             *
 * 3. Re-insert all cache octs [ngridmax+1 .. ngridmax+ifree_cache-1].    */
extern "C" void mtl_reset_hash(int ifree, int ngridmax, int ifree_cache,
                                int hash_size)
{
    /* Step 1: zero the hash table (GPU-side blit, ordered on the queue). */
    {
        id<MTLCommandBuffer>      cmd  = [s_queue commandBuffer];
        id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
        [blit fillBuffer:s_hash_key range:NSMakeRange(0, (NSUInteger)hash_size * sizeof(long)) value:0];
        [blit fillBuffer:s_hash_val range:NSMakeRange(0, (NSUInteger)hash_size * sizeof(int))  value:0];
        [blit endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];
    }
    /* Step 2: insert real octs [1 .. ifree-1]. */
    {
        int num = ifree - 1;
        if (num > 0) {
            NSUInteger tg = 128, nblk = ((NSUInteger)num + tg - 1) / tg;
            id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:s_pso_insert_hash_all];
            [enc setBuffer:s_grid         offset:0 atIndex:0];
            [enc setBuffer:s_hash_key     offset:0 atIndex:1];
            [enc setBuffer:s_hash_val     offset:0 atIndex:2];
            [enc setBuffer:s_ckey_max_dev offset:0 atIndex:3];
            [enc setBuffer:s_key_off_dev  offset:0 atIndex:4];
            [enc setBytes:&hash_size length:sizeof(int) atIndex:5];
            int head = 1;
            [enc setBytes:&head      length:sizeof(int) atIndex:6];
            [enc setBytes:&num       length:sizeof(int) atIndex:7];
            [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
    }
    /* Step 3: insert cache octs [ngridmax+1 .. ngridmax+ifree_cache-1]. */
    {
        int cache_head = ngridmax + 1;
        int cache_num  = ifree_cache - 1;
        if (cache_num > 0) {
            NSUInteger tg = 128, nblk = ((NSUInteger)cache_num + tg - 1) / tg;
            id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:s_pso_insert_hash_all];
            [enc setBuffer:s_grid         offset:0 atIndex:0];
            [enc setBuffer:s_hash_key     offset:0 atIndex:1];
            [enc setBuffer:s_hash_val     offset:0 atIndex:2];
            [enc setBuffer:s_ckey_max_dev offset:0 atIndex:3];
            [enc setBuffer:s_key_off_dev  offset:0 atIndex:4];
            [enc setBytes:&hash_size   length:sizeof(int) atIndex:5];
            [enc setBytes:&cache_head  length:sizeof(int) atIndex:6];
            [enc setBytes:&cache_num   length:sizeof(int) atIndex:7];
            [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
        }
    }
}

/* --- init_global_swap_table_kernel: identity permutation --------------- */

extern "C" void mtl_init_swap_table(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_init_swap_table];
    [enc setBuffer:s_swap_global offset:0 atIndex:0];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:1];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:2];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- init_prefix_sum_level_kernel: bit = (lev != ilevel) -------------- */

extern "C" void mtl_init_prefix_level(int head_idx, int num_octs, int ilevel)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_prefix_level];
    [enc setBuffer:s_grid        offset:0 atIndex:0];
    [enc setBuffer:s_swap_global offset:0 atIndex:1];
    [enc setBuffer:s_prefix_sum  offset:0 atIndex:2];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
    [enc setBytes:&ilevel   length:sizeof(int) atIndex:5];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- init_prefix_sum_bit_kernel: bit ibit of Hilbert key -------------- */

extern "C" void mtl_init_prefix_bit(int head_idx, int num_octs, int ibit)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_prefix_bit];
    [enc setBuffer:s_grid        offset:0 atIndex:0];
    [enc setBuffer:s_swap_global offset:0 atIndex:1];
    [enc setBuffer:s_prefix_sum  offset:0 atIndex:2];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
    [enc setBytes:&ibit     length:sizeof(int) atIndex:5];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- compute_local_swap_table_kernel: LSD scatter ---------------------- */

extern "C" void mtl_compute_local_swap(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_local_swap];
    [enc setBuffer:s_swap_local  offset:0 atIndex:0];
    [enc setBuffer:s_swap_global offset:0 atIndex:1];
    [enc setBuffer:s_prefix_sum  offset:0 atIndex:2];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- update_global_swap_table_kernel: apply local swap to global ------- */

extern "C" void mtl_update_global_swap(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_global_swap];
    [enc setBuffer:s_swap_global offset:0 atIndex:0];
    [enc setBuffer:s_swap_local  offset:0 atIndex:1];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- scan_phase_fenced: one scan pass in a new encoder with fence I/O ---
 * Waits on fence_in before dispatching, updates fence_out after.
 * fence_in and fence_out may be the same object (sequential chaining).   */
static void scan_phase_fenced(id<MTLComputePipelineState> pso,
                              id<MTLBuffer> data, id<MTLBuffer> psums,
                              int offset, int n,
                              id<MTLCommandBuffer> cmd,
                              id<MTLFence> fence_in, id<MTLFence> fence_out)
{
    NSUInteger tg = 256, nblk = ((NSUInteger)n + tg - 1) / tg;
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc waitForFence:fence_in];
    [enc setComputePipelineState:pso];
    [enc setBuffer:data  offset:0 atIndex:0];
    [enc setBuffer:psums offset:0 atIndex:1];
    [enc setBytes:&offset length:sizeof(int) atIndex:2];
    [enc setBytes:&n      length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc updateFence:fence_out];
    [enc endEncoding];
}

/* --- mtl_hilbert_sort_level: one command buffer per bit pass -----------
 * Each pass contains 5–7 compute encoders connected by s_sort_fence.
 * Reduces commit/waitUntilCompleted from 6*num_bits to num_bits — a 6×
 * reduction — while guaranteeing memory visibility via Metal fences.     */
extern "C" void mtl_hilbert_sort_level(int head_idx, int num_octs, int num_bits)
{
    if (num_octs <= 0 || num_bits <= 0) return;

    int offset    = head_idx - 1;
    int n         = num_octs;
    int nblk_sort = (n + 127) / 128;
    int nblk_scan = (n + 255) / 256;
    bool need_fixup = (nblk_scan > 1);

    for (int ibit = 0; ibit < num_bits; ibit++) {
        id<MTLCommandBuffer> cmd = [s_queue commandBuffer];

        /* 1. init_prefix_bit — no fence_in (first encoder in this buffer) */
        {
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:s_pso_prefix_bit];
            [enc setBuffer:s_grid        offset:0 atIndex:0];
            [enc setBuffer:s_swap_global offset:0 atIndex:1];
            [enc setBuffer:s_prefix_sum  offset:0 atIndex:2];
            [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
            [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
            [enc setBytes:&ibit     length:sizeof(int) atIndex:5];
            [enc dispatchThreadgroups:{(NSUInteger)nblk_sort,1,1}
                 threadsPerThreadgroup:{128,1,1}];
            [enc updateFence:s_sort_fence];
            [enc endEncoding];
        }
        /* 2. prefix scan: phase 1 */
        scan_phase_fenced(s_pso_scan_block, s_prefix_sum, s_partial_sums,
                          offset, n, cmd, s_sort_fence, s_sort_fence);
        if (need_fixup) {
            int nblk_scan2 = (nblk_scan + 255) / 256;
            /* phase 2: scan partial_sums; block totals → partial_sums_2 (separate buffer) */
            scan_phase_fenced(s_pso_scan_block, s_partial_sums, s_partial_sums_2,
                              0, nblk_scan, cmd, s_sort_fence, s_sort_fence);
            if (nblk_scan2 > 1) {
                /* nblk_scan > 256: need a third level.
                 * single-block scan of partial_sums_2; block totals → partial_sums_3. */
                scan_phase_fenced(s_pso_scan_block, s_partial_sums_2, s_partial_sums_3,
                                  0, nblk_scan2, cmd, s_sort_fence, s_sort_fence);
                /* fixup partial_sums using partial_sums_2 */
                scan_phase_fenced(s_pso_scan_fixup, s_partial_sums, s_partial_sums_2,
                                  0, nblk_scan, cmd, s_sort_fence, s_sort_fence);
            }
            /* phase 3: fixup prefix_sum using partial_sums */
            scan_phase_fenced(s_pso_scan_fixup, s_prefix_sum, s_partial_sums,
                              offset, n, cmd, s_sort_fence, s_sort_fence);
        }
        /* 3. compute_local_swap */
        {
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc waitForFence:s_sort_fence];
            [enc setComputePipelineState:s_pso_local_swap];
            [enc setBuffer:s_swap_local  offset:0 atIndex:0];
            [enc setBuffer:s_swap_global offset:0 atIndex:1];
            [enc setBuffer:s_prefix_sum  offset:0 atIndex:2];
            [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
            [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
            [enc dispatchThreadgroups:{(NSUInteger)nblk_sort,1,1}
                 threadsPerThreadgroup:{128,1,1}];
            [enc updateFence:s_sort_fence];
            [enc endEncoding];
        }
        /* 4. update_global_swap — last encoder, no outgoing fence needed */
        {
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc waitForFence:s_sort_fence];
            [enc setComputePipelineState:s_pso_global_swap];
            [enc setBuffer:s_swap_global offset:0 atIndex:0];
            [enc setBuffer:s_swap_local  offset:0 atIndex:1];
            [enc setBytes:&head_idx length:sizeof(int) atIndex:2];
            [enc setBytes:&num_octs length:sizeof(int) atIndex:3];
            [enc dispatchThreadgroups:{(NSUInteger)nblk_sort,1,1}
                 threadsPerThreadgroup:{128,1,1}];
            [enc endEncoding];
        }

        [cmd commit];
        [cmd waitUntilCompleted];
    }
}

/* --- sort_gather_grid_kernel: pack grid metadata into flag2 scratch ---- */

extern "C" void mtl_sort_gather_grid(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gather_grid];
    [enc setBuffer:s_flag2       offset:0 atIndex:0];
    [enc setBuffer:s_grid        offset:0 atIndex:1];
    [enc setBuffer:s_swap_global offset:0 atIndex:2];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- sort_scatter_grid_kernel: unpack flag2 → grid, recompute hkey ---- */

extern "C" void mtl_sort_scatter_grid(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_scatter_grid];
    [enc setBuffer:s_grid  offset:0 atIndex:0];
    [enc setBuffer:s_flag2 offset:0 atIndex:1];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- sort_gather_flag_kernel: flag2[:,oct] = flag1[:,swap_global[oct]] - */

extern "C" void mtl_sort_gather_flag(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gather_flag];
    [enc setBuffer:s_flag2       offset:0 atIndex:0];
    [enc setBuffer:s_flag1       offset:0 atIndex:1];
    [enc setBuffer:s_swap_global offset:0 atIndex:2];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- sort_scatter_flag_kernel: flag1[:,oct] = flag2[:,oct] ------------- */

extern "C" void mtl_sort_scatter_flag(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_scatter_flag];
    [enc setBuffer:s_flag1 offset:0 atIndex:0];
    [enc setBuffer:s_flag2 offset:0 atIndex:1];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- sort_gather_hydro_kernel: unew[:,oct] = uold[:,swap_global[oct]] -- */

extern "C" void mtl_sort_gather_hydro(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, total = (NSUInteger)num_octs * 8;
    NSUInteger nblk = (total + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gather_hydro];
    [enc setBuffer:s_unew        offset:0 atIndex:0];
    [enc setBuffer:s_uold        offset:0 atIndex:1];
#ifdef MHD
    [enc setBuffer:s_bnew        offset:0 atIndex:2];
    [enc setBuffer:s_bold        offset:0 atIndex:3];
    [enc setBuffer:s_swap_global offset:0 atIndex:4];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:5];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:6];
#else
    [enc setBuffer:s_swap_global offset:0 atIndex:2];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
#endif
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- gravity sort_gather/scatter functions --- */
extern "C" void mtl_sort_gather_force(int head_idx, int num_octs, int idim)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gather_force];
    [enc setBuffer:s_nref        offset:0 atIndex:0]; // s_nref as scratch
    [enc setBuffer:s_f_grav      offset:0 atIndex:1];
    [enc setBuffer:s_swap_global offset:0 atIndex:2];
    [enc setBytes:&idim     length:sizeof(int) atIndex:3];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:4];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:5];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

extern "C" void mtl_sort_scatter_force(int head_idx, int num_octs, int idim)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_scatter_force];
    [enc setBuffer:s_f_grav      offset:0 atIndex:0];
    [enc setBuffer:s_nref        offset:0 atIndex:1];
    [enc setBytes:&idim     length:sizeof(int) atIndex:2];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

extern "C" void mtl_sort_gather_phi(int head_idx, int num_octs, int is_phi_old)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gather_phi];
    [enc setBuffer:s_nref        offset:0 atIndex:0];
    [enc setBuffer:(is_phi_old ? s_phi_old : s_phi) offset:0 atIndex:1];
    [enc setBuffer:s_swap_global offset:0 atIndex:2];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

extern "C" void mtl_sort_scatter_phi(int head_idx, int num_octs, int is_phi_old)
{
    if (num_octs <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_octs + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_scatter_phi];
    [enc setBuffer:(is_phi_old ? s_phi_old : s_phi) offset:0 atIndex:0];
    [enc setBuffer:s_nref        offset:0 atIndex:1];
    [enc setBytes:&head_idx length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- mtl_blit_unew_to_uold: copy unew[head..] → uold (mirrors cudaMemcpy) */

extern "C" void mtl_blit_unew_to_uold(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger oct_floats = (NSUInteger)s_nvar * 8;
    NSUInteger src_off    = (NSUInteger)(head_idx - 1) * oct_floats * sizeof(float);
    NSUInteger blit_len   = (NSUInteger)num_octs * oct_floats * sizeof(float);
    id<MTLCommandBuffer>      cmd  = [s_queue commandBuffer];
    id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
    [blit copyFromBuffer:s_unew sourceOffset:src_off
               toBuffer:s_uold destinationOffset:src_off size:blit_len];
    [blit endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
#ifdef MHD
    NSUInteger boct = 6 * 8 * sizeof(float);
    NSUInteger boff = (NSUInteger)(head_idx - 1) * boct;
    id<MTLCommandBuffer> bcmd = [s_queue commandBuffer];
    id<MTLBlitCommandEncoder> bblit = [bcmd blitCommandEncoder];
    [bblit copyFromBuffer:s_bnew sourceOffset:boff toBuffer:s_bold destinationOffset:boff size:(NSUInteger)num_octs * boct];
    [bblit endEncoding];
    [bcmd commit];
    [bcmd waitUntilCompleted];
#endif
}

/* --- update_nbor_prefix_kernel: build one nbor column + prefix_sum ----- */

extern "C" void mtl_update_nbor_prefix(int head_idx, int num_subgrids,
                                       int hash_size, int input_ind)
{
    if (num_subgrids <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_subgrids + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_nbor_prefix];
    [enc setBuffer:s_nbor             offset:0 atIndex:0];
    [enc setBuffer:s_grid             offset:0 atIndex:1];
    [enc setBuffer:s_hash_key         offset:0 atIndex:2];
    [enc setBuffer:s_hash_val         offset:0 atIndex:3];
    [enc setBuffer:s_ckey_max_dev     offset:0 atIndex:4];
    [enc setBuffer:s_key_off_dev      offset:0 atIndex:5];
    [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:6];
    [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:7];
    [enc setBuffer:s_periodic_dev     offset:0 atIndex:8];
    [enc setBuffer:s_prefix_sum       offset:0 atIndex:9];
    [enc setBytes:&hash_size     length:sizeof(int) atIndex:10];
    [enc setBytes:&head_idx      length:sizeof(int) atIndex:11];
    [enc setBytes:&num_subgrids  length:sizeof(int) atIndex:12];
    [enc setBytes:&input_ind     length:sizeof(int) atIndex:13];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- compute_cache_swap_table_kernel: select subgrids missing a neighbour */

extern "C" void mtl_compute_cache_swap(int head_idx, int num_subgrids)
{
    if (num_subgrids <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)num_subgrids + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_cache_swap];
    [enc setBuffer:s_swap_local  offset:0 atIndex:0];
    [enc setBuffer:s_prefix_sum  offset:0 atIndex:1];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:2];
    [enc setBytes:&num_subgrids length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- make_cache_octs_kernel: create ghost octs in the cache region ----- */

extern "C" void mtl_make_cache_octs(int head_idx, int num_subgrids,
                                     int hash_size, int ngridmax,
                                     int ifree_cache, int new_noct,
                                     int input_ind, int interpol_var,
                                     int interpol_type, float smallr)
{
    if (new_noct <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)new_noct + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_make_cache];
    [enc setBuffer:s_grid             offset:0 atIndex:0];
    [enc setBuffer:s_flag1            offset:0 atIndex:1];
    [enc setBuffer:s_uold             offset:0 atIndex:2];
#ifdef MHD
    [enc setBuffer:s_bold             offset:0 atIndex:3];
    [enc setBuffer:s_swap_local       offset:0 atIndex:4];
    [enc setBuffer:s_father           offset:0 atIndex:5];
    [enc setBuffer:s_nbor             offset:0 atIndex:6];
    [enc setBuffer:s_hash_key         offset:0 atIndex:7];
    [enc setBuffer:s_hash_val         offset:0 atIndex:8];
    [enc setBuffer:s_ckey_max_dev     offset:0 atIndex:9];
    [enc setBuffer:s_key_off_dev      offset:0 atIndex:10];
    [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:11];
    [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:12];
    [enc setBuffer:s_periodic_dev     offset:0 atIndex:13];
    [enc setBytes:&hash_size   length:sizeof(int) atIndex:14];
    [enc setBytes:&ngridmax    length:sizeof(int) atIndex:15];
    [enc setBytes:&ifree_cache length:sizeof(int) atIndex:16];
    [enc setBytes:&new_noct    length:sizeof(int) atIndex:17];
    [enc setBytes:&input_ind   length:sizeof(int) atIndex:18];
    [enc setBytes:&interpol_var length:sizeof(int) atIndex:19];
    [enc setBytes:&interpol_type length:sizeof(int) atIndex:20];
    [enc setBytes:&smallr      length:sizeof(float) atIndex:21];
#ifdef GRAV
    [enc setBuffer:s_f_grav offset:0 atIndex:22];
    [enc setBuffer:s_phi offset:0 atIndex:23];
    [enc setBuffer:s_phi_old offset:0 atIndex:24];
#endif
#else
    [enc setBuffer:s_swap_local       offset:0 atIndex:3];
    [enc setBuffer:s_father           offset:0 atIndex:4];
    [enc setBuffer:s_nbor             offset:0 atIndex:5];
    [enc setBuffer:s_hash_key         offset:0 atIndex:6];
    [enc setBuffer:s_hash_val         offset:0 atIndex:7];
    [enc setBuffer:s_ckey_max_dev     offset:0 atIndex:8];
    [enc setBuffer:s_key_off_dev      offset:0 atIndex:9];
    [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:10];
    [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:11];
    [enc setBuffer:s_periodic_dev     offset:0 atIndex:12];
    [enc setBytes:&hash_size   length:sizeof(int) atIndex:13];
    [enc setBytes:&ngridmax    length:sizeof(int) atIndex:14];
    [enc setBytes:&ifree_cache length:sizeof(int) atIndex:15];
    [enc setBytes:&new_noct    length:sizeof(int) atIndex:16];
    [enc setBytes:&input_ind   length:sizeof(int) atIndex:17];
#ifdef GRAV
    [enc setBuffer:s_f_grav    offset:0 atIndex:18];
    [enc setBuffer:s_phi       offset:0 atIndex:19];
    [enc setBuffer:s_phi_old   offset:0 atIndex:20];
#endif
#endif
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- mtl_nbor_scan: update_nbor_prefix + all scan phases in one cmd buf -
 * Replaces separate mtl_update_nbor_prefix + mtl_prefix_scan + mtl_get_prefix_total.
 * Returns the inclusive-scan total (= number of subgrids missing neighbour ind).
 * Uses s_sort_fence for inter-encoder ordering within the command buffer.  */
extern "C" int mtl_nbor_scan(int head_idx, int num_subgrids,
                              int hash_size, int input_ind)
{
    if (num_subgrids <= 0) return 0;

    int offset    = head_idx - 1;   /* 0-based for scan */
    int n         = num_subgrids;
    int nblk_nbor = (n + 127) / 128;
    int nblk_scan = (n + 255) / 256;
    bool need_fixup = (nblk_scan > 1);

    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];

    /* 1. update_nbor_prefix_kernel — first encoder, no waitForFence */
    {
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_nbor_prefix];
        [enc setBuffer:s_nbor             offset:0 atIndex:0];
        [enc setBuffer:s_grid             offset:0 atIndex:1];
        [enc setBuffer:s_hash_key         offset:0 atIndex:2];
        [enc setBuffer:s_hash_val         offset:0 atIndex:3];
        [enc setBuffer:s_ckey_max_dev     offset:0 atIndex:4];
        [enc setBuffer:s_key_off_dev      offset:0 atIndex:5];
        [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:6];
        [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:7];
        [enc setBuffer:s_periodic_dev     offset:0 atIndex:8];
        [enc setBuffer:s_prefix_sum       offset:0 atIndex:9];
        [enc setBytes:&hash_size    length:sizeof(int) atIndex:10];
        [enc setBytes:&head_idx     length:sizeof(int) atIndex:11];
        [enc setBytes:&num_subgrids length:sizeof(int) atIndex:12];
        [enc setBytes:&input_ind    length:sizeof(int) atIndex:13];
        [enc dispatchThreadgroups:{(NSUInteger)nblk_nbor,1,1}
             threadsPerThreadgroup:{128,1,1}];
        [enc updateFence:s_sort_fence];
        [enc endEncoding];
    }
    /* 2. prefix scan: phase 1 */
    scan_phase_fenced(s_pso_scan_block, s_prefix_sum, s_partial_sums,
                      offset, n, cmd, s_sort_fence, s_sort_fence);
    if (need_fixup) {
        int nblk_scan2 = (nblk_scan + 255) / 256;
        scan_phase_fenced(s_pso_scan_block, s_partial_sums, s_partial_sums_2,
                          0, nblk_scan, cmd, s_sort_fence, s_sort_fence);
        if (nblk_scan2 > 1) {
            scan_phase_fenced(s_pso_scan_block, s_partial_sums_2, s_partial_sums_3,
                              0, nblk_scan2, cmd, s_sort_fence, s_sort_fence);
            scan_phase_fenced(s_pso_scan_fixup, s_partial_sums, s_partial_sums_2,
                              0, nblk_scan, cmd, s_sort_fence, s_sort_fence);
        }
        scan_phase_fenced(s_pso_scan_fixup, s_prefix_sum, s_partial_sums,
                          offset, n, cmd, s_sort_fence, s_sort_fence);
    }

    [cmd commit];
    [cmd waitUntilCompleted];

    /* Read total from shared memory — GPU done, no sync needed. */
    return ((int *)s_prefix_sum.contents)[offset + n - 1];
}

/* --- mtl_cache_fill: compute_cache_swap + make_cache_octs + insert_hash -
 * Replaces three separate commit/wait calls with one.
 * Uses s_sort_fence for inter-encoder ordering.                           */
extern "C" void mtl_cache_fill(int head_idx, int num_subgrids,
                                int hash_size, int input_ind,
                                int ngridmax, int ifree_cache, int new_noct,
                                int interpol_var, int interpol_type, float smallr)
{
    if (new_noct <= 0) return;

    int nblk_sg    = (num_subgrids + 127) / 128;
    int nblk_cache = (new_noct     + 127) / 128;

    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];

    /* 1. compute_cache_swap_table_kernel */
    {
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_cache_swap];
        [enc setBuffer:s_swap_local  offset:0 atIndex:0];
        [enc setBuffer:s_prefix_sum  offset:0 atIndex:1];
        [enc setBytes:&head_idx     length:sizeof(int) atIndex:2];
        [enc setBytes:&num_subgrids length:sizeof(int) atIndex:3];
        [enc dispatchThreadgroups:{(NSUInteger)nblk_sg,1,1}
             threadsPerThreadgroup:{128,1,1}];
        [enc updateFence:s_sort_fence];
        [enc endEncoding];
    }
    /* 2. make_cache_octs_kernel */
    {
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc waitForFence:s_sort_fence];
        [enc setComputePipelineState:s_pso_make_cache];
        [enc setBuffer:s_grid             offset:0 atIndex:0];
        [enc setBuffer:s_flag1            offset:0 atIndex:1];
        [enc setBuffer:s_uold             offset:0 atIndex:2];
#ifdef MHD
        [enc setBuffer:s_bold             offset:0 atIndex:3];
        [enc setBuffer:s_swap_local       offset:0 atIndex:4];
        [enc setBuffer:s_father           offset:0 atIndex:5];
        [enc setBuffer:s_nbor             offset:0 atIndex:6];
        [enc setBuffer:s_hash_key         offset:0 atIndex:7];
        [enc setBuffer:s_hash_val         offset:0 atIndex:8];
        [enc setBuffer:s_ckey_max_dev     offset:0 atIndex:9];
        [enc setBuffer:s_key_off_dev      offset:0 atIndex:10];
        [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:11];
        [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:12];
        [enc setBuffer:s_periodic_dev     offset:0 atIndex:13];
        [enc setBytes:&hash_size    length:sizeof(int) atIndex:14];
        [enc setBytes:&ngridmax     length:sizeof(int) atIndex:15];
        [enc setBytes:&ifree_cache  length:sizeof(int) atIndex:16];
        [enc setBytes:&new_noct     length:sizeof(int) atIndex:17];
        [enc setBytes:&input_ind    length:sizeof(int) atIndex:18];
        [enc setBytes:&interpol_var length:sizeof(int) atIndex:19];
        [enc setBytes:&interpol_type length:sizeof(int) atIndex:20];
        [enc setBytes:&smallr       length:sizeof(float) atIndex:21];
#ifdef GRAV
        [enc setBuffer:s_f_grav     offset:0 atIndex:22];
        [enc setBuffer:s_phi        offset:0 atIndex:23];
        [enc setBuffer:s_phi_old    offset:0 atIndex:24];
#endif
#else
        [enc setBuffer:s_swap_local       offset:0 atIndex:3];
        [enc setBuffer:s_father           offset:0 atIndex:4];
        [enc setBuffer:s_nbor             offset:0 atIndex:5];
        [enc setBuffer:s_hash_key         offset:0 atIndex:6];
        [enc setBuffer:s_hash_val         offset:0 atIndex:7];
        [enc setBuffer:s_ckey_max_dev     offset:0 atIndex:8];
        [enc setBuffer:s_key_off_dev      offset:0 atIndex:9];
        [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:10];
        [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:11];
        [enc setBuffer:s_periodic_dev     offset:0 atIndex:12];
        [enc setBytes:&hash_size    length:sizeof(int) atIndex:13];
        [enc setBytes:&ngridmax     length:sizeof(int) atIndex:14];
        [enc setBytes:&ifree_cache  length:sizeof(int) atIndex:15];
        [enc setBytes:&new_noct     length:sizeof(int) atIndex:16];
        [enc setBytes:&input_ind    length:sizeof(int) atIndex:17];
#ifdef GRAV
        [enc setBuffer:s_f_grav     offset:0 atIndex:18];
        [enc setBuffer:s_phi        offset:0 atIndex:19];
        [enc setBuffer:s_phi_old    offset:0 atIndex:20];
#endif
#endif
        [enc dispatchThreadgroups:{(NSUInteger)nblk_cache,1,1}
             threadsPerThreadgroup:{128,1,1}];
        [enc updateFence:s_sort_fence];
        [enc endEncoding];
    }
    /* 3. insert_hash_cache_kernel */
    {
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc waitForFence:s_sort_fence];
        [enc setComputePipelineState:s_pso_insert_hash_cache];
        [enc setBuffer:s_grid         offset:0 atIndex:0];
        [enc setBuffer:s_hash_key     offset:0 atIndex:1];
        [enc setBuffer:s_hash_val     offset:0 atIndex:2];
        [enc setBuffer:s_ckey_max_dev offset:0 atIndex:3];
        [enc setBuffer:s_key_off_dev  offset:0 atIndex:4];
        [enc setBytes:&hash_size    length:sizeof(int) atIndex:5];
        [enc setBytes:&ngridmax     length:sizeof(int) atIndex:6];
        [enc setBytes:&ifree_cache  length:sizeof(int) atIndex:7];
        [enc setBytes:&new_noct     length:sizeof(int) atIndex:8];
        [enc dispatchThreadgroups:{(NSUInteger)nblk_cache,1,1}
             threadsPerThreadgroup:{128,1,1}];
        [enc endEncoding];   /* last encoder — no outgoing fence needed */
    }

    [cmd commit];
    [cmd waitUntilCompleted];
}

/* --- insert_hash_cache_kernel: insert new cache octs into hash table --- */

extern "C" void mtl_insert_hash_cache_r(int hash_size, int ngridmax,
                                         int ifree_cache, int new_noct)
{
    if (new_noct <= 0) return;
    NSUInteger tg = 128, nblk = ((NSUInteger)new_noct + tg - 1) / tg;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_insert_hash_cache];
    [enc setBuffer:s_grid         offset:0 atIndex:0];
    [enc setBuffer:s_hash_key     offset:0 atIndex:1];
    [enc setBuffer:s_hash_val     offset:0 atIndex:2];
    [enc setBuffer:s_ckey_max_dev offset:0 atIndex:3];
    [enc setBuffer:s_key_off_dev  offset:0 atIndex:4];
    [enc setBytes:&hash_size   length:sizeof(int) atIndex:5];
    [enc setBytes:&ngridmax    length:sizeof(int) atIndex:6];
    [enc setBytes:&ifree_cache length:sizeof(int) atIndex:7];
    [enc setBytes:&new_noct    length:sizeof(int) atIndex:8];
    [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* =======================================================================
 * Gravity / Poisson bridge functions
 * ======================================================================= */

/* -----------------------------------------------------------------------
 * mtl_alloc_grav — allocate all gravity AMR + MG device buffers.
 * Called once from metal_allocate_grav in metal_runner.f90.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_alloc_grav(int ngridmax, int ncachemax,
                                int ngridmax_mg, int ncachemax_mg,
                                int hash_size_mg_in)
{
    s_ngridmax_mg  = ngridmax_mg;
    s_ncachemax_mg = ncachemax_mg;
    s_hash_size_mg = hash_size_mg_in;

    int ntotal    = ngridmax + ncachemax;
    int ntotal_mg = ngridmax_mg + ncachemax_mg;   /* CUDA: grid/phi/f indexed up to ngridmax_mg+ncachemax_mg */

    NSUInteger rho_bytes       = (NSUInteger)ntotal    * 8 * sizeof(float);
    NSUInteger nref_bytes      = (NSUInteger)ntotal    * 8 * sizeof(int);
    NSUInteger phi_bytes       = (NSUInteger)ntotal    * 8 * sizeof(float);
    NSUInteger f_grav_bytes    = (NSUInteger)ntotal    * 24 * sizeof(float);  /* 8 cells × 3 dims */
    NSUInteger multipole_bytes = (NSUInteger)ngridmax  * 4 * sizeof(float);
    NSUInteger scalar_bytes    = 4 * sizeof(float);  /* multipole_tot writes [m,mx,my,mz] */

    NSUInteger grid_mg_bytes   = (NSUInteger)ntotal_mg    * sizeof(oct_t);        /* CUDA: ngridmax_mg+ncachemax_mg */
    NSUInteger phi_mg_bytes    = (NSUInteger)ntotal_mg    * 8 * sizeof(float);
    NSUInteger f_mg_bytes      = (NSUInteger)ntotal_mg    * 24 * sizeof(float);
    NSUInteger nbor_mg_bytes   = (NSUInteger)ngridmax_mg  * 27 * sizeof(int);     /* CUDA: ngridmax_mg only, no cache */
    NSUInteger father_mg_bytes = (NSUInteger)(ngridmax + ngridmax_mg + ncachemax_mg) * sizeof(int);
    NSUInteger hkey_mg_bytes   = (NSUInteger)hash_size_mg_in * sizeof(long);
    NSUInteger hval_mg_bytes   = (NSUInteger)hash_size_mg_in * sizeof(int);

    s_rho           = [s_device newBufferWithLength:rho_bytes       options:MTLResourceStorageModeShared];
    s_nref          = [s_device newBufferWithLength:nref_bytes      options:MTLResourceStorageModeShared];
    s_phi           = [s_device newBufferWithLength:phi_bytes       options:MTLResourceStorageModeShared];
    s_phi_old       = [s_device newBufferWithLength:phi_bytes       options:MTLResourceStorageModeShared];
    s_f_grav        = [s_device newBufferWithLength:f_grav_bytes    options:MTLResourceStorageModeShared];
    s_multipole_buf = [s_device newBufferWithLength:multipole_bytes options:MTLResourceStorageModeShared];
    s_scalar_buf    = [s_device newBufferWithLength:scalar_bytes    options:MTLResourceStorageModeShared];

    s_grid_mg       = [s_device newBufferWithLength:grid_mg_bytes   options:MTLResourceStorageModeShared];
    s_phi_mg        = [s_device newBufferWithLength:phi_mg_bytes    options:MTLResourceStorageModeShared];
    s_f_mg          = [s_device newBufferWithLength:f_mg_bytes      options:MTLResourceStorageModeShared];
    s_nbor_mg       = [s_device newBufferWithLength:nbor_mg_bytes   options:MTLResourceStorageModeShared];
    s_father_mg     = [s_device newBufferWithLength:father_mg_bytes options:MTLResourceStorageModeShared];
    s_hash_key_mg   = [s_device newBufferWithLength:hkey_mg_bytes   options:MTLResourceStorageModeShared];
    s_hash_val_mg   = [s_device newBufferWithLength:hval_mg_bytes   options:MTLResourceStorageModeShared];

    memset(s_hash_key_mg.contents, 0, hkey_mg_bytes);
    memset(s_hash_val_mg.contents, 0, hval_mg_bytes);
}

/* -----------------------------------------------------------------------
 * Upload / download helpers for gravity arrays.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_upload_rho(void *rho_host, int head_idx, int num_octs)
{
    size_t off = (size_t)(head_idx - 1) * 8 * sizeof(float);
    size_t nb  = (size_t)num_octs * 8 * sizeof(float);
    memcpy((char *)s_rho.contents + off, rho_host, nb);
}

/* -----------------------------------------------------------------------
 * mtl_run_scan — inclusive prefix scan on s_prefix_sum[head_idx-1 .. head_idx+num_octs-2].
 * Returns the inclusive total (last element after scan).
 * Reuses existing scan PSOs; no nbor_prefix marking step.
 * ----------------------------------------------------------------------- */
extern "C" int mtl_run_scan(int head_idx, int num_octs)
{
    if (num_octs <= 0) return 0;

    int offset    = head_idx - 1;   /* 0-based */
    int n         = num_octs;
    int nb0       = (n   + 255) / 256;
    int nb1       = (nb0 + 255) / 256;
    bool need_fixup = (nb0 > 1);

    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];

    /* Phase 1: within-block scan; first encoder — no waitForFence */
    {
        NSUInteger nblk = (NSUInteger)nb0;
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_scan_block];
        [enc setBuffer:s_prefix_sum   offset:0 atIndex:0];
        [enc setBuffer:s_partial_sums offset:0 atIndex:1];
        [enc setBytes:&offset length:sizeof(int) atIndex:2];
        [enc setBytes:&n      length:sizeof(int) atIndex:3];
        [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{256,1,1}];
        [enc updateFence:s_sort_fence];
        [enc endEncoding];
    }

    if (need_fixup) {
        /* Phase 2: scan the block totals */
        scan_phase_fenced(s_pso_scan_block, s_partial_sums, s_partial_sums_2,
                          0, nb0, cmd, s_sort_fence, s_sort_fence);
        if (nb1 > 1) {
            /* 3rd level needed */
            scan_phase_fenced(s_pso_scan_block, s_partial_sums_2, s_partial_sums_3,
                              0, nb1, cmd, s_sort_fence, s_sort_fence);
            scan_phase_fenced(s_pso_scan_fixup, s_partial_sums, s_partial_sums_2,
                              0, nb0, cmd, s_sort_fence, s_sort_fence);
        }
        /* Phase 3: add block offsets back */
        scan_phase_fenced(s_pso_scan_fixup, s_prefix_sum, s_partial_sums,
                          offset, n, cmd, s_sort_fence, s_sort_fence);
    }

    [cmd commit];
    [cmd waitUntilCompleted];

    return ((int *)s_prefix_sum.contents)[offset + n - 1];
}

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
 * rho / density kernels
 * ----------------------------------------------------------------------- */

extern "C" void mtl_reset_rho(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_rho];
    [enc setBuffer:s_rho      offset:0 atIndex:0];
    [enc setBuffer:s_nref     offset:0 atIndex:1];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:3];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_multipole_leaf(int head_idx, int num_octs, float scale)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_multipole_leaf];
    [enc setBuffer:s_unew        offset:0 atIndex:0];  /* output: unew holds monopole+dipoles */
    [enc setBuffer:s_uold        offset:0 atIndex:1];  /* input:  gas density ivar=1           */
    [enc setBuffer:s_grid        offset:0 atIndex:2];  /* input:  oct ckey/refined              */
    [enc setBytes:&head_idx      length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs      length:sizeof(int)   atIndex:4];
    [enc setBytes:&scale         length:sizeof(float) atIndex:5];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_multipole_upload(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_multipole_upload];
    [enc setBuffer:s_grid     offset:0 atIndex:0];
    [enc setBuffer:s_father   offset:0 atIndex:1];
    [enc setBuffer:s_unew     offset:0 atIndex:2];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:4];
    DISPATCH_1D_128(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_multipole_tot(int head_idx, int num_octs, float *tot4)
{
    if (num_octs <= 0) { tot4[0]=tot4[1]=tot4[2]=tot4[3]=0.0f; return; }

    /* scalar_buf holds 4 floats; pre-zero before dispatch */
    memset(s_scalar_buf.contents, 0, 4 * sizeof(float));

    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_multipole_tot];
    [enc setBuffer:s_unew       offset:0 atIndex:0];
    [enc setBuffer:s_scalar_buf offset:0 atIndex:1];
    [enc setBytes:&head_idx        length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs        length:sizeof(int) atIndex:3];
    DISPATCH_1D_256_OCT(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    memcpy(tot4, s_scalar_buf.contents, 4 * sizeof(float));
}

extern "C" void mtl_deposit_rho(int head_idx, int num_octs,
                                float dx, float vol_loc,
                                float m_refine, float mass_sph,
                                float var_cut_refine, int ivar_refine,
                                int ngridmax)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_deposit_rho];
    [enc setBuffer:s_grid          offset:0 atIndex:0];
    [enc setBuffer:s_uold          offset:0 atIndex:1];
    [enc setBuffer:s_unew          offset:0 atIndex:2];
    [enc setBuffer:s_rho           offset:0 atIndex:3];
    [enc setBuffer:s_nref          offset:0 atIndex:4];
    [enc setBuffer:s_nbor          offset:0 atIndex:5];
    [enc setBytes:&head_idx        length:sizeof(int)   atIndex:6];
    [enc setBytes:&num_octs        length:sizeof(int)   atIndex:7];
    [enc setBytes:&dx              length:sizeof(float) atIndex:8];
    [enc setBytes:&vol_loc         length:sizeof(float) atIndex:9];
    [enc setBytes:&m_refine        length:sizeof(float) atIndex:10];
    [enc setBytes:&mass_sph        length:sizeof(float) atIndex:11];
    [enc setBytes:&var_cut_refine  length:sizeof(float) atIndex:12];
    [enc setBytes:&ivar_refine     length:sizeof(int)   atIndex:13];
    [enc setBytes:&ngridmax        length:sizeof(int)   atIndex:14];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_init_phi(int head_idx, int num_octs, int ngridmax, float tfrac)
{
    if (num_octs <= 0) return;

    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_make_initial_phi];
    [enc setBuffer:s_grid     offset:0 atIndex:0];
    [enc setBuffer:s_father   offset:0 atIndex:1];
    [enc setBuffer:s_nbor     offset:0 atIndex:2];
    [enc setBuffer:s_phi      offset:0 atIndex:3];
    [enc setBuffer:s_phi_old  offset:0 atIndex:4];
    [enc setBytes:&tfrac      length:sizeof(float) atIndex:5];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:6];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:7];
    [enc setBytes:&ngridmax   length:sizeof(int)   atIndex:8];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_save_phi_old_fine(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_save_phi_old];
    [enc setBuffer:s_phi_old  offset:0 atIndex:0];
    [enc setBuffer:s_phi      offset:0 atIndex:1];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:3];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * gradient / force
 * ----------------------------------------------------------------------- */

extern "C" void mtl_gradient_phi_fine(int head_idx, int num_octs, float dx)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gradient_phi];
    [enc setBuffer:s_phi      offset:0 atIndex:0];
    [enc setBuffer:s_f_grav   offset:0 atIndex:1];
    [enc setBuffer:s_nbor     offset:0 atIndex:2];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:4];
    [enc setBytes:&dx         length:sizeof(float) atIndex:5];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_gradient_phi_mg(int head_idx, int num_octs, float dx)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gradient_phi];
    [enc setBuffer:s_phi_mg   offset:0 atIndex:0];
    [enc setBuffer:s_f_mg     offset:0 atIndex:1];
    [enc setBuffer:s_nbor_mg  offset:0 atIndex:2];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:4];
    [enc setBytes:&dx         length:sizeof(float) atIndex:5];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * Scalar reduction bridge helpers
 * ----------------------------------------------------------------------- */

extern "C" void mtl_cmp_epot(int head_idx, int num_octs, float *epot_out)
{
    *epot_out = 0.0f;
    if (num_octs <= 0) return;
    memset(s_scalar_buf.contents, 0, sizeof(float));
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_cmp_epot];
    [enc setBuffer:s_grid       offset:0 atIndex:0];
    [enc setBuffer:s_f_grav     offset:0 atIndex:1];
    [enc setBuffer:s_scalar_buf offset:0 atIndex:2];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:3];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:4];
    DISPATCH_1D_256_OCT(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    memcpy(epot_out, s_scalar_buf.contents, sizeof(float));
}

extern "C" void mtl_cmp_rhomax(int head_idx, int num_octs, float *rhomax_out)
{
    *rhomax_out = 0.0f;
    if (num_octs <= 0) return;
    memset(s_scalar_buf.contents, 0, sizeof(float));
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_cmp_rhomax];
    [enc setBuffer:s_rho        offset:0 atIndex:0];
    [enc setBuffer:s_scalar_buf offset:0 atIndex:1];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:3];
    DISPATCH_1D_256_OCT(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    memcpy(rhomax_out, s_scalar_buf.contents, sizeof(float));
}

/* -----------------------------------------------------------------------
 * Mask operations
 * ----------------------------------------------------------------------- */

extern "C" void mtl_reset_mask_fine(int head_idx, int num_octs, float mask_val)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_mask_mg];
    [enc setBuffer:s_f_grav   offset:0 atIndex:0];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:1];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:2];
    [enc setBytes:&mask_val   length:sizeof(float) atIndex:3];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_reset_mask_mg(int head_idx, int num_octs, float mask_val)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_mask_mg];
    [enc setBuffer:s_f_mg     offset:0 atIndex:0];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:1];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:2];
    [enc setBytes:&mask_val   length:sizeof(float) atIndex:3];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_restrict_mask_fine(int head_idx, int head_father, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_restrict_mask_mg];
    [enc setBuffer:s_grid       offset:0 atIndex:0];
    [enc setBuffer:s_father_mg  offset:0 atIndex:1];
    [enc setBuffer:s_f_grav     offset:0 atIndex:2];
    [enc setBuffer:s_f_mg       offset:0 atIndex:3];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:4];
    [enc setBytes:&head_father  length:sizeof(int) atIndex:5];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:6];
    DISPATCH_1D_128(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_restrict_mask_mg(int head_idx, int head_father, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_restrict_mask_mg];
    [enc setBuffer:s_grid_mg    offset:0 atIndex:0];
    [enc setBuffer:s_father_mg  offset:0 atIndex:1];
    [enc setBuffer:s_f_mg       offset:0 atIndex:2];
    [enc setBuffer:s_f_mg       offset:0 atIndex:3];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:4];
    [enc setBytes:&head_father  length:sizeof(int) atIndex:5];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:6];
    DISPATCH_1D_128(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_volume_to_mask_fine(int head_idx, int num_octs, float *mask_max_out)
{
    *mask_max_out = -1.0f;
    if (num_octs <= 0) return;
    memset(s_scalar_buf.contents, 0, sizeof(float));
    *(float *)s_scalar_buf.contents = -1.0f;   /* init to -1 for max reduction */
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_volume_to_mask];
    [enc setBuffer:s_f_grav     offset:0 atIndex:0];
    [enc setBuffer:s_scalar_buf offset:0 atIndex:1];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:3];
    DISPATCH_1D_256_OCT(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    memcpy(mask_max_out, s_scalar_buf.contents, sizeof(float));
}

extern "C" void mtl_volume_to_mask_mg(int head_idx, int num_octs, float *mask_max_out)
{
    *mask_max_out = -1.0f;
    if (num_octs <= 0) return;
    *(float *)s_scalar_buf.contents = -1.0f;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_volume_to_mask];
    [enc setBuffer:s_f_mg       offset:0 atIndex:0];
    [enc setBuffer:s_scalar_buf offset:0 atIndex:1];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:3];
    DISPATCH_1D_256_OCT(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    memcpy(mask_max_out, s_scalar_buf.contents, sizeof(float));
}

/* -----------------------------------------------------------------------
 * Residual, RHS, phi reset
 * ----------------------------------------------------------------------- */

extern "C" void mtl_cmp_residual_fine(int head_idx, int num_octs, int ngridmax,
                                       float fourpi, float offset_val, float oneoverdx2)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_cmp_residual];
    [enc setBuffer:s_phi      offset:0 atIndex:0];
    [enc setBuffer:s_f_grav   offset:0 atIndex:1];
    [enc setBuffer:s_nbor     offset:0 atIndex:2];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:4];
    [enc setBytes:&ngridmax   length:sizeof(int)   atIndex:5];
    [enc setBytes:&fourpi     length:sizeof(float) atIndex:6];
    [enc setBytes:&offset_val length:sizeof(float) atIndex:7];
    [enc setBytes:&oneoverdx2 length:sizeof(float) atIndex:8];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_cmp_residual_mg(int head_idx, int num_octs, int ngridmax_mg_loc,
                                     float fourpi, float offset_val, float oneoverdx2)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_cmp_residual];
    [enc setBuffer:s_phi_mg   offset:0 atIndex:0];
    [enc setBuffer:s_f_mg     offset:0 atIndex:1];
    [enc setBuffer:s_nbor_mg  offset:0 atIndex:2];
    [enc setBytes:&head_idx       length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs       length:sizeof(int)   atIndex:4];
    [enc setBytes:&ngridmax_mg_loc length:sizeof(int)  atIndex:5];
    [enc setBytes:&fourpi         length:sizeof(float) atIndex:6];
    [enc setBytes:&offset_val     length:sizeof(float) atIndex:7];
    [enc setBytes:&oneoverdx2     length:sizeof(float) atIndex:8];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_gauss_seidel_fine(int head_idx, int num_octs, int ngridmax,
                                       float dx2, int safe, int redstep)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gauss_seidel];
    [enc setBuffer:s_phi      offset:0 atIndex:0];
    [enc setBuffer:s_f_grav   offset:0 atIndex:1];
    [enc setBuffer:s_nbor     offset:0 atIndex:2];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:4];
    [enc setBytes:&ngridmax   length:sizeof(int)   atIndex:5];
    [enc setBytes:&dx2        length:sizeof(float) atIndex:6];
    [enc setBytes:&safe       length:sizeof(int)   atIndex:7];
    [enc setBytes:&redstep    length:sizeof(int)   atIndex:8];
    DISPATCH_2D_4_32(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_gauss_seidel_mg(int head_idx, int num_octs, int ngridmax_mg_loc,
                                     float dx2, int safe, int redstep)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gauss_seidel];
    [enc setBuffer:s_phi_mg   offset:0 atIndex:0];
    [enc setBuffer:s_f_mg     offset:0 atIndex:1];
    [enc setBuffer:s_nbor_mg  offset:0 atIndex:2];
    [enc setBytes:&head_idx       length:sizeof(int)   atIndex:3];
    [enc setBytes:&num_octs       length:sizeof(int)   atIndex:4];
    [enc setBytes:&ngridmax_mg_loc length:sizeof(int)  atIndex:5];
    [enc setBytes:&dx2            length:sizeof(float) atIndex:6];
    [enc setBytes:&safe           length:sizeof(int)   atIndex:7];
    [enc setBytes:&redstep        length:sizeof(int)   atIndex:8];
    DISPATCH_2D_4_32(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_reset_phi_fine(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_phi_mg];
    [enc setBuffer:s_phi      offset:0 atIndex:0];
    [enc setBuffer:s_f_grav   offset:0 atIndex:1];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:3];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_reset_phi_mg(int head_idx, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_phi_mg];
    [enc setBuffer:s_phi_mg   offset:0 atIndex:0];
    [enc setBuffer:s_f_mg     offset:0 atIndex:1];
    [enc setBytes:&head_idx   length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs   length:sizeof(int) atIndex:3];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_reset_phi_val_fine(int head_idx, int num_octs, float phi_val)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_phi_val];
    [enc setBuffer:s_phi      offset:0 atIndex:0];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:1];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:2];
    [enc setBytes:&phi_val    length:sizeof(float) atIndex:3];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_reset_phi_val_mg(int head_idx, int num_octs, float phi_val)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_phi_val];
    [enc setBuffer:s_phi_mg   offset:0 atIndex:0];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:1];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:2];
    [enc setBytes:&phi_val    length:sizeof(float) atIndex:3];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_reset_rhs_fine(int head_idx, int num_octs, int ngridmax,
                                    float fourpi, float offset_val, float oneoverdx2)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_rhs_mg];
    [enc setBuffer:s_phi      offset:0 atIndex:0];
    [enc setBuffer:s_rho      offset:0 atIndex:1];
    [enc setBuffer:s_f_grav   offset:0 atIndex:2];
    [enc setBuffer:s_nbor     offset:0 atIndex:3];
    [enc setBytes:&head_idx   length:sizeof(int)   atIndex:4];
    [enc setBytes:&num_octs   length:sizeof(int)   atIndex:5];
    [enc setBytes:&ngridmax   length:sizeof(int)   atIndex:6];
    [enc setBytes:&fourpi     length:sizeof(float) atIndex:7];
    [enc setBytes:&offset_val length:sizeof(float) atIndex:8];
    [enc setBytes:&oneoverdx2 length:sizeof(float) atIndex:9];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_reset_rhs_mg(int head_idx, int num_octs, int ngridmax_mg_loc,
                                  float fourpi, float offset_val, float oneoverdx2)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_reset_rhs_mg];
    [enc setBuffer:s_phi_mg   offset:0 atIndex:0];
    [enc setBuffer:s_f_mg     offset:0 atIndex:1];   /* rho slot reuses f2 in MG path */
    [enc setBuffer:s_f_mg     offset:0 atIndex:2];
    [enc setBuffer:s_nbor_mg  offset:0 atIndex:3];
    [enc setBytes:&head_idx       length:sizeof(int)   atIndex:4];
    [enc setBytes:&num_octs       length:sizeof(int)   atIndex:5];
    [enc setBytes:&ngridmax_mg_loc length:sizeof(int)  atIndex:6];
    [enc setBytes:&fourpi         length:sizeof(float) atIndex:7];
    [enc setBytes:&offset_val     length:sizeof(float) atIndex:8];
    [enc setBytes:&oneoverdx2     length:sizeof(float) atIndex:9];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * Residual norm
 * ----------------------------------------------------------------------- */

extern "C" void mtl_residual_norm_fine(int head_idx, int num_octs, float *norm_out)
{
    *norm_out = 0.0f;
    if (num_octs <= 0) return;
    memset(s_scalar_buf.contents, 0, sizeof(float));
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_residual_norm];
    [enc setBuffer:s_f_grav     offset:0 atIndex:0];
    [enc setBuffer:s_scalar_buf offset:0 atIndex:1];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:3];
    DISPATCH_1D_256_OCT(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    memcpy(norm_out, s_scalar_buf.contents, sizeof(float));
}

extern "C" void mtl_rhs_norm_fine(int head_idx, int num_octs, float *norm_out)
{
    *norm_out = 0.0f;
    if (num_octs <= 0) return;
    memset(s_scalar_buf.contents, 0, sizeof(float));
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_rhs_norm];
    [enc setBuffer:s_f_grav     offset:0 atIndex:0];
    [enc setBuffer:s_scalar_buf offset:0 atIndex:1];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:3];
    DISPATCH_1D_256_OCT(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    memcpy(norm_out, s_scalar_buf.contents, sizeof(float));
}

extern "C" void mtl_residual_norm_mg(int head_idx, int num_octs, float *norm_out)
{
    *norm_out = 0.0f;
    if (num_octs <= 0) return;
    memset(s_scalar_buf.contents, 0, sizeof(float));
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_residual_norm];
    [enc setBuffer:s_f_mg       offset:0 atIndex:0];
    [enc setBuffer:s_scalar_buf offset:0 atIndex:1];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:2];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:3];
    DISPATCH_1D_256_OCT(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    memcpy(norm_out, s_scalar_buf.contents, sizeof(float));
}

/* -----------------------------------------------------------------------
 * Restrict residual
 * ----------------------------------------------------------------------- */

extern "C" void mtl_restrict_residual_fine(int head_idx, int head_father, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_restrict_residual];
    [enc setBuffer:s_grid       offset:0 atIndex:0];
    [enc setBuffer:s_father_mg  offset:0 atIndex:1];
    [enc setBuffer:s_f_grav     offset:0 atIndex:2];
    [enc setBuffer:s_f_mg       offset:0 atIndex:3];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:4];
    [enc setBytes:&head_father  length:sizeof(int) atIndex:5];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:6];
    DISPATCH_1D_128(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_restrict_residual_mg(int head_idx, int head_father, int num_octs)
{
    if (num_octs <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_restrict_residual];
    [enc setBuffer:s_grid_mg    offset:0 atIndex:0];
    [enc setBuffer:s_father_mg  offset:0 atIndex:1];
    [enc setBuffer:s_f_mg       offset:0 atIndex:2];
    [enc setBuffer:s_f_mg       offset:0 atIndex:3];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:4];
    [enc setBytes:&head_father  length:sizeof(int) atIndex:5];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:6];
    DISPATCH_1D_128(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

/* -------------------------------------------------------------------------
 * mtl_debug_rho — print min/max/mean/stddev of s_rho for cells belonging
 * to grids [head_grid .. head_grid+num_octs-1] (1-based, Fortran convention).
 * twotondim cells per grid (typically 8 for NDIM=3).
 * Called from Fortran: call mtl_debug_rho(head, noct, twotondim, ilevel)
 * ------------------------------------------------------------------------- */
extern "C" void mtl_debug_rho(int head_grid, int num_octs, int twotondim, int ilevel)
{
    if (num_octs <= 0 || s_rho == nil) {
        printf("[MTL] mtl_debug_rho: no grids at level %d\n", ilevel);
        return;
    }

    // s_rho is MTLResourceStorageModeShared — readable directly on the CPU.
    const float *rho = (const float *)s_rho.contents;

    int    n       = num_octs * twotondim;
    int    base    = (head_grid - 1) * twotondim;  // 0-based flat index of first cell

    double sum     = 0.0;
    double sum_sq  = 0.0;
    double rho_min = (double)rho[base];
    double rho_max = (double)rho[base];
    int    n_nonzero = 0;

    for (int i = 0; i < n; i++) {
        double v = (double)rho[base + i];
        sum    += v;
        sum_sq += v * v;
        if (v < rho_min) rho_min = v;
        if (v > rho_max) rho_max = v;
        if (v != 0.0) n_nonzero++;
    }

    double mean   = sum / n;
    double var    = sum_sq / n - mean * mean;
    double stddev = (var > 0.0) ? sqrt(var) : 0.0;

    printf("[MTL] rho stats level %d: ncells=%d nonzero=%d  min=%e  max=%e  mean=%e  stddev=%e\n",
           ilevel, n, n_nonzero, rho_min, rho_max, mean, stddev);
}

/* -------------------------------------------------------------------------
 * mtl_debug_xp — print min/max/mean of s_xp for particles [head_idx .. head_idx+num_parts-1]
 * (1-based, Fortran convention) for each spatial dimension.
 * Called from Fortran: call mtl_debug_xp(head_idx, num_parts, ilevel)
 * ------------------------------------------------------------------------- */
extern "C" void mtl_debug_xp(int head_idx, int num_parts, int ilevel)
{
    if (num_parts <= 0 || s_xp == nil) {
        printf("[MTL] mtl_debug_xp: no particles at level %d\n", ilevel);
        return;
    }

    // s_xp is MTLResourceStorageModeShared — readable directly on the CPU.
    const float *xp = (const float *)s_xp.contents;
    long leading    = (long)s_npartmax;
    int  base       = head_idx - 1;  // 0-based index of first particle

    const char *dim_names[3] = {"x", "y", "z"};
    for (int idim = 0; idim < 3; idim++) {
        long off = (long)idim * leading;
        double xmin = (double)xp[base + off];
        double xmax = xmin;
        double xsum = 0.0;
        for (int i = 0; i < num_parts; i++) {
            double v = (double)xp[base + i + off];
            if (v < xmin) xmin = v;
            if (v > xmax) xmax = v;
            xsum += v;
        }
        double xmean = xsum / num_parts;
        printf("[MTL] xp(%s) stats level %d: npart=%d  min=%e  max=%e  mean=%e\n",
               dim_names[idim], ilevel, num_parts, xmin, xmax, xmean);
    }
}

/* -----------------------------------------------------------------------
 * Interpolate / correct
 * ----------------------------------------------------------------------- */

extern "C" void mtl_interpolate_correct_fine(int head_idx, int head_father, int num_octs)
{
    if (num_octs <= 0) return;
    /* threadgroup memory: 27*16 ints per array, two arrays */
    NSUInteger tg_bytes = 2 * 27 * 16 * sizeof(int);
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_interpolate_correct];
    [enc setBuffer:s_grid       offset:0 atIndex:0];
    [enc setBuffer:s_father_mg  offset:0 atIndex:1];
    [enc setBuffer:s_nbor_mg    offset:0 atIndex:2];
    [enc setBuffer:s_phi        offset:0 atIndex:3];
    [enc setBuffer:s_phi_mg     offset:0 atIndex:4];
    [enc setBuffer:s_f_grav     offset:0 atIndex:5];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:6];
    [enc setBytes:&head_father  length:sizeof(int) atIndex:7];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:8];
    [enc setThreadgroupMemoryLength:tg_bytes / 2 atIndex:0];   /* igrid_nbor */
    [enc setThreadgroupMemoryLength:tg_bytes / 2 atIndex:1];   /* icell_nbor */
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_interpolate_correct_mg(int head_idx, int head_father, int num_octs)
{
    if (num_octs <= 0) return;
    NSUInteger tg_bytes = 2 * 27 * 16 * sizeof(int);
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_interpolate_correct];
    [enc setBuffer:s_grid_mg    offset:0 atIndex:0];
    [enc setBuffer:s_father_mg  offset:0 atIndex:1];
    [enc setBuffer:s_nbor_mg    offset:0 atIndex:2];
    [enc setBuffer:s_phi_mg     offset:0 atIndex:3];
    [enc setBuffer:s_phi_mg     offset:0 atIndex:4];
    [enc setBuffer:s_f_mg       offset:0 atIndex:5];
    [enc setBytes:&head_idx     length:sizeof(int) atIndex:6];
    [enc setBytes:&head_father  length:sizeof(int) atIndex:7];
    [enc setBytes:&num_octs     length:sizeof(int) atIndex:8];
    [enc setThreadgroupMemoryLength:tg_bytes / 2 atIndex:0];
    [enc setThreadgroupMemoryLength:tg_bytes / 2 atIndex:1];
    DISPATCH_2D_8_16(enc, num_octs);
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

/* -----------------------------------------------------------------------
 * mtl_build_mg_fine / mtl_build_mg_mg
 *
 * Build one MG level (AMR source or MG source):
 *   1. init_prefix_sum_mg   — mark first oct per parent group
 *   2. mtl_run_scan         — inclusive scan → new_noct (total distinct parents)
 *   3. compute_father_swap  — build swap table
 *   4. make_father_octs     — fill new MG octs
 *   5. insert_hash_mg       — insert new MG octs into hash table
 *   6. update_father_array  — link each source oct to its MG parent
 *   7. update_nbor_array_mg — build MG nbor table from hash
 *
 * ifine:       fine AMR level (Fortran 1-based)
 * ilevel:      coarse MG level = ifine - 1
 * head_idx:    first source oct (1-based)
 * num_octs:    number of source octs
 * head_father: first MG parent slot (1-based index into s_grid_mg / s_father_mg)
 * new_noct:    (out) number of new MG parent octs created
 * ----------------------------------------------------------------------- */

/* head_father: offset into s_father_mg for source octs (= CUDA's head_father)
 * head_mg:     starting index of new MG octs in s_grid_mg (= CUDA's head_mg)
 * For the fine case both equal 1; for MG-MG they differ. */
static void build_mg_common(id<MTLBuffer> grid_src, int head_idx, int num_octs,
                             int head_father, int head_mg, int *new_noct_out)
{
    /* Step 1: init prefix sum marks */
    {
        NSUInteger tg  = 128;
        NSUInteger nblk = ((NSUInteger)num_octs + tg - 1) / tg;
        id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_init_prefix_sum_mg];
        [enc setBuffer:grid_src      offset:0 atIndex:0];
        [enc setBuffer:s_prefix_sum  offset:0 atIndex:1];
        [enc setBytes:&head_idx      length:sizeof(int) atIndex:2];
        [enc setBytes:&num_octs      length:sizeof(int) atIndex:3];
        [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    }

    /* Step 2: prefix scan → new_noct */
    int new_noct = mtl_run_scan(head_idx, num_octs);
    *new_noct_out = new_noct;
    if (new_noct > 0) {

      if (head_mg + new_noct - 1 > s_ngridmax_mg) {
        fprintf(stderr, "No more grid memory, increase ngridmax for MG\n");
        fprintf(stderr, "New multigrid octs: %d, head_mg: %d, ngridmax_mg: %d\n", new_noct, head_mg, s_ngridmax_mg);
        exit(1);
      }

      /* Step 3: compute father swap table */
      {
        NSUInteger tg  = 128;
        NSUInteger nblk = ((NSUInteger)num_octs + tg - 1) / tg;
        id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_compute_father_swap];
        [enc setBuffer:s_swap_local  offset:0 atIndex:0];
        [enc setBuffer:s_prefix_sum  offset:0 atIndex:1];
        [enc setBytes:&head_idx      length:sizeof(int) atIndex:2];
        [enc setBytes:&num_octs      length:sizeof(int) atIndex:3];
        [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
      }

      /* Step 4: make father octs in MG grid — use head_mg (new octs position in s_grid_mg) */
      {
        NSUInteger tg  = 128;
        NSUInteger nblk = ((NSUInteger)new_noct + tg - 1) / tg;
        id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_make_father_octs];
        [enc setBuffer:grid_src      offset:0 atIndex:0];
        [enc setBuffer:s_grid_mg     offset:0 atIndex:1];
        [enc setBuffer:s_phi_mg      offset:0 atIndex:2];
        [enc setBuffer:s_f_mg        offset:0 atIndex:3];
        [enc setBuffer:s_swap_local  offset:0 atIndex:4];
        [enc setBytes:&head_mg       length:sizeof(int) atIndex:5];
        [enc setBytes:&new_noct      length:sizeof(int) atIndex:6];
        [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
      }

      /* Step 5: insert new MG octs into hash — use head_mg */
      {
        NSUInteger tg  = 128;
        NSUInteger nblk = ((NSUInteger)new_noct + tg - 1) / tg;
        id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_insert_hash_all];   /* reuse AMR insert_hash_all */
        [enc setBuffer:s_grid_mg        offset:0 atIndex:0];
        [enc setBuffer:s_hash_key_mg    offset:0 atIndex:1];
        [enc setBuffer:s_hash_val_mg    offset:0 atIndex:2];
        [enc setBuffer:s_ckey_max_dev   offset:0 atIndex:3];
        [enc setBuffer:s_key_off_dev    offset:0 atIndex:4];
        [enc setBytes:&s_hash_size_mg   length:sizeof(int) atIndex:5];
        [enc setBytes:&head_mg          length:sizeof(int) atIndex:6];
        [enc setBytes:&new_noct         length:sizeof(int) atIndex:7];
        [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
      }
    }

    /* Step 6: update father_mg for each source oct — use head_father (s_father_mg slot) */
    {
        NSUInteger tg  = 128;
        NSUInteger nblk = ((NSUInteger)num_octs + tg - 1) / tg;
        id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_update_father_array];
        [enc setBuffer:s_father_mg      offset:0 atIndex:0];
        [enc setBuffer:grid_src         offset:0 atIndex:1];
        [enc setBuffer:s_hash_key_mg    offset:0 atIndex:2];
        [enc setBuffer:s_hash_val_mg    offset:0 atIndex:3];
        [enc setBuffer:s_ckey_max_dev   offset:0 atIndex:4];
        [enc setBuffer:s_key_off_dev    offset:0 atIndex:5];
        [enc setBytes:&s_hash_size_mg   length:sizeof(int) atIndex:6];
        [enc setBytes:&head_idx         length:sizeof(int) atIndex:7];
        [enc setBytes:&head_father      length:sizeof(int) atIndex:8];
        [enc setBytes:&num_octs         length:sizeof(int) atIndex:9];
        [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
    }

    /* Step 7: build MG nbor array for the new MG octs — use head_mg */
    if (new_noct > 0) {
      {
        NSUInteger tg  = 128;
        NSUInteger nblk = ((NSUInteger)new_noct + tg - 1) / tg;
        id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_update_nbor_array_mg];
        [enc setBuffer:s_nbor_mg        offset:0 atIndex:0];
        [enc setBuffer:s_grid_mg        offset:0 atIndex:1];
        [enc setBuffer:s_hash_key_mg    offset:0 atIndex:2];
        [enc setBuffer:s_hash_val_mg    offset:0 atIndex:3];
        [enc setBuffer:s_ckey_max_dev   offset:0 atIndex:4];
        [enc setBuffer:s_key_off_dev    offset:0 atIndex:5];
        [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:6];
        [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:7];
        [enc setBuffer:s_periodic_dev   offset:0 atIndex:8];
        [enc setBytes:&s_hash_size_mg   length:sizeof(int) atIndex:9];
        [enc setBytes:&head_mg          length:sizeof(int) atIndex:10];
        [enc setBytes:&new_noct         length:sizeof(int) atIndex:11];
        [enc dispatchThreadgroups:{nblk,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
      }
    }
}

/* -----------------------------------------------------------------------
 * mtl_clean_mg_hashes — zero the MG hash table to remove stale tombstones.
 * Called by metal_clean_mg before each new multigrid cycle.
 * ----------------------------------------------------------------------- */
extern "C" void mtl_clean_mg_hashes(void)
{
    if (s_hash_key_mg) memset(s_hash_key_mg.contents, 0, s_hash_key_mg.length);
    if (s_hash_val_mg) memset(s_hash_val_mg.contents, 0, s_hash_val_mg.length);
}

extern "C" void mtl_build_mg_fine(int ifine, int ilevel,
                                   int head_idx, int num_octs,
                                   int head_father, int head_mg, int *new_noct)
{
    (void)ifine; (void)ilevel;
    build_mg_common(s_grid, head_idx, num_octs, head_father, head_mg, new_noct);
}

extern "C" void mtl_build_mg_mg(int ifine, int ilevel,
                                 int head_idx, int num_octs,
                                 int head_father, int head_mg, int *new_noct)
{
    (void)ifine; (void)ilevel;
    build_mg_common(s_grid_mg, head_idx, num_octs, head_father, head_mg, new_noct);
}

struct CoolingParams {
    int table_n1;
    int table_n2;
    int head_idx;
    int num_octs;
    float dlog_nH;
    float dlog_T2;
    float X_frac;
    float dt;
    float scale_T2;
    float scale_nH;
    float z_ave;
    float T2max;
    float gamma;
    float smallr;
    float smallc2;
    int eos_type;
    float eos_T2;
    float eos_nH;
    float eos_index;
    int imetal;
    int cooling;
    int metal;
    int self_shielding;
    int isothermal;
};

static void upload_array(__strong id<MTLBuffer>& buf, const double* src, NSUInteger size_bytes, NSUInteger count) {
    if (buf && buf.length != size_bytes) {
        buf = nil;
    }
    if (!buf) {
        buf = [s_device newBufferWithLength:size_bytes options:MTLResourceStorageModeShared];
    }
    float* dst = (float*)buf.contents;
    for (NSUInteger i = 0; i < count; i++) {
        dst[i] = (float)src[i];
    }
}

extern "C" void mtl_upload_cooling_table(
    int n1, int n2,
    const double* nH_tbl, const double* T2_tbl,
    const double* cool, const double* heat,
    const double* cool_com, const double* heat_com,
    const double* metal, const double* cool_prime,
    const double* heat_prime, const double* cool_com_prime,
    const double* heat_com_prime, const double* metal_prime
) {
    s_table_n1 = n1;
    s_table_n2 = n2;
    s_table_dlog_nH = (float)(n1 - 1) / (float)(nH_tbl[n1 - 1] - nH_tbl[0]);
    s_table_dlog_T2 = (float)(n2 - 1) / (float)(T2_tbl[n2 - 1] - T2_tbl[0]);

    NSUInteger axis1_bytes = n1 * sizeof(float);
    NSUInteger axis2_bytes = n2 * sizeof(float);
    NSUInteger table_bytes = (NSUInteger)n1 * n2 * sizeof(float);

    upload_array(s_nH_tbl_d, nH_tbl, axis1_bytes, n1);
    upload_array(s_T2_tbl_d, T2_tbl, axis2_bytes, n2);

    NSUInteger total_cells = (NSUInteger)n1 * n2;
    upload_array(s_cool_d, cool, table_bytes, total_cells);
    upload_array(s_heat_d, heat, table_bytes, total_cells);
    upload_array(s_cool_com_d, cool_com, table_bytes, total_cells);
    upload_array(s_heat_com_d, heat_com, table_bytes, total_cells);
    upload_array(s_metal_d, metal, table_bytes, total_cells);
    upload_array(s_cool_prime_d, cool_prime, table_bytes, total_cells);
    upload_array(s_heat_prime_d, heat_prime, table_bytes, total_cells);
    upload_array(s_cool_com_prime_d, cool_com_prime, table_bytes, total_cells);
    upload_array(s_heat_com_prime_d, heat_com_prime, table_bytes, total_cells);
    upload_array(s_metal_prime_d, metal_prime, table_bytes, total_cells);

    s_table_uploaded = YES;
}

extern "C" void mtl_cooling(
    int head_idx, int num_octs,
    float gamma, float smallr, float smallc2,
    double dtcool, int eos_type, double eos_T2,
    double eos_nH, double eos_index, double scale_T2, double scale_nH,
    int cooling, int metal, int imetal, double z_ave,
    int self_shielding, double X_frac, double T2max, int isothermal
) {
    if (num_octs <= 0) return;

    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_cooling];

    CoolingParams params;
    params.table_n1 = s_table_n1;
    params.table_n2 = s_table_n2;
    params.head_idx = head_idx - 1; // head_idx is 1-based, offset in grid is 0-based
    params.num_octs = num_octs;
    params.dlog_nH = s_table_dlog_nH;
    params.dlog_T2 = s_table_dlog_T2;
    params.X_frac = (float)X_frac;
    params.dt = (float)dtcool;
    params.scale_T2 = (float)scale_T2;
    params.scale_nH = (float)scale_nH;
    params.z_ave = (float)z_ave;
    params.T2max = (float)T2max;
    params.gamma = gamma;
    params.smallr = smallr;
    params.smallc2 = smallc2;
    params.eos_type = eos_type;
    params.eos_T2 = (float)eos_T2;
    params.eos_nH = (float)eos_nH;
    params.eos_index = (float)eos_index;
    params.imetal = imetal;
    params.cooling = cooling;
    params.metal = metal;
    params.self_shielding = self_shielding;
    params.isothermal = isothermal;

    [enc setBuffer:s_uold offset:0 atIndex:0];
    [enc setBuffer:nil offset:0 atIndex:1];
    [enc setBuffer:s_grid offset:0 atIndex:2];
    [enc setBytes:&params length:sizeof(CoolingParams) atIndex:3];

    if (cooling && s_table_uploaded) {
        [enc setBuffer:s_nH_tbl_d offset:0 atIndex:4];
        [enc setBuffer:s_T2_tbl_d offset:0 atIndex:5];
        [enc setBuffer:s_cool_d offset:0 atIndex:6];
        [enc setBuffer:s_heat_d offset:0 atIndex:7];
        [enc setBuffer:s_cool_com_d offset:0 atIndex:8];
        [enc setBuffer:s_heat_com_d offset:0 atIndex:9];
        [enc setBuffer:s_metal_d offset:0 atIndex:10];
        [enc setBuffer:s_cool_prime_d offset:0 atIndex:11];
        [enc setBuffer:s_heat_prime_d offset:0 atIndex:12];
        [enc setBuffer:s_cool_com_prime_d offset:0 atIndex:13];
        [enc setBuffer:s_heat_com_prime_d offset:0 atIndex:14];
        [enc setBuffer:s_metal_prime_d offset:0 atIndex:15];
    }

    NSUInteger tgWidth = 8;
    NSUInteger tgHeight = 16;
    MTLSize threadgroupsPerGrid = MTLSizeMake(1, (num_octs + tgHeight - 1) / tgHeight, 1);
    MTLSize threadsPerThreadgroup = MTLSizeMake(tgWidth, tgHeight, 1);

    [enc dispatchThreadgroups:threadgroupsPerGrid threadsPerThreadgroup:threadsPerThreadgroup];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];
}

/* =======================================================================
 * Particle bridge functions
 * ======================================================================= */

extern "C" void mtl_alloc_part(int npartmax)
{
    s_npartmax = npartmax;
    NSUInteger xp_bytes = (NSUInteger)npartmax * 3 * sizeof(float);
    NSUInteger vp_bytes = (NSUInteger)npartmax * 3 * sizeof(float);
    NSUInteger mp_bytes = (NSUInteger)npartmax * sizeof(float);
    NSUInteger idp_bytes = (NSUInteger)npartmax * sizeof(int);
    NSUInteger levelp_bytes = (NSUInteger)npartmax * sizeof(int);
    NSUInteger sortp_bytes = (NSUInteger)npartmax * sizeof(int);
    NSUInteger isp_swap_bytes = (NSUInteger)npartmax * sizeof(int);

    s_xp              = [s_device newBufferWithLength:xp_bytes        options:MTLResourceStorageModeShared];
    s_vp              = [s_device newBufferWithLength:vp_bytes        options:MTLResourceStorageModeShared];
    s_mp              = [s_device newBufferWithLength:mp_bytes        options:MTLResourceStorageModeShared];
    s_idp             = [s_device newBufferWithLength:idp_bytes       options:MTLResourceStorageModeShared];
    s_levelp          = [s_device newBufferWithLength:levelp_bytes    options:MTLResourceStorageModeShared];
    s_sortp           = [s_device newBufferWithLength:sortp_bytes     options:MTLResourceStorageModeShared];
    s_xp_swap         = [s_device newBufferWithLength:xp_bytes        options:MTLResourceStorageModeShared];
    s_isp_swap        = [s_device newBufferWithLength:isp_swap_bytes  options:MTLResourceStorageModeShared];
    s_idp_swap        = [s_device newBufferWithLength:(NSUInteger)npartmax * sizeof(long) options:MTLResourceStorageModeShared];
    s_prefix_sum_part = [s_device newBufferWithLength:sortp_bytes     options:MTLResourceStorageModeShared];
    s_part_scalar_buf = [s_device newBufferWithLength:8               options:MTLResourceStorageModeShared];
    s_multipole_q_part_buf = [s_device newBufferWithLength:4 * sizeof(float) options:MTLResourceStorageModeShared];

    // Allocate partial sum buffers for particles
    NSUInteger partial_bytes  = (NSUInteger)((npartmax + 255) / 256) * sizeof(int);
    NSUInteger partial2_bytes = (NSUInteger)((npartmax / 256 + 255) / 256) * sizeof(int);
    NSUInteger partial3_bytes = (NSUInteger)((npartmax / (256 * 256) + 255) / 256) * sizeof(int);
    s_partial_sums_part  = [s_device newBufferWithLength:partial_bytes  options:MTLResourceStorageModeShared];
    s_partial_sums_part2 = [s_device newBufferWithLength:partial2_bytes options:MTLResourceStorageModeShared];
    s_partial_sums_part3 = [s_device newBufferWithLength:partial3_bytes options:MTLResourceStorageModeShared];
    s_partial_sums_part4 = [s_device newBufferWithLength:sizeof(int)    options:MTLResourceStorageModeShared];
}

extern "C" void mtl_upload_part(void* xp, void* vp, void* mp, void* levelp, void* sortp, void* idp, int npart)
{
    if (npart <= 0) return;
    memcpy(s_xp.contents, xp, (size_t)s_npartmax * 3 * sizeof(float));
    memcpy(s_vp.contents, vp, (size_t)s_npartmax * 3 * sizeof(float));
    memcpy(s_mp.contents, mp, (size_t)s_npartmax * sizeof(float));
    memcpy(s_levelp.contents, levelp, (size_t)s_npartmax * sizeof(int));
    memcpy(s_sortp.contents, sortp, (size_t)s_npartmax * sizeof(int));
    if (idp) {
        memcpy(s_idp.contents, idp, (size_t)s_npartmax * sizeof(int));
    }
}

extern "C" void mtl_download_part(void* xp, void* vp, void* mp, void* levelp, void* sortp, void* idp, int npart)
{
    if (npart <= 0) return;
    memcpy(xp, s_xp.contents, (size_t)s_npartmax * 3 * sizeof(float));
    memcpy(vp, s_vp.contents, (size_t)s_npartmax * 3 * sizeof(float));
    memcpy(mp, s_mp.contents, (size_t)s_npartmax * sizeof(float));
    memcpy(levelp, s_levelp.contents, (size_t)s_npartmax * sizeof(int));
    memcpy(sortp, s_sortp.contents, (size_t)s_npartmax * sizeof(int));
    if (idp) {
        memcpy(idp, s_idp.contents, (size_t)s_npartmax * sizeof(int));
    }
}

/* mtl_multipole_q_part — GPU reduction of total mass + centre-of-mass moments.
 * Mirrors multipole_q_kernel (gpu_part.cuf) called at ilevel==levelmin.
 * Returns q_out[0..3] = { sum(mp), sum(mp*x), sum(mp*y), sum(mp*z) }.
 */
extern "C" void mtl_multipole_q_part(int head_idx, int num_parts, long leading, float *q_out)
{
    if (num_parts <= 0) { q_out[0]=q_out[1]=q_out[2]=q_out[3]=0.0f; return; }

    // Zero the accumulator buffer (shared memory, so host write is fine)
    memset(s_multipole_q_part_buf.contents, 0, 4 * sizeof(float));

    NSUInteger tg = 256;
    NSUInteger nb = ((NSUInteger)num_parts + tg - 1) / tg;

    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_multipole_q_part];
    [enc setBuffer:s_xp                   offset:0 atIndex:0];
    [enc setBuffer:s_mp                   offset:0 atIndex:1];
    [enc setBuffer:s_multipole_q_part_buf offset:0 atIndex:2];
    [enc setBytes:&head_idx               length:sizeof(int)  atIndex:3];
    [enc setBytes:&num_parts              length:sizeof(int)  atIndex:4];
    [enc setBytes:&leading                length:sizeof(long) atIndex:5];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    memcpy(q_out, s_multipole_q_part_buf.contents, 4 * sizeof(float));
}

extern "C" void mtl_kick_drift_part(int action_part, int ilevel, int head_idx, int num_parts,
				    float skip1, float skip2, float skip3, float dx_loc,
                                    float *box_size, int *periodic, float *dtnew, float *dtold)
{
    if (num_parts <= 0) return;
    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_kick_drift_part];
    [enc setBuffer:s_xp              offset:0 atIndex:0];
    [enc setBuffer:s_vp              offset:0 atIndex:1];
    [enc setBuffer:s_levelp          offset:0 atIndex:2];
    [enc setBuffer:s_f_grav          offset:0 atIndex:3];
    [enc setBuffer:s_hash_key        offset:0 atIndex:4];
    [enc setBuffer:s_hash_val        offset:0 atIndex:5];
    [enc setBuffer:s_ckey_max_dev    offset:0 atIndex:6];
    [enc setBuffer:s_key_off_dev     offset:0 atIndex:7];
    [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:8];
    [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:9];
    [enc setBytes:box_size           length:3 * sizeof(float) atIndex:10];
    [enc setBytes:periodic           length:3 * sizeof(int)   atIndex:11];
    [enc setBytes:dtnew              length:(ilevel + 2) * sizeof(float) atIndex:12];
    [enc setBytes:dtold              length:(ilevel + 2) * sizeof(float) atIndex:13];
    [enc setBytes:&s_hash_size       length:sizeof(int) atIndex:14];
    [enc setBytes:&s_ngridmax        length:sizeof(int) atIndex:15];
    [enc setBytes:&skip1             length:sizeof(float) atIndex:16];
    [enc setBytes:&skip2             length:sizeof(float) atIndex:17];
    [enc setBytes:&skip3             length:sizeof(float) atIndex:18];
    [enc setBytes:&dx_loc            length:sizeof(float) atIndex:19];
    [enc setBytes:&action_part       length:sizeof(int) atIndex:20];
    [enc setBytes:&ilevel            length:sizeof(int) atIndex:21];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:22];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:23];
    long leading = s_npartmax;
    [enc setBytes:&leading           length:sizeof(long) atIndex:24];

    NSUInteger tg = 128;
    NSUInteger nb = (num_parts + tg - 1) / tg;
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_newdt_part(int head_idx, int num_parts, float *vmax_out, float *ekin_out)
{
    if (num_parts <= 0) {
        *vmax_out = 0.0f;
        *ekin_out = 0.0f;
        return;
    }
    memset(s_part_scalar_buf.contents, 0, 8);

    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_newdt_part];
    [enc setBuffer:s_vp              offset:0 atIndex:0];
    [enc setBuffer:s_mp              offset:0 atIndex:1];
    [enc setBuffer:s_part_scalar_buf offset:0 atIndex:2];
    [enc setBuffer:s_part_scalar_buf offset:4 atIndex:3];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:4];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:5];
    long leading = s_npartmax;
    [enc setBytes:&leading           length:sizeof(long) atIndex:6];

    NSUInteger tg = 256;
    NSUInteger nb = (num_parts + tg - 1) / tg;
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];

    float *res = (float *)s_part_scalar_buf.contents;
    *vmax_out = res[0];
    *ekin_out = res[1];
}

extern "C" void mtl_split_part(int head_idx, int num_parts, int ilevel,
			       float skip1, float skip2, float skip3, float dx_loc,
			       int *n_fine_out)
{
    if (num_parts <= 0) {
        *n_fine_out = 0;
        return;
    }

    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_bucket_part];
    [enc setBuffer:s_xp              offset:0 atIndex:0];
    [enc setBuffer:s_isp_swap        offset:0 atIndex:1];
    [enc setBuffer:s_grid            offset:0 atIndex:2];
    [enc setBuffer:s_hash_key        offset:0 atIndex:3];
    [enc setBuffer:s_hash_val        offset:0 atIndex:4];
    [enc setBuffer:s_ckey_max_dev    offset:0 atIndex:5];
    [enc setBuffer:s_key_off_dev     offset:0 atIndex:6];
    [enc setBytes:&s_hash_size       length:sizeof(int) atIndex:7];
    [enc setBytes:&skip1             length:sizeof(float) atIndex:8];
    [enc setBytes:&skip2             length:sizeof(float) atIndex:9];
    [enc setBytes:&skip3             length:sizeof(float) atIndex:10];
    [enc setBytes:&dx_loc            length:sizeof(float) atIndex:11];
    [enc setBytes:&ilevel            length:sizeof(int) atIndex:12];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:13];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:14];
    long leading = s_npartmax;
    [enc setBytes:&leading           length:sizeof(long) atIndex:15];

    NSUInteger tg = 128;
    NSUInteger nb = (num_parts + tg - 1) / tg;
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];

    // init prefix sum
    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_init_ps_hilbert_part];
    [enc setBuffer:s_isp_swap        offset:0 atIndex:0];
    [enc setBuffer:s_sortp           offset:0 atIndex:1];
    [enc setBuffer:s_prefix_sum_part offset:0 atIndex:2];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:3];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];

    // Compute prefix sum using the generic helper
    mtl_prefix_scan_buffer(s_prefix_sum_part, s_partial_sums_part, s_partial_sums_part2,
			   s_partial_sums_part3, s_partial_sums_part4, head_idx - 1, num_parts);

    int n_fine = ((int *)s_prefix_sum_part.contents)[(head_idx - 1) + num_parts - 1];
    *n_fine_out = n_fine;
    int n_coarse = num_parts - n_fine;

    // write swap global partition
    cmd = [s_queue commandBuffer];
    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_write_swap_partition];
    [enc setBuffer:s_isp_swap        offset:0 atIndex:0];
    [enc setBuffer:s_sortp           offset:0 atIndex:1];
    [enc setBuffer:s_prefix_sum_part offset:0 atIndex:2];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:3];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:4];
    [enc setBytes:&n_coarse          length:sizeof(int) atIndex:5];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];

    // write back to sortp
    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_write_sortp_part];
    [enc setBuffer:s_sortp           offset:0 atIndex:0];
    [enc setBuffer:s_isp_swap        offset:0 atIndex:1];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:2];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];

    // Shuffle data (gather/scatter for real columns and real 1D, int 1D)
    for (int idim = 1; idim <= 3; idim++) {
        enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_gather_real_col];
        [enc setBuffer:s_xp_swap         offset:0 atIndex:0];
        [enc setBuffer:s_xp              offset:0 atIndex:1];
        [enc setBuffer:s_sortp           offset:0 atIndex:2];
        [enc setBytes:&head_idx          length:sizeof(int) atIndex:3];
        [enc setBytes:&num_parts         length:sizeof(int) atIndex:4];
        [enc setBytes:&leading           length:sizeof(long) atIndex:5];
        [enc setBytes:&idim              length:sizeof(int) atIndex:6];
        [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding];

        enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_scatter_real_col];
        [enc setBuffer:s_xp              offset:0 atIndex:0];
        [enc setBuffer:s_xp_swap         offset:0 atIndex:1];
        [enc setBytes:&head_idx          length:sizeof(int) atIndex:2];
        [enc setBytes:&num_parts         length:sizeof(int) atIndex:3];
        [enc setBytes:&leading           length:sizeof(long) atIndex:4];
        [enc setBytes:&idim              length:sizeof(int) atIndex:5];
        [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding];
    }

    for (int idim = 1; idim <= 3; idim++) {
        enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_gather_real_col];
        [enc setBuffer:s_xp_swap         offset:0 atIndex:0];
        [enc setBuffer:s_vp              offset:0 atIndex:1];
        [enc setBuffer:s_sortp           offset:0 atIndex:2];
        [enc setBytes:&head_idx          length:sizeof(int) atIndex:3];
        [enc setBytes:&num_parts         length:sizeof(int) atIndex:4];
        [enc setBytes:&leading           length:sizeof(long) atIndex:5];
        [enc setBytes:&idim              length:sizeof(int) atIndex:6];
        [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding];

        enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_scatter_real_col];
        [enc setBuffer:s_vp              offset:0 atIndex:0];
        [enc setBuffer:s_xp_swap         offset:0 atIndex:1];
        [enc setBytes:&head_idx          length:sizeof(int) atIndex:2];
        [enc setBytes:&num_parts         length:sizeof(int) atIndex:3];
        [enc setBytes:&leading           length:sizeof(long) atIndex:4];
        [enc setBytes:&idim              length:sizeof(int) atIndex:5];
        [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding];
    }

    // Masses
    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gather_real_1d];
    [enc setBuffer:s_xp_swap         offset:0 atIndex:0];
    [enc setBuffer:s_mp              offset:0 atIndex:1];
    [enc setBuffer:s_sortp           offset:0 atIndex:2];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:3];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];

    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_scatter_real_1d];
    [enc setBuffer:s_mp              offset:0 atIndex:0];
    [enc setBuffer:s_xp_swap         offset:0 atIndex:1];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:2];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];

    // IDs (idp_swap as tmp, but s_idp is 32-bit int)
    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gather_i4_1d];
    [enc setBuffer:s_idp_swap        offset:0 atIndex:0];
    [enc setBuffer:s_idp             offset:0 atIndex:1];
    [enc setBuffer:s_sortp           offset:0 atIndex:2];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:3];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];

    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_scatter_i4_1d];
    [enc setBuffer:s_idp             offset:0 atIndex:0];
    [enc setBuffer:s_idp_swap        offset:0 atIndex:1];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:2];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];

    // Levels (levelp)
    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_gather_i4_1d];
    [enc setBuffer:s_idp_swap        offset:0 atIndex:0];
    [enc setBuffer:s_levelp          offset:0 atIndex:1];
    [enc setBuffer:s_sortp           offset:0 atIndex:2];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:3];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:4];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];

    enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_scatter_i4_1d];
    [enc setBuffer:s_levelp          offset:0 atIndex:0];
    [enc setBuffer:s_idp_swap        offset:0 atIndex:1];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:2];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:3];
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding];

    [cmd commit]; [cmd waitUntilCompleted];
}

extern "C" void mtl_sort_part(int head_idx, int num_parts, int level, float shift, float dx_inv, float *skip)
{
    if (num_parts <= 0) return;

    NSUInteger tg = 128;
    NSUInteger nb = (num_parts + tg - 1) / tg;

    id<MTLCommandBuffer> cmd = [s_queue commandBuffer];

    // Initialize s_sortp to identity permutation
    {
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_init_swap_table];
        [enc setBuffer:s_sortp           offset:0 atIndex:0];
        [enc setBytes:&head_idx          length:sizeof(int) atIndex:1];
        [enc setBytes:&num_parts         length:sizeof(int) atIndex:2];
        [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding];
    }

    {
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:s_pso_hkey_part];
        [enc setBuffer:s_xp              offset:0 atIndex:0];
        [enc setBuffer:s_idp_swap        // idp_swap used to hold 64-bit long keys
                                         offset:0 atIndex:1];
        [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:2];
        [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:3];
        [enc setBytes:skip               length:3 * sizeof(float) atIndex:4];
        [enc setBytes:&dx_inv            length:sizeof(float) atIndex:5];
        [enc setBytes:&head_idx          length:sizeof(int) atIndex:6];
        [enc setBytes:&num_parts         length:sizeof(int) atIndex:7];
        long leading = s_npartmax;
        [enc setBytes:&leading           length:sizeof(long) atIndex:8];
        [enc setBytes:&level             length:sizeof(int) atIndex:9];
        [enc setBytes:&shift             length:sizeof(float) atIndex:10];
        [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
        [enc endEncoding];
    }

    // Radix sort via bit loop
    for (int ibit = 0; ibit < 3 * level; ibit++) {
        {
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:s_pso_prefix_part_bit];
            [enc setBuffer:s_idp_swap        offset:0 atIndex:0];
            [enc setBuffer:s_sortp           offset:0 atIndex:1];
            [enc setBuffer:s_prefix_sum_part offset:0 atIndex:2];
            [enc setBytes:&head_idx          length:sizeof(int) atIndex:3];
            [enc setBytes:&num_parts         length:sizeof(int) atIndex:4];
            [enc setBytes:&ibit              length:sizeof(int) atIndex:5];
            [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
            [enc endEncoding];
        }

        mtl_prefix_scan_buffer_cb(s_prefix_sum_part, s_partial_sums_part, s_partial_sums_part2,
				  s_partial_sums_part3, s_partial_sums_part4, head_idx - 1, num_parts, cmd);

        {
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:s_pso_local_swap];
            [enc setBuffer:s_isp_swap        offset:0 atIndex:0];
            [enc setBuffer:s_sortp           offset:0 atIndex:1];
            [enc setBuffer:s_prefix_sum_part offset:0 atIndex:2];
            [enc setBytes:&head_idx          length:sizeof(int) atIndex:3];
            [enc setBytes:&num_parts         length:sizeof(int) atIndex:4];
            [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
            [enc endEncoding];
        }

        {
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:s_pso_global_swap];
            [enc setBuffer:s_sortp           offset:0 atIndex:0];
            [enc setBuffer:s_isp_swap        offset:0 atIndex:1];
            [enc setBytes:&head_idx          length:sizeof(int) atIndex:2];
            [enc setBytes:&num_parts         length:sizeof(int) atIndex:3];
            [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
            [enc endEncoding];
        }
    }

    [cmd commit];
    [cmd waitUntilCompleted];
}

extern "C" void mtl_cic_part_medium(
    int head_idx, int num_parts,
    float skip1, float skip2, float skip3,
    float dx_loc, float vol_loc, float mass_sph,
    int star, float m_refine_at_level, float mass_cut_refine, int ilevel)
{
    if (num_parts <= 0) return;

    id<MTLCommandBuffer>         cmd = [s_queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:s_pso_cic_part_medium];
    [enc setBuffer:s_sortp           offset:0 atIndex:0];
    [enc setBuffer:s_hash_key        offset:0 atIndex:1];
    [enc setBuffer:s_hash_val        offset:0 atIndex:2];
    [enc setBuffer:s_ckey_max_dev    offset:0 atIndex:3];
    [enc setBuffer:s_key_off_dev     offset:0 atIndex:4];
    [enc setBuffer:s_box_ckey_min_dev offset:0 atIndex:5];
    [enc setBuffer:s_box_ckey_max_dev offset:0 atIndex:6];
    [enc setBuffer:s_xp              offset:0 atIndex:7];
    [enc setBuffer:s_mp              offset:0 atIndex:8];
    [enc setBuffer:s_rho             offset:0 atIndex:9];
    [enc setBuffer:s_nref            offset:0 atIndex:10];
    [enc setBytes:&s_hash_size       length:sizeof(int) atIndex:11];
    [enc setBytes:&skip1             length:sizeof(float) atIndex:12];
    [enc setBytes:&skip2             length:sizeof(float) atIndex:13];
    [enc setBytes:&skip3             length:sizeof(float) atIndex:14];
    [enc setBytes:&dx_loc            length:sizeof(float) atIndex:15];
    [enc setBytes:&vol_loc           length:sizeof(float) atIndex:16];
    [enc setBytes:&mass_sph          length:sizeof(float) atIndex:17];
    [enc setBytes:&star              length:sizeof(int) atIndex:18];
    [enc setBytes:&m_refine_at_level length:sizeof(float) atIndex:19];
    [enc setBytes:&mass_cut_refine   length:sizeof(float) atIndex:20];
    [enc setBytes:&ilevel            length:sizeof(int) atIndex:21];
    long leading = s_npartmax;
    [enc setBytes:&leading           length:sizeof(long) atIndex:22];
    [enc setBytes:&head_idx          length:sizeof(int) atIndex:23];
    [enc setBytes:&num_parts         length:sizeof(int) atIndex:24];

    NSUInteger tg = 256;
    NSUInteger nb = (num_parts + tg - 1) / tg;
    [enc dispatchThreadgroups:{nb,1,1} threadsPerThreadgroup:{tg,1,1}];
    [enc endEncoding]; [cmd commit]; [cmd waitUntilCompleted];
}

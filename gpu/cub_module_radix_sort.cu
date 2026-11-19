#include "cub_module_radix_sort.h"
#include <iostream>
#include <cub/cub.cuh>
#include <cstdint>

/// working space for device scan
void *g_cub_module_radix_sort_workspace = nullptr;
size_t g_cub_module_radix_sort_work_size = 0;
int g_cub_module_radix_sort_verbose = 0;


#define CUDA_CHECK(call) \
{ \
  cudaError_t ierr = call; \
  if (cudaSuccess != ierr) \
  { \
    std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << std::endl \
              << "ERROR: " << cudaGetErrorString(ierr) << std::endl; \
    abort(); \
  } \
}

///  checks the alighnment of the passed ptr and aborts if it is not correctly aligned
#if !defined(NDEBUG)
// @returns true if the pointer is aligned to it's type's size
template<typename T>
bool is_aligned(T *ptr)
{
  if (0 == uintptr_t(ptr) % sizeof(T))
    return true;
  return false;
}

#define CHECK_ALIGN(var) \
if (!is_aligned(var)) \
{ \
  std::cerr << "ERROR: " #var " is not " << sizeof(decltype(var)) \
            << " byte aligned." << std::endl; \
  abort(); \
}

#else
#define CHECK_ALIGN(var)
#endif

// --------------------------------------------------------------------------
void cub_module_radix_sort_allocate_workspace(int32_t num_elem)
{
  // check if we have enough work space
  size_t work_size = 0;
#if defined(OPT_MEM)
  // DoubleBuffer overload: CUB ping-pongs between the caller's buffers and
  // does not allocate its own copy of the keys and values
  cub::DoubleBuffer<int64_t> d_keys(nullptr, nullptr);
  cub::DoubleBuffer<int32_t> d_vals(nullptr, nullptr);
  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(nullptr, work_size, d_keys, d_vals, num_elem))
#else
  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(nullptr, work_size,
                                             (int64_t*)nullptr, (int64_t*)nullptr,
                                             (int32_t*)nullptr, (int32_t*)nullptr,
                                             num_elem))
#endif
  if (g_cub_module_radix_sort_work_size < work_size)
  {
    // free old work space
    if (g_cub_module_radix_sort_workspace)
    {
      cudaFree(g_cub_module_radix_sort_workspace);
      g_cub_module_radix_sort_workspace = nullptr;
      g_cub_module_radix_sort_work_size = 0;
    }

    // allocate new work space
    g_cub_module_radix_sort_work_size = work_size;
    CUDA_CHECK(cudaMalloc(&g_cub_module_radix_sort_workspace,
                          g_cub_module_radix_sort_work_size))
    if (0 < g_cub_module_radix_sort_verbose)
      std::cerr << "allocated " << work_size << " bytes on the device" << std::endl;
  }
}

/// input an array of keys to sort
void cub_module_radix_sort(int64_t *keys_in, int64_t *keys_out,
                           int32_t *idx_in, int32_t *idx_out,
                           int32_t idx0, int32_t num_elem,
                           int32_t start_bit, int32_t end_bit)
{
  if (1 < g_cub_module_radix_sort_verbose)
    std::cerr << "cub_module_radix_sort" << std::endl;

  // check for misaligned buffers.
  CHECK_ALIGN(keys_in)
  CHECK_ALIGN(keys_out)
  CHECK_ALIGN(idx_in)
  CHECK_ALIGN(idx_out)

  // allocate workspace
  cub_module_radix_sort_allocate_workspace(num_elem);

#if defined(OPT_MEM)
  // in-place sort: keys_in and idx_in are overwritten
  cub::DoubleBuffer<int64_t> d_keys(keys_in, keys_out);
  cub::DoubleBuffer<int32_t> d_vals(idx_in, idx_out);
  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(g_cub_module_radix_sort_workspace,
                                             g_cub_module_radix_sort_work_size,
                                             d_keys, d_vals, num_elem, start_bit, end_bit))
  // the sorted keys and indices may end up in the input buffers
  if (d_keys.Current() != keys_out)
    CUDA_CHECK(cudaMemcpyAsync(keys_out, d_keys.Current(), num_elem * sizeof(int64_t),
                               cudaMemcpyDeviceToDevice))
  if (d_vals.Current() != idx_out)
    CUDA_CHECK(cudaMemcpyAsync(idx_out, d_vals.Current(), num_elem * sizeof(int32_t),
                               cudaMemcpyDeviceToDevice))
#else
  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(g_cub_module_radix_sort_workspace,
                                             g_cub_module_radix_sort_work_size,
                                             keys_in, keys_out, idx_in, idx_out,
                                             num_elem, start_bit, end_bit))
#endif
}

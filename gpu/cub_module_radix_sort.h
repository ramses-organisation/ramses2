#ifndef cub_module_radix_sort_h
#define cub_module_radix_sort_h

#include <cstdint>

extern "C"
{
/** sorts keys and generates a permutation array.
 * @param[in] keys_in keys to sort (overwritten when compiled with OPT_MEM)
 * @param[out] keys_out array for the sorted keys
 * @param[in] idx_in an empty array of the same length as keys_in (overwritten when compiled with OPT_MEM)
 * @param[out] idx_out array for the sorted indices. this is the result
 * @param[in] head_idx the 1 based index into the keys
 * @param[in] num_elem the number of keys to sort
 * @param[in] start_bit the first bit to use in the radix sort
 * @param[in] end_bit the last bit to use in the radix sort
 */
void cub_module_radix_sort(int64_t *keys_in, int64_t *keys_out,
                           int32_t *idx_in, int32_t *idx_out,
                           int32_t head_idx,  int32_t num_elem,
                           int32_t start_bit, int32_t end_bit);

/** allocate the workspace for a sort of num_elem.
 * note: the sort call will always check if there is enough space and allocate
 * as needed.
 */
void cub_module_radix_sort_allocate_workspace(int32_t num_elem);
}

#endif

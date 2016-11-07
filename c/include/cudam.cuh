#ifndef __CUDAM_CUH__
#define __CUDAM_CUH__
#define HIGH_BIT 0X8000000000000000

typedef uint64_t* cudam_t;

__host__ __device__ cudam_t cudam_init();
__host__ __device__ int cudam_cmp(cudam_t a, cudam_t b);
__host__ __device__ void cudam_add(cudam_t c, cudam_t a, cudam_t b);
__host__ __device__ void cudam_add_ui(cudam_t c, cudam_t a, uint64_t b);
__host__ __device__ void cudam_sll(cudam_t a, int shift);
__host__ __device__ void cudam_set(cudam_t a, cudam_t b);
__host__ __device__ void cudam_zero(cudam_t a);
__host__ __device__ void cudam_mult_ui(cudam_t c, cudam_t a, uint64_t b);
__host__ cudam_t cudam_from_str(char * str);
__host__ __device__ void cudam_sub(cudam_t c, cudam_t a, cudam_t b);
__host__ __device__ void cudam_mod(cudam_t r, cudam_t a, cudam_t m);
__host__ __device__ void cudam_multm(cudam_t c,cudam_t a, cudam_t b, cudam_t p);
__host__ __device__ void cudam_mult(cudam_t c, cudam_t a, cudam_t b);
__host__ __device__ int cudam_is_zero(cudam_t a);

#endif

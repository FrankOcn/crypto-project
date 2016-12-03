#ifndef __VEC_MPZ_H__
#define __VEC_MPZ_H__
#include "stdint.h"
#include <gmp.h>

struct vec_mpz_t {
  mpz_t * elems;
  int size;
  int capacity;
};

struct vec_mpz_t* vec_mpz_init();
struct vec_mpz_t* vec_mpz_init_size(size_t);
void vec_mpz_insert(struct vec_mpz_t*, mpz_t);
void vec_mpz_free(struct vec_mpz_t* vec);
struct vec_mpz_t* vec_mpz_product_tree(struct vec_mpz_t* vec, const unsigned int BATCH_SIZE);
struct vec_mpz_t* vec_mpz_remainder(mpz_t n, struct vec_mpz_t* vec, const unsigned int BATCH_SIZE);
struct vec_mpz_t* remainders_in(struct vec_mpz_t* remainder_tree, const unsigned int BATCH_SIZE);
void compute_smooth(struct vec_mpz_t* remainders, struct vec_mpz_t* batch, int* smooth_bool_arr, const unsigned int BATCH_SIZE);
int* batch_smooth_parts(mpz_t z, struct vec_mpz_t* primes_product_tree, struct vec_mpz_t* batch, const unsigned int BATCH_SIZE);

#endif

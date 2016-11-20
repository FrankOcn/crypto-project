#ifndef __VEC_MPZ_H__
#define __VEC_MPZ_H__
#include "stdint.h"
#include "gmp.h"

struct vec_mpz_t {
  mpz_t * elems;
  int size;
  int capacity;
};

struct vec_mpz_t* vec_mpz_init();
struct vec_mpz_t* vec_mpz_init_size(size_t);
void vec_mpz_insert(struct vec_mpz_t*, mpz_t);
void vec_mpz_free(struct vec_mpz_t* vec);
struct vec_mpz_t* vec_mpz_product_tree(struct vec_mpz_t* vec);

#endif

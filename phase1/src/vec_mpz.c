#include "vec_mpz.h"
#include "stdlib.h"
#include "gmp.h"

struct vec_mpz_t* vec_mpz_init() {
  struct vec_mpz_t* vec = (struct vec_mpz_t*)malloc(sizeof(struct vec_mpz_t));
  vec->capacity = 8;
  vec->size = 0;
  vec->elems = (mpz_t *)malloc(vec->capacity * sizeof(mpz_t));
  int i;
  for (i = 0; i < vec->capacity; ++i ) {
    mpz_init(vec->elems[i]);
  }
  return vec;
}

void vec_mpz_insert(struct vec_mpz_t* vec, mpz_t elem) {
  int i;
  if (vec->size == vec->capacity) {
    vec->capacity = vec->capacity * 2;
    mpz_t* elems = (mpz_t *)malloc(vec->capacity * sizeof(mpz_t));
    for (i = 0; i < vec->size; ++i) {
      mpz_init_set(elems[i], vec->elems[i]);
      mpz_clear(vec->elems[i]);
    }
    free(vec->elems);
    vec->elems = elems;
  }
  mpz_set(vec->elems[vec->size], elem);
  vec->size++;
}

void vec_mpz_free(struct vec_mpz_t* vec) {
  int i;
  for (i = 0; i < vec->size; ++i) {
    mpz_clear(vec->elems[i]);
  }
  free(vec->elems);
  free(vec);
}

void vec_mpz_product_tree(struct vec_mpz_t* vec) {
  int i = 0;
  mpz_t tmp;
  mpz_init(tmp);
  for (i = 0; i < vec->size - 1; i += 2) {
    mpz_mul(tmp, vec->elems[i], vec->elems[i+1]);
    vec_mpz_insert(vec, tmp);
  }
  mpz_clear(tmp);
}
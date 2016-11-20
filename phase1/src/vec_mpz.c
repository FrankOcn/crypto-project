#include "vec_mpz.h"
#include "stdlib.h"
#include "gmp.h"

struct vec_mpz_t* vec_mpz_init() {
  struct vec_mpz_t* vec = (struct vec_mpz_t*)malloc(sizeof(struct vec_mpz_t));
  vec->capacity = 32;
  vec->size = 0;
  vec->elems = (mpz_t *)malloc(vec->capacity * sizeof(mpz_t));
  int i;
  for (i = 0; i < 32; ++i ) {
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
      mpz_init(elems[i]);
      mpz_set(elems[i], vec->elems[i]);
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

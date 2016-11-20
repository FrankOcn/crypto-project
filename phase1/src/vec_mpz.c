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

struct vec_mpz_t* vec_mpz_init_size(size_t capacity) {
  struct vec_mpz_t* vec = (struct vec_mpz_t*)malloc(sizeof(struct vec_mpz_t));
  vec->capacity = capacity;
  vec->size = capacity;
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
    for (i = vec->size; i < vec->capacity; ++i) {
      mpz_init(elems[i]);
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

struct vec_mpz_t* vec_mpz_product_tree(struct vec_mpz_t* vec) {
  int i;
  struct vec_mpz_t* product_tree = vec_mpz_init_size((vec->size * 2) - 1);
  for (i = 0; i < vec->size; ++i) {
    mpz_set(product_tree->elems[i + vec->size - 1], vec->elems[i]);
  }
  for (i = product_tree->size - vec->size - 1; i >= 0; --i) {
    mpz_mul(
      product_tree->elems[i],
      product_tree->elems[i*2 + 1],
      product_tree->elems[i*2 + 2]
    );
  }
  return product_tree;
}

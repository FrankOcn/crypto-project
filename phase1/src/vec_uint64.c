#include "vec_uint64.h"
#include "stdlib.h"
#include "stdint.h"

struct vec_uint64_t* vec_uint64_init() {
  struct vec_uint64_t* vec = (struct vec_uint64_t*)malloc(sizeof(struct vec_uint64_t));
  vec->capacity = 8;
  vec->size = 0;
  vec->elems = (uint64_t *)malloc(vec->capacity * sizeof(uint64_t));
  return vec;
}

void vec_uint64_insert(struct vec_uint64_t* vec, uint64_t elem) {
  int i;
  if (vec->size == vec->capacity) {
    vec->capacity = vec->capacity * 2;
    uint64_t* elems = (uint64_t *)malloc(vec->capacity * sizeof(uint64_t));
    for (i = 0; i < vec->size; ++i) {
      elems[i] = vec->elems[i];
    }
    free(vec->elems);
    vec->elems = elems;
  }
  vec->elems[vec->size] = elem;
  vec->size++;
}

void vec_uint64_free(struct vec_uint64_t* vec) {
  free(vec->elems);
  free(vec);
}

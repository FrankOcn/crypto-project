#ifndef __VEC_UINT64_H__
#define __VEC_UINT64_H__
#include "stdint.h"

struct vec_uint64_t {
  uint64_t * elems;
  int size;
  int capacity;
};

struct vec_uint64_t* vec_uint64_init();
void vec_uint64_insert(struct vec_uint64_t*, uint64_t);
void vec_uint64_free(struct vec_uint64_t* vec);

#endif

#ifndef __HASHTABLE_UINT64_H__
#define __HASHTABLE_UINT64_H__
#include "stdint.h"

struct hashtable_entry_uint64_t {
  uint64_t key;
  uint64_t value;
};

struct hashtable_uint64_t {
  struct hashtable_entry_uint64_t ** elems;
  int size;
  int capacity;
};

struct hashtable_uint64_t* hashtable_uint64_init();
void hashtable_uint64_insert(struct hashtable_uint64_t*, uint64_t, uint64_t);
uint64_t* hashtable_uint64_retrieve(struct hashtable_uint64_t*, uint64_t);
void hashtable_uint64_free(struct hashtable_uint64_t*);
#endif

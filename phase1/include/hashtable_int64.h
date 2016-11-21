#ifndef __HASHTABLE_INT64_H__
#define __HASHTABLE_INT64_H__
#include "stdint.h"

struct hashtable_entry_int64_t {
  int64_t key;
  int64_t value;
};

struct hashtable_int64_t {
  struct hashtable_entry_int64_t ** elems;
  int size;
  int capacity;
};

struct hashtable_int64_t* hashtable_int64_init();
void hashtable_int64_insert(struct hashtable_int64_t*, int64_t, int64_t);
int64_t* hashtable_int64_retrieve(struct hashtable_int64_t*, int64_t);
void hashtable_int64_free(struct hashtable_int64_t*);
#endif

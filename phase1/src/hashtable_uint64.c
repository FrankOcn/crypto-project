#include "hashtable_uint64.h"
#include "stdlib.h"
#include "stdint.h"

struct hashtable_uint64_t* hashtable_uint64_init() {
  int i;
  struct hashtable_uint64_t* ht = (struct hashtable_uint64_t*)malloc(sizeof(struct hashtable_uint64_t));
  ht->capacity = 8;
  ht->size = 0;
  ht->elems = (struct hashtable_entry_uint64_t **)malloc(ht->capacity * sizeof(struct hashtable_entry_uint64_t*));
  for (i = 0; i < ht->size; ++i) {
    ht->elems[i] = 0;
  }
  return ht;
}

void hashtable_uint64_free(struct hashtable_uint64_t* ht) {
  int i;
  for (i = 0; i < ht->capacity; ++i) {
    if (ht->elems[i] != 0) free(ht->elems[i]);
  }
  free(ht->elems);
  free(ht);
}

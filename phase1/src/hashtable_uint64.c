#include "hashtable_uint64.h"
#include "stdlib.h"
#include "stdint.h"

struct hashtable_uint64_t* hashtable_uint64_init() {
  int i;
  struct hashtable_uint64_t* ht = (struct hashtable_uint64_t*)malloc(sizeof(struct hashtable_uint64_t));
  ht->capacity = 10;
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

void hashtable_uint64_insert(struct hashtable_uint64_t* table, uint64_t key, uint64_t value) {
  int i;
  size_t location, probe, new_capacity;
  if (table->size >= table->capacity/2) {
    new_capacity = table->capacity * 2;
    struct hashtable_entry_uint64_t** elems = (struct hashtable_entry_uint64_t **)malloc(new_capacity * sizeof(struct hashtable_entry_uint64_t*));
    for (i = 0; i < table->capacity; ++i) {
      if (table->elems[i] != 0) {
        location = table->elems[i]->key % new_capacity;
        probe = 0;
        while (elems[(location + probe)%new_capacity] != 0) {
          probe = probe == 0 ? 1 : 2*probe;
        }
        elems[(location + probe)%new_capacity] = table->elems[i];
        table->elems[i] = 0;
      }
    }
    free(table->elems);
    table->elems = elems;
    table->capacity = new_capacity;
  }
  location = key % table->capacity;
  probe = 0;
  while (table->elems[(location + probe)%table->capacity] != 0) {
    probe = probe == 0 ? 1 : 2*probe;
  }
  table->elems[(location + probe)%table->capacity] = (struct hashtable_entry_uint64_t *)malloc(sizeof(struct hashtable_entry_uint64_t*));
  table->elems[(location + probe)%table->capacity]->key = key;
  table->elems[(location + probe)%table->capacity]->value = value;
  table->size += 1;
}

uint64_t* hashtable_uint64_retrieve(struct hashtable_uint64_t* table, uint64_t key) {
  size_t location = key % table->capacity;
  size_t probe = 0;
  while (table->elems[location + probe] != 0 && table->elems[location + probe]->key != key) {
    probe = probe == 0 ? 1 : 2*probe;
  }
  if (table->elems[location + probe] == 0) return 0;
  return &(table->elems[location + probe]->value);
}

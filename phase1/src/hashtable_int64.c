#include "hashtable_int64.h"
#include "stdlib.h"
#include "stdint.h"

struct hashtable_int64_t* hashtable_int64_init() {
  int i;
  struct hashtable_int64_t* ht = (struct hashtable_int64_t*)malloc(sizeof(struct hashtable_int64_t));
  ht->capacity = 10;
  ht->size = 0;
  ht->elems = (struct hashtable_entry_int64_t **)malloc(ht->capacity * sizeof(struct hashtable_entry_int64_t*));
  for (i = 0; i < ht->capacity; ++i) {
    ht->elems[i] = 0;
  }
  return ht;
}

void hashtable_int64_free(struct hashtable_int64_t* ht) {
  int i;
  for (i = 0; i < ht->capacity; ++i) {
    if (ht->elems[i] != 0) free(ht->elems[i]);
  }
  free(ht->elems);
  free(ht);
}

void hashtable_int64_insert(struct hashtable_int64_t* table, int64_t key, int64_t value) {
  int i;
  size_t location, probe, new_capacity;
  if (table->size >= table->capacity/2) {
    new_capacity = table->capacity * 2;
    struct hashtable_entry_int64_t** elems = (struct hashtable_entry_int64_t **)malloc(new_capacity * sizeof(struct hashtable_entry_int64_t*));
    for (i = 0; i < new_capacity; ++i) {
      elems[i] = 0;
    }
    for (i = 0; i < table->capacity; ++i) {
      if (table->elems[i] != 0) {
        location = table->elems[i]->key % new_capacity;
        probe = 0;
        while (elems[(location + (probe*probe))%new_capacity] != 0) {
          probe += 1;
        }
        elems[(location + (probe*probe))%new_capacity] = table->elems[i];
        table->elems[i] = 0;
      }
    }
    free(table->elems);
    table->elems = elems;
    table->capacity = new_capacity;
  }
  location = key % table->capacity;
  probe = 0;
  while (table->elems[(location + (probe*probe))%table->capacity] != 0) {
    probe += 1;
  }
  table->elems[(location + (probe*probe))%table->capacity] = (struct hashtable_entry_int64_t *)malloc(sizeof(struct hashtable_entry_int64_t*));
  table->elems[(location + (probe*probe))%table->capacity]->key = key;
  table->elems[(location + (probe*probe))%table->capacity]->value = value;
  table->size += 1;
}

int64_t* hashtable_int64_retrieve(struct hashtable_int64_t* table, int64_t key) {
  size_t location = key % table->capacity;
  size_t probe = 0;
  while (table->elems[(location + (probe*probe))%table->capacity] != 0 && table->elems[(location + (probe*probe))% table->capacity]->key != key) {
    probe += 1;
  }
  if (table->elems[(location + (probe*probe))%table->capacity] == 0) return 0;
  return &(table->elems[(location + (probe*probe))%table->capacity]->value);
}

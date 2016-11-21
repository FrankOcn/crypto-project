#include "vec_mpz.h"
#include "stdlib.h"
#include "math.h"
#include "gmp.h"

unsigned int log_2(int v)
{
  if (v == 1) return 0;
  unsigned int ret = 0;
  while (v)
  {
    v >>= 1;
    ret++;
  }
  return ret;
}

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

// Given integer n, and sequence of integers X[]:
// returns vec of [ (n mod X[0]), (n mod X[1]), (n mod X[2])..  ]
struct vec_mpz_t* vec_mpz_remainder_tree(mpz_t n, struct vec_mpz_t* vec) 
{
  int i;
  struct vec_mpz_t* remainder_tree = vec_mpz_init_size(vec->size);
  
  // gmp_printf("%Zd mod %Zd\n", n, vec->elems[0]);
  mpz_mod(remainder_tree->elems[0], n, vec->elems[0]);
  for (i = 1; i < remainder_tree->size; i++)
  {
    if (i % 2 == 0)
    {
      // gmp_printf("%Zd mod %Zd\n", remainder_tree->elems[i/2 -1], vec->elems[i]);
      mpz_mod(remainder_tree->elems[i], remainder_tree->elems[i/2 - 1], vec->elems[i]);
    }
    else
    {
      // gmp_printf("%Zd mod %Zd\n", remainder_tree->elems[i/2], vec->elems[i]);
      mpz_mod(remainder_tree->elems[i], remainder_tree->elems[i/2], vec->elems[i]);
    }
  }
  return remainder_tree;
}

struct vec_mpz_t* remainders_in(struct vec_mpz_t* remainder_tree, int batch_size)
{
  struct vec_mpz_t* remainders = vec_mpz_init_size(batch_size);
  int j = remainder_tree->size - 1;
  for (int i = 0; i < batch_size; i++)
  {
    mpz_set(remainders->elems[i], remainder_tree->elems[j - i]);
  }
  return remainders;
}

// iterate e until 2^(2^e) > x_k
// then push remainder[i]^(2^2)
struct vec_mpz_t* compute_squares(struct vec_mpz_t* remainders, struct vec_mpz_t* batch)
{
  struct vec_mpz_t* squares = vec_mpz_init();
  mpz_t temp, two_sq_sq;
  mpz_init(temp);
  mpz_init(two_sq_sq);
  uint64_t two_sq;

  for (int i = 0; i < remainders->size; i++)
  {
    two_sq = 1;
    mpz_ui_pow_ui(two_sq_sq, 2, two_sq);
    while (mpz_cmp(two_sq_sq, batch->elems[i]) < 0)
    {
      two_sq <<= 1;
      mpz_ui_pow_ui(two_sq_sq, 2, two_sq);
    }
    mpz_powm_ui(temp, remainders->elems[i], two_sq, batch->elems[i]);
    vec_mpz_insert(squares, temp);
  }

  return squares;
}

struct vec_mpz_t* compute_smooth_parts(struct vec_mpz_t* batch, struct vec_mpz_t* squares)
{
  mpz_t temp;
  mpz_init(temp);
  struct vec_mpz_t* smooth_parts = vec_mpz_init_size(batch->size);
  for (int i = 0; i < batch->size; i++)
  {
    mpz_gcd(smooth_parts->elems[i], batch->elems[i], squares->elems[i]);
  }
  return smooth_parts;
}

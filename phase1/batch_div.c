#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"
#include "gmp.h"
#include "pthread.h"
#include "time.h"
#include "vec_mpz.h"
#include "vec_uint64.h"
#include "crypto_utils.h"

int main(int argc, char ** argv)
{
  int i;
  const unsigned int BATCH_SIZE = 8;
  mpz_t a, z;
  mpz_init(a);
  mpz_init(z); // product of all primes in factor base

  struct vec_mpz_t* batch = vec_mpz_init(); // residues 
  struct vec_mpz_t* squares;
  struct vec_mpz_t* primes = vec_mpz_init(); // primes in factor base
  struct vec_mpz_t* primes_product_tree = vec_mpz_init(); // prime prod tree
  struct vec_mpz_t* primes_remainder_tree;
  struct vec_mpz_t* remainders;
  struct vec_mpz_t* smooth_parts = vec_mpz_init_size(BATCH_SIZE);
  int* smooth_bools;

  FILE *fp;
  fp = fopen("./primes.txt", "r");
  while (mpz_inp_str(a, fp, 10))
  {
    vec_mpz_insert(primes, a);
  }
  fclose(fp);

  fp = fopen("./batch.txt", "r");
  while (mpz_inp_str(a, fp, 10))
  {
    vec_mpz_insert(batch, a);
  }
  fclose(fp);

  primes_product_tree = vec_mpz_product_tree(primes);
  mpz_set(z, primes_product_tree->elems[0]);

  primes_remainder_tree = vec_mpz_remainder_tree(z, batch);
  remainders = remainders_in(primes_remainder_tree, BATCH_SIZE);
  squares = compute_squares(remainders, batch);

  smooth_parts = compute_smooth_parts(batch, squares);

  /*
  printf("THE SMOOTH PARTS OF:\n");
  for (i = 0; i < BATCH_SIZE; i++)
  {
    gmp_printf("%Zd\n", batch->elems[i]);
  }
  printf("WITH FACTOR BASE:\n");
  for (i = 0; i < primes->size; i++)
  {
    gmp_printf("%Zd | ", primes->elems[i]);
  }
  printf("\n ARE:\n");
  */
  smooth_bools = smooth_check(smooth_parts);
  for (i = 0; i < BATCH_SIZE; i++)
  {
    gmp_printf("%Zd | %d\n", batch->elems[i], smooth_bools[i]);
  }

  printf("\n");

  vec_mpz_free(batch);
  vec_mpz_free(squares);
  vec_mpz_free(primes);
  vec_mpz_free(primes_product_tree);
  vec_mpz_free(primes_remainder_tree);
  vec_mpz_free(remainders);
  vec_mpz_free(smooth_parts);
}

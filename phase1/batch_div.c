#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"
#include "gmp.h"
#include "pthread.h"
#include "time.h"
#include "vec_mpz.h"
#include "vec_uint64.h"
#include "vec_uint64.h"
#include "crypto_utils.h"
#include "hashtable_int64.h"

int main(int argc, char ** argv)
{
  int i;
  const unsigned int BATCH_SIZE = 1000;
  mpz_t a, z, base, modulus; // z = total product of factor base
  mpz_inits(a, z, base, modulus);

  struct vec_mpz_t* r_batch = vec_mpz_init_size(BATCH_SIZE); // r values
  struct vec_mpz_t* s_batch = vec_mpz_init_size(BATCH_SIZE); // s values
  struct vec_mpz_t* primes = vec_mpz_init(); // primes in factor base
  int* smooth_bools;

  FILE *fp;
  FILE *output;
  output = fopen(filename, "w+");
  fp = fopen("./primes.txt", "r");
  while (mpz_inp_str(a, fp, 10))
  {
    vec_mpz_insert(primes, a);
  }
  fclose(fp);


  smooth_bools = batch_smooth_parts(primes, s_batch);

  for (i = 0; i < BATCH_SIZE; i++)
  {
    gmp_printf("%Zd | %d\n", batch->elems[i], smooth_bools[i]);
  }
  
  printf("\n");

  mpz_clear(a);
  mpz_clear(z);
  vec_mpz_free(batch);
  vec_mpz_free(primes);
}

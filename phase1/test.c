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
  clock_t start, diff;
  unsigned long msec;

  mpz_t slow_prod, fast_prod, temp;
  mpz_init(slow_prod);
  mpz_init(fast_prod);
  mpz_init(temp);
  struct vec_mpz_t* primes = vec_mpz_init();
  struct vec_mpz_t* product_tree;

  mpz_set_ui(slow_prod, 1);
  mpz_set_ui(fast_prod, 1);

  FILE *fp;
  fp = fopen("./primes.txt", "r");
  int j = 0;
  start = clock();
  while (mpz_inp_str(temp, fp, 10))
  {
    j++;
    printf("%d\n", j);
    vec_mpz_insert(primes, temp);
  }
  fclose(fp);

  start = clock();
  for (int i = 0; i < primes->size; i++)
  {
    mpz_set(temp, slow_prod);
    mpz_mul(slow_prod, temp, primes->elems[i]);
  }
  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  printf(" --- %lu seconds --- \n", msec);

  start = clock();
  product_tree = vec_mpz_product_tree(primes, primes->size); 
  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  printf(" --- %lu seconds --- \n", msec);

}

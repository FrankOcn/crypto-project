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

struct shared_data
{
  mpz_t modulus;
  mpz_t base;
  mpz_t z; // z = total product of factor base
  struct vec_mpz_t* primes;
  struct vec_mpz_t* primes_product_tree;
};

int main(int argc, char ** argv)
{
  // struct shared_data* data;

  int i;
  const unsigned int BATCH_SIZE = 100000;
  int NO_SMOOTH_FOUND = 1;

  clock_t start, diff;
  unsigned long msec;

  // data stuff
  mpz_t modulus, base, z;
  struct vec_mpz_t* primes = vec_mpz_init();
  struct vec_mpz_t* primes_product_tree;

  mpz_t a, pow, test_num, test_num_r, test_num_s; // z = total product of factor base
  mpz_inits(a, z, base, modulus, pow, test_num, test_num_r, test_num_s, NULL);

  struct vec_mpz_t* temp = vec_mpz_init();
  struct vec_mpz_t* r_batch = vec_mpz_init_size(BATCH_SIZE); // r values
  struct vec_mpz_t* s_batch = vec_mpz_init_size(BATCH_SIZE); // s values
  struct vec_mpz_t* pow_batch = vec_mpz_init_size(BATCH_SIZE);
  primes = vec_mpz_init(); // primes in factor base
  int* r_smooth_bools;
  int* s_smooth_bools;

  gmp_randstate_t rand_state;
  gmp_randinit_mt(rand_state);
  gmp_randseed_ui(rand_state, (int)time(NULL));

  mpz_set_str(base, "2233690997702420376898541331161315870090615501504568211724152", 10);
  mpz_set_str(modulus, "3217014639198601771090467299986349868436393574029172456674199", 10);
 

  FILE *fp;
  FILE *output;
  output = fopen("./output.txt", "w+");
  fp = fopen("./primes_80k.txt", "r");
  while (mpz_inp_str(a, fp, 10))
  {
    vec_mpz_insert(primes, a);
  }
  fclose(fp);

  primes_product_tree = vec_mpz_product_tree(primes, primes->size);
  mpz_set(z, primes_product_tree->elems[0]);

  printf("Starting process.. \n");
  start = clock();
  while (NO_SMOOTH_FOUND)
  {
    r_batch->size = 0;
    s_batch->size = 0;
    pow_batch->size = 0;

    for (i = 0; i < BATCH_SIZE; i++)
    {
      mpz_urandomm(pow, rand_state, modulus);
      mpz_powm(test_num, base, pow, modulus);
      eea_bounded_mpz(test_num, modulus, test_num_r, test_num_s, temp->elems);
      mpz_swap(r_batch->elems[i], test_num_r);
      mpz_swap(s_batch->elems[i], test_num_s);
      mpz_swap(pow_batch->elems[i], pow);
    }
    r_smooth_bools = batch_smooth_parts(z, primes_product_tree, r_batch, BATCH_SIZE); 

    for (i = 0; i < BATCH_SIZE; i++)
    {
      if (r_smooth_bools[i] == 1)
      {
        s_smooth_bools = batch_smooth_parts(z, primes_product_tree, s_batch, BATCH_SIZE);
        // both r and s values are smooth
        if (s_smooth_bools[i] == 1)
        {
          char r_val[100];
          char s_val[100];
          char pow_val[100];
          mpz_get_str(r_val, 10, r_batch->elems[i]);
          mpz_get_str(s_val, 10, s_batch-> elems[i]);
          mpz_get_str(pow_val, 10, pow_batch->elems[i]);
          fprintf(output, "%s:", r_val);
          fprintf(output, "%s:", s_val);
          fprintf(output, "%s:", pow_val);
          fprintf(output, "\n");
          diff = clock() - start;
          msec = diff * 1000 / CLOCKS_PER_SEC; 
          printf("--- !!!FOUND ONE!!! ---\n");
          printf("--- %lu seconds ---\n", msec);
        }
      }
    }
  }

  mpz_clear(a);
  mpz_clear(z);
  mpz_clear(base);
  mpz_clear(modulus);
  mpz_clear(pow);
  mpz_clear(test_num);
  mpz_clear(test_num_r);
  mpz_clear(test_num_s);
  vec_mpz_free(temp);
  vec_mpz_free(r_batch);
  vec_mpz_free(s_batch);
  vec_mpz_free(primes);
  vec_mpz_free(primes_product_tree);
}

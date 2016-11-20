#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"
#include "gmp.h"
#include "pthread.h"
#include "time.h"
#include "vec_mpz.h"
#include "vec_uint64.h"
#include "crypto_utils.h"

struct shared_data {
  mpz_t modulus;
  mpz_t base;
  struct vec_uint64_t* primes;
};

void *find_smooth_function( void * );

int main(int argc, char ** argv) {

  int i, n, numThreads = 0, c = 0;

  uint64_t a;

  FILE *fp;

  struct shared_data data;

  if (argc != 2) {
    printf("Please pass number of threads as cli option.\n");
    return 1;
  }

  numThreads = atoi(argv[1]);

  printf("Loading prime factor base\n");

  data.primes = vec_uint64_init();

  mpz_init(data.modulus);
  mpz_set_str(
    data.modulus,
    "3217014639198601771090467299986349868436393574029172456674199",
    10
  );

  mpz_init(data.base);
  mpz_set_str(
    data.base,
    "1483470583978676987274572113759546747043713044116589831659689",
    10
  );

  fp = fopen("./primes.txt", "r");

  while (fscanf(fp, "%llu", &a) > 0) {
    vec_uint64_insert(data.primes, a);
  }

  fclose(fp);

  printf(
    "Done loading prime factor base: %llu primes loaded.\n",
    data.primes->size
  );

  printf("Starting with %d threads.\n", numThreads);

  pthread_t* threads = malloc(numThreads * sizeof(pthread_t));

  for (i = 0; i < numThreads; ++i) {
    pthread_create(&threads[0], NULL, find_smooth_function, &data);
  }

  for (i = 0; i < numThreads; ++i) {
    void* status;
    pthread_join(threads[i], &status);
  }
  mpz_clear(data.modulus);
  vec_uint64_free(data.primes);
}

void *find_smooth_function( void *ptr ) {
  int i, msec, attempts = 0;

  clock_t start, diff;

  struct shared_data * data = (struct shared_data*)ptr;

  gmp_randstate_t rand_state;
  gmp_randinit_mt(rand_state);
  gmp_randseed_ui(rand_state, pthread_self());

  uint64_t largest = data->primes->elems[data->primes->size - 1];

  mpz_t pow, test_num, factor_num;
  mpz_inits(pow, test_num, factor_num, 0);

  char filename[128];

  sprintf(filename, "./smooth_thread_%u", pthread_self());

  FILE *smoothFile;

  smoothFile = fopen(filename, "w+");

  uint64_t *primes = data->primes->elems;
  int number_of_primes = data->primes->size;
  uint64_t prime;
  start = clock();


  while (1) {

    mpz_urandomm(pow, rand_state, data->modulus);
    mpz_powm(test_num, data->base, pow, data->modulus);
    mpz_set(factor_num, test_num);

    if (attempts % 1000 == 0) {
      diff = clock() - start;
      msec = diff * 1000 / CLOCKS_PER_SEC;
      printf("Attempt count: %d Time: %d\n", attempts, msec);
      start = clock();
    }
    // Some tests pass at this point, break it down (trial division)

    struct vec_uint64_t *factors = vec_uint64_init();
    struct vec_uint64_t *powers = vec_uint64_init();
    int power;

    for (i = 0; i < number_of_primes; i++) {
      power = 0;
      prime = primes[i];
      while (mpz_divisible_ui_p(factor_num, prime) != 0) {
        mpz_divexact_ui(factor_num, factor_num, prime);
        power += 1;
      }
      if (power != 0) {
        vec_uint64_insert(factors, prime);
        vec_uint64_insert(powers, power);
      }
    }

    if (mpz_cmp_ui(factor_num, 1) == 0) {
      gmp_printf("Fully factored: %Zd\n", test_num);
      char powStr[100];

      mpz_get_str(powStr, 10, pow);

      fprintf(smoothFile, "%s:", powStr);

      for (i = 0; i < factors->size; ++i) {
        fprintf(smoothFile, "(%lu, %lu)", factors->elems[i], powers->elems[i]);
      }

      fprintf(smoothFile, "\n");
      fflush(smoothFile);
    }

    vec_uint64_free(factors);
    vec_uint64_free(powers);
    attempts += 1;
  }


}

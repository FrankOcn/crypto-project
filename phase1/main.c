#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"
#include "gmp.h"
#include "pthread.h"
#include "time.h"
#include "vec_mpz.h"
#include "vec_uint64.h"
#include "crypto_utils.h"
#include "hashtable_int64.h"

struct shared_data {
  mpz_t modulus;
  mpz_t base;
  struct vec_mpz_t* primes;
  struct vec_mpz_t* prime_product_tree;
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

  data.primes = vec_mpz_init();
  data.prime_product_tree = vec_mpz_init();

  mpz_init(data.modulus);
  mpz_set_str(
    data.modulus,
    "3217014639198601771090467299986349868436393574029172456674199",
    10
  );

  mpz_init(data.base);
  mpz_set_str(
    data.base,
    "983323641496181394191925968825033998345778072524604244950047",
    10
  );

  fp = fopen("./primes.txt", "r");

  mpz_t tmp;
  mpz_init(tmp);

  while (fscanf(fp, "%llu", &a) > 0) {
    mpz_set_ui(tmp, a);
    vec_mpz_insert(data.primes, tmp);
  }

  mpz_clear(tmp);

  data.prime_product_tree = vec_mpz_product_tree(data.primes);


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
  mpz_clear(data.base);
  vec_mpz_free(data.primes);
}

void *find_smooth_function( void *ptr ) {
  int i, msec, failed, attempts = 0;

  clock_t start, diff;

  struct shared_data * data = (struct shared_data*)ptr;

  gmp_randstate_t rand_state;
  gmp_randinit_mt(rand_state);
  gmp_randseed_ui(rand_state, pthread_self());

  mpz_t largest_prime;
  mpz_init_set(largest_prime, data->primes->elems[0]);

  // Initialize some test mpz_t for temp use in crypto util functions.
  struct vec_mpz_t* tmp = vec_mpz_init();

  mpz_t pow, test_num, test_num_r, test_num_s,
    pollard_num, pollard_factor, factor_num;

  mpz_inits(pow, test_num, test_num_r, test_num_s,
    pollard_num, pollard_factor, factor_num, 0);

  char filename[128];

  sprintf(filename, "./smooth_thread_%u", pthread_self());

  FILE *smoothFile;

  smoothFile = fopen(filename, "w+");

  mpz_t *primes = data->primes->elems;
  int number_of_primes = data->primes->size;
  mpz_t prime;
  mpz_init(prime);
  start = clock();


  while (1) {

    mpz_urandomm(pow, rand_state, data->modulus);
    mpz_powm(test_num, data->base, pow, data->modulus);

    if (attempts % 1000 == 0) {
      diff = clock() - start;
      msec = diff * 1000 / CLOCKS_PER_SEC;
      printf("Attempt count: %d Time: %d\n", attempts, msec);
      start = clock();
    }

    eea_bounded_mpz( test_num, data->modulus, test_num_r, test_num_s, tmp->elems);

    mpz_set(pollard_num, test_num_s);
    if (mpz_sgn(pollard_num) < 0) {
      mpz_abs(pollard_num, pollard_num);
    }
    if (
      miller_rabin(test_num_r, 10, rand_state, tmp->elems) == 1 ||
      miller_rabin(pollard_num, 10, rand_state, tmp->elems) == 1
    ) {
      attempts += 1;
      continue;
    }

    mpz_set(pollard_num, test_num_r);

    /*for (i = 0; i < 150; ++i) {
      mpz_set(prime, primes[i]);
      while (mpz_divisible_p(pollard_num, prime) != 0) {
        mpz_divexact(pollard_num, pollard_num, prime);
      }
    }*/

    failed = 0;

    while (failed == 0 && mpz_cmp(pollard_num, largest_prime) > 0) {
      if (pollard_rho(pollard_factor, pollard_num, 4000, tmp->elems) == 1) {
        mpz_divexact(pollard_num, pollard_num, pollard_factor);
        if (mpz_cmp(pollard_num, largest_prime) > 0 && miller_rabin(pollard_num, 10, rand_state, tmp->elems) == 1) {
          failed = 1;
        }
      } else {
        failed = 1;
      }
    }

    if (failed == 1) {
      attempts += 1;
      continue;
    }

    gmp_printf("Passed (r): %Zd\n", test_num_r);
    gmp_printf("Testing (s): %Zd\n", test_num_s);

    mpz_set(pollard_num, test_num_s);
    if (mpz_sgn(pollard_num) < 0) {
      mpz_abs(pollard_num, pollard_num);
    }
    /*for (i = 0; i < 150; ++i) {
      mpz_set(prime, primes[i]);
      while (mpz_divisible_p(pollard_num, prime) != 0) {
        mpz_divexact(pollard_num, pollard_num, prime);
      }
    }*/

    failed = 0;
    while (failed == 0 && mpz_cmp(pollard_num, largest_prime) > 0) {
      if (pollard_rho(pollard_factor, pollard_num, 4000, tmp->elems) == 1) {
        mpz_divexact(pollard_num, pollard_num, pollard_factor);
      } else {
        failed = 1;
      }
    }

    if (failed == 1) {
      attempts += 1;
      continue;
    }

    gmp_printf("Attempting\n");
    mpz_set(factor_num, test_num_r);
    struct hashtable_int64_t *table = hashtable_int64_init();
    int power;

    for (i = 0; i < number_of_primes; i++) {
      power = 0;
      mpz_set(prime, primes[i]);
      while (mpz_divisible_p(factor_num, prime) != 0) {
        mpz_divexact(factor_num, factor_num, prime);
        power += 1;
      }
      if (power != 0) {
        hashtable_int64_insert(table, mpz_get_ui(prime), power);
      }
    }

    if (mpz_cmp_ui(factor_num, 1) != 0) {
      attempts += 1;
      hashtable_int64_free(table);
      continue;
    }

    gmp_printf("Fully factored (r): %Zd\n", test_num_r);
    mpz_set(factor_num, test_num_s);

    if (mpz_sgn(factor_num) < 0) {
      hashtable_int64_insert(table, -1, 1);
      mpz_abs(factor_num, factor_num);
    }

    for (i = 0; i < number_of_primes; i++) {
      power = 0;
      mpz_set(prime, primes[i]);
      while (mpz_divisible_p(factor_num, prime) != 0) {
        mpz_divexact(factor_num, factor_num, prime);
        power += 1;
      }
      if (power != 0) {
        int64_t *saved = hashtable_int64_retrieve(table, mpz_get_ui(prime));
        if (saved == 0) {
          hashtable_int64_insert(table, mpz_get_ui(prime), -power);
        } else {
          *saved = *saved - power;
        }
      }
    }

    if (mpz_cmp_ui(factor_num, 1) == 0) {
      gmp_printf("Fully factored (s): %Zd\n", test_num_s);
      printf("------------------FUCK YES----------------\n");
      char powStr[100];

      mpz_get_str(powStr, 10, pow);

      fprintf(smoothFile, "%s:", powStr);

      for (i = 0; i < table->capacity; ++i) {
        if (table->elems[i] != 0) {
          fprintf(smoothFile, "(%d, %d)", table->elems[i]->key, table->elems[i]->value);
        }
      }

      fprintf(smoothFile, "\n");
      fflush(smoothFile);
    }

    hashtable_int64_free(table);
    attempts += 1;
  }


}

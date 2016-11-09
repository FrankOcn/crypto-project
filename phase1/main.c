#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>
#include <pthread.h>
#include <time.h>

#define MOD_ELEM_0 135291589527
#define MOD_ELEM_1 18446744073705357312U
#define MOD_ELEM_2 9223371487099224063
#define MOD_ELEM_3 512

struct thread_data {
  mpz_t mP;
  struct vec_t* primes;
  mpz_t primeProd;
  mpz_t fastFilterMax;
};

struct vec_t {
  uint64_t * elems;
  int size;
  int capacity;
};

void insert(struct vec_t* vec, uint64_t elem) {
  int i;
  if (vec->size == vec->capacity) {
    vec->capacity = vec->capacity * 2;
    uint64_t* newElems = (uint64_t *)malloc(vec->capacity * sizeof(uint64_t));
    for (i = 0; i < vec->size; ++i) {
      newElems[i] = vec->elems[i];
    }
    free(vec->elems);
    vec->elems = newElems;
  }
  vec->elems[vec->size] = elem;
  vec->size++;
}

void pollard_rho_g(mpz_t x, mpz_t c) {
  mpz_powm_ui(x, x, 2, c);
  mpz_add_ui(x, x, 1);
  mpz_mod(x, x, c);
}

int pollard_rho(mpz_t d, mpz_t n) {
  mpz_t x, y, a;
  mpz_init(x);
  mpz_init(y);
  mpz_init(a);
  mpz_set_ui(x, 2);
  mpz_set_ui(y, 2);
  mpz_set_ui(d, 1);
  int iterations = 0;
  while (mpz_cmp_ui(d, 1) == 0 && iterations < 5000) {
    pollard_rho_g(x, n);
    pollard_rho_g(y, n);
    pollard_rho_g(y, n);
    mpz_sub(a, x, y);
    mpz_abs(a, a);
    mpz_gcd(d, a, n);
    iterations += 1;
  }
  if (iterations == 5000) {
    return -1;
  }

  return 1;
}


int miller_rabin(mpz_t n, int testLimit, gmp_randstate_t rand_state) {
  int i, j, r;
  mpz_t d, x, n_1, n_4, a;
  mpz_init(d);
  mpz_init(x);
  mpz_init(n_1);
  mpz_init(n_4);
  mpz_init(a);
  mpz_set(d, n);
  mpz_set(n_1, n);
  mpz_set(n_4, n);
  mpz_sub_ui(d, d, 1);
  mpz_sub_ui(n_1, n_1, 1);
  mpz_sub_ui(n_4, n_4, 4);

  while (mpz_even_p(d) != 0) {
    r += 1;
    mpz_tdiv_q_2exp(d, d, 1);
  }

  int continueOuter;
  for (i = 0; i < testLimit; ++i) {
    continueOuter = 0;
    mpz_urandomm(a, rand_state, n_4);
    mpz_add_ui(a, a, 2);
    mpz_powm(x, a, d, n);
    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_1) == 0) {
      continue;
    }
    for (j = 0; j < r - 1; ++j) {
      mpz_powm_ui(x, x, 2, n);
      if (mpz_cmp_ui(x, 1) == 0) {
        return -1;
      }
      if (mpz_cmp(x, n_1) == 1) {
        continueOuter = 1;
        break;
      }
    }
    if (continueOuter == 0) {
      return -1;
    }
  }
  return 1;
}

int fast_filter(mpz_t n, mpz_t max, gmp_randstate_t rand_state) {
  mpz_t f, n_1;
  mpz_init(f);
  mpz_init(n_1);
  mpz_set(n_1, n);

  while (mpz_cmp(n_1, max) > 0) {
    if (pollard_rho(f, n_1) != 1) {
      mpz_clear(f);
      mpz_clear(n_1);
      return -1;
    }
    //gmp_printf("%Zd plrd\n", n);
    mpz_divexact(n_1, n_1, f);
  }
  mpz_clear(f);
  mpz_clear(n_1);
  return 1;
}

struct vec_t* initVec() {
  struct vec_t* newVec = (struct vec_t*)malloc(sizeof(struct vec_t));
  newVec->capacity = 32;
  newVec->size = 0;
  newVec->elems = (uint64_t *)malloc(newVec->capacity * sizeof(uint64_t));
  return newVec;
}

void clear(struct vec_t* vec) {
  free(vec->elems);
  vec->capacity = 64;
  vec->size = 0;
  vec->elems = (uint64_t *)malloc(vec->capacity * sizeof(uint64_t));
}


void *find_smooth_function( void *ptr ) {
  int i;

  clock_t start, diff;

  gmp_randstate_t rand_state;
  struct thread_data * td;
  td = (struct thread_data*)ptr;
  mpz_t testNum, pow, phi, base, modResult, gcd, softTest, aBase, poll;

  mpz_t quotient, remainder;

  gmp_randinit_default(rand_state);
  gmp_randseed_ui(rand_state, pthread_self());
  mpz_init(testNum);
  mpz_init(softTest);
  mpz_init(modResult);
  mpz_init(aBase);
  mpz_init(pow);
  mpz_init(poll);
  mpz_init(phi);
  mpz_init(gcd);
  mpz_init(base);
  mpz_set_ui(base, 5);
  mpz_set(phi, td->mP);
  mpz_sub_ui(phi, phi, 1);
  mpz_divexact_ui(phi, phi, 2);


  uint64_t smallPrimeProd = 1;

  char filename[128];

  sprintf(filename, "./smooth_thread_%u", pthread_self());

  FILE *smoothFile;

  smoothFile = fopen(filename, "w+");

  struct vec_t *factors = initVec();
  struct vec_t *powers = initVec();
  uint64_t prime;
  int attemptCount = 0;

  mpz_t f, n_1;
  mpz_init(f);
  mpz_init(n_1);

  while (1) {

    factors->size = 0;
    powers->size = 0;

    mpz_urandomm(pow, rand_state, td->mP);

    mpz_mod(pow, pow, phi);

    mpz_powm(testNum, base, pow, td->mP);

    if (miller_rabin(testNum, 100, rand_state) == 1) {
      continue;
    }
    mpz_set(n_1, testNum);
    int failed = 0;

    /*while (failed == 0 && mpz_cmp(n_1, td->fastFilterMax) > 0) {
      if (pollard_rho(f, n_1) != 1) {
        failed = 1;
      }
      //gmp_printf("%Zd plrd\n", n);
      mpz_divexact(n_1, n_1, f);
    }
*/

    //printf("running filter\n");
    if (failed != 1) {
      //printf("filter passed!\n");
      int power = 0;
      for (i = 0; i < td->primes->size; ++i) {
        power = 0;
        prime = td->primes->elems[i];
        while (mpz_divisible_ui_p(testNum, prime) != 0) {
          mpz_divexact_ui(testNum, testNum, prime);
          power += 1;
        }
        if (power != 0) {
          insert(factors, prime);
          insert(powers, power);
        }
      }

      /*for (i = 0; i < factors->size; ++i) {
        printf("%d: %d\n", factors->elems[i], powers->elems[i]);
      }*/
      if (mpz_cmp_ui(testNum, 1) == 0) {
        char powStr[100];

        mpz_get_str(powStr, 10, pow);

        fprintf(smoothFile, "%s:", powStr);

        for (i = 0; i < factors->size; ++i) {
          fprintf(smoothFile, "(%lu, %lu)", factors->elems[i], powers->elems[i]);
        }

        fprintf(smoothFile, "\n");
      }
    } else {
      //printf("filter failed!\n");
    }

    attemptCount += 1;
    if (attemptCount % 100 == 0) {
      printf("Attempts: %lu\n", attemptCount);
    }
  }
  mpz_clear(testNum);
  mpz_clear(pow);
  mpz_clear(phi);
  mpz_clear(f);
  mpz_clear(n_1);

}

int main(int argc, char ** argv) {

  int i;

  uint64_t* modulus = malloc(sizeof(uint64_t) * 4);
  modulus[0] = MOD_ELEM_0;
  modulus[1] = MOD_ELEM_1;
  modulus[2] = MOD_ELEM_2;
  modulus[3] = MOD_ELEM_3;

  printf("Loading prime factor base\n");

  FILE *fp;

  struct thread_data td;
  td.primes = initVec();
  mpz_t mQrTest, mPTest, mInverse, primePow1, primePow2;
  mpz_init(td.mP);
  mpz_init(td.primeProd);
  mpz_init(td.fastFilterMax);
  mpz_init(mQrTest);
  mpz_init(mPTest);
  mpz_init(mInverse);
  mpz_init(primePow1);
  mpz_init(primePow2);
  mpz_set_ui(td.primeProd, 1);
  mpz_import(td.mP, 4, -1, sizeof(uint64_t), 0, 0, modulus);
  mpz_set(mQrTest, td.mP);
  mpz_add_ui(mQrTest, mQrTest, 1);
  mpz_divexact_ui(mQrTest, mQrTest, 4);

  fp = fopen("./primes.txt", "r");
  int c = 0;
  uint64_t a;
  while (c != -1) {
    a = 0;
    do {
      c = fgetc(fp);
      if (c >= '0' && c <= '9') {
        a *= 10;
        uint64_t n = c - '0';
        a += n;
      }
    } while (c >= '0' && c <= '9');
    if (a != 0) {
      /*mpz_set_ui(mPTest, a);
      mpz_powm(mPTest, mPTest, mQrTest, td.mP);
      mpz_powm_ui(mPTest, mPTest, 2, td.mP);
      if (mpz_cmp_d(mPTest, a) == 0) {*/
        insert(td.primes, a);
      //}
    }
  }

  mpz_ui_pow_ui(td.fastFilterMax, td.primes->elems[td.primes->size -1], 4);

  fclose(fp);

  printf("Done loading prime factor base: %llu primes loaded.\n", td.primes->size);

  printf("%Zd\n", td.primeProd);

  int numThreads = 0;

  for (i = 0; argv[1][i] != '\0'; ++i) {

    if (argv[1][i] >= '0' && argv[1][i] <= '9') {
      numThreads *= 10;
      int n = argv[1][i] - '0';
      numThreads += n;
    }
  }

  printf("Starting with %d threads.\n", numThreads);

  pthread_t* threads = malloc(numThreads * sizeof(pthread_t));

  for (i = 0; i < numThreads; ++i) {
    pthread_create(&threads[0], NULL, find_smooth_function, &td);
  }

  for (i = 0; i < numThreads; ++i) {
    void* status;
    pthread_join(threads[i], &status);
  }
  mpz_clear(mQrTest);
  mpz_clear(mPTest);
  clear(td.primes);
}

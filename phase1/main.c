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
  struct vec_t* smooth;
  mpz_t primeProd;
  mpz_t pollardTest;
  mpz_t pollardTestPrime;
  uint64_t smallNonQuad;
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

int miller_rabin(mpz_t n, int testLimit, gmp_randstate_t rand_state) {
  int i, j, r;
  // There's a bug if n <= 4, so putting the answers here for good measure.
  if (mpz_cmp_ui(n, 1) == 0 || mpz_cmp_ui(n, 4) == 0) {
    return -1;
  }
  if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) {
    return 1;
  }
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

  if (mpz_cmp_ui(d, 1) == 0) {
    mpz_clear(d);
    mpz_clear(x);
    mpz_clear(n_1);
    mpz_clear(n_4);
    mpz_clear(a);
    return -1;
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
        mpz_clear(d);
        mpz_clear(x);
        mpz_clear(n_1);
        mpz_clear(n_4);
        mpz_clear(a);
        return -1;
      }
      if (mpz_cmp(x, n_1) == 1) {
        continueOuter = 1;
        break;
      }
    }
    if (continueOuter == 0) {
      mpz_clear(d);
      mpz_clear(x);
      mpz_clear(n_1);
      mpz_clear(n_4);
      mpz_clear(a);
      return -1;
    }
  }
  mpz_clear(d);
  mpz_clear(x);
  mpz_clear(n_1);
  mpz_clear(n_4);
  mpz_clear(a);
  return 1;
}

int pollard_rho(mpz_t d, mpz_t n, int limit) {
  int i;
  mpz_t x, y, diff;
  mpz_init(x);
  mpz_init(y);
  mpz_init(diff);
  mpz_set_ui(x, 2);
  mpz_set_ui(y, 2);
  mpz_set_ui(d, 1);
  for (i = 0; i < limit && mpz_cmp_ui(d, 1) == 0; ++i) {
    mpz_powm_ui(x, x, 2, n);
    mpz_add_ui(x, x, 1);
    mpz_powm_ui(y, y, 2, n);
    mpz_add_ui(y, y, 1);
    mpz_powm_ui(y, y, 2, n);
    mpz_add_ui(y, y, 1);
    mpz_sub(diff, x, y);
    mpz_abs(diff, diff);
    mpz_gcd(d, diff, n);
  }
  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(diff);
  return i == limit ? -1 : 1;
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

  int msec = 0;

  gmp_randstate_t rand_state;
  struct thread_data * td;
  td = (struct thread_data*)ptr;
  mpz_t testNum, pow, phi, base, modResult, gcd, softTest, aBase, poll;

  mpz_t quotient, remainder, num;

  gmp_randinit_mt(rand_state);
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
  mpz_init(quotient);
  mpz_init(remainder);
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
  int power;

  start = clock();

  while (1) {

    factors->size = 0;
    powers->size = 0;

    mpz_urandomm(pow, rand_state, phi);

    mpz_powm(testNum, base, pow, td->mP);

    //mpz_set_ui(testNum, 676250);

    /*mpz_gcd_ui(gcd, testNum, td->smallNonQuad);

    if (mpz_cmp_ui(gcd, 1) != 0) {
      attemptCount += 1;
      if (attemptCount % 100 == 0) {
        printf("Small Quad Attempts: %lu\n", attemptCount);
      }
      continue;
    }*/

    //mpz_set(softTest, testNum);
    // Factor out small values
    /*for (int i = 0; i < 500; ++i) {
      prime = td->primes->elems[i];
      mpz_set_ui(remainder, 0);
      mpz_set(quotient, softTest);
      while (mpz_cmp_ui(remainder, 0) == 0) {
        mpz_set(softTest, quotient);
        mpz_tdiv_qr_ui(quotient, remainder, quotient, prime);
      }
    }*/

    /*mpz_t d, q;
    mpz_init(d);
    mpz_init(q);


    if (pollard_rho(d, softTest, 4000) == -1) {
      attemptCount += 1;
      if (attemptCount % 100 == 0) {
        gmp_printf("Pollard rho attempts: %Zd, %Zd\n", softTest, d);
        printf("Attempts: %lu\n", attemptCount);
      }
      continue;
    }

    mpz_divexact(q, testNum, d);

    if (miller_rabin(q, 100, rand_state) == 1) {
      attemptCount += 1;
      if (attemptCount % 100 == 0) {
        gmp_printf("Pollard rho prime attempts: %Zd, %Zd\n", softTest, d);
        printf("Attempts: %lu\n", attemptCount);
      }
      continue;
    }*/

    mpz_set(softTest, testNum);

    // Factor out powers of 2.
    while (mpz_even_p(softTest) != 0) {
      mpz_tdiv_q_2exp(softTest, softTest, 1);
    }
    mpz_gcd(gcd, testNum, td->pollardTest);

    mpz_divexact(softTest, testNum, gcd);

    int largest = td->primes->elems[td->primes->size - 1];
    largest *= largest;
    int failed = 0, iterations = 0;

    while (mpz_cmp_ui(softTest, largest) > 0 && miller_rabin(softTest, 10, rand_state) != 1) {
      // Pollard's p-1
      iterations += 1;
      mpz_powm(modResult, td->pollardTestPrime, td->pollardTest, softTest);
      mpz_sub_ui(modResult, modResult, 1);
      mpz_gcd(gcd, softTest, modResult);

      // We're not adjusting any parameters. If we find nothing, throw it away.
      if (mpz_cmp_ui(gcd, 1) == 0 || mpz_cmp(gcd, softTest) == 0) {
        //printf("Could not factor. Whatever.\n");
        failed = 1;
        break;
      }

      // If the factor we found is a prime larger than the largest prime in our
      // list, we're done
      if (miller_rabin(gcd, 10, rand_state) == 1 && mpz_cmp_ui(gcd, largest) > 0) {
        //printf("Prime!.\n");
        failed = 1;
        break;
      }

      mpz_divexact(softTest, softTest, gcd);

    }

    if (mpz_cmp_ui(softTest, largest) > 0 || failed == 1) {
      attemptCount += 1;
      if (attemptCount % 10000 == 0) {
        diff = clock() - start;
        msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Attempts: %lu, Time: %d\n, Iterations: %d\n", attemptCount, msec, iterations);
        start = clock();
      }
      continue;
    }

    int found = 0, j;
    for (i = 0; i < td->primes->size; i++ ) {
      power = 0;
      prime = td->primes->elems[i];
      mpz_set_ui(remainder, 0);
      mpz_set(quotient, testNum);
      while (mpz_cmp_ui(remainder, 0) == 0) {
        mpz_tdiv_qr_ui(quotient, remainder, quotient, prime);
        if (mpz_cmp_ui(remainder, 0) == 0) {
          mpz_divexact_ui(testNum, testNum, prime);
          power += 1;
        }
      }
      if (power != 0) {
        found = 1;
        insert(factors, prime);
        insert(powers, power);
      }

    }

    if (mpz_cmp_ui(testNum, 1) == 0) {
      char powStr[100];

      mpz_get_str(powStr, 10, pow);

      fprintf(smoothFile, "%s:", powStr);

      for (i = 0; i < factors->size; ++i) {
        fprintf(smoothFile, "(%lu, %lu)", factors->elems[i], powers->elems[i]);
      }

      fprintf(smoothFile, "\n");
    }

    attemptCount += 1;
    if (attemptCount % 1000 == 0) {
      diff = clock() - start;
      msec = diff * 10000 / CLOCKS_PER_SEC;
      printf("Attempts: %lu, Time: %d\n", attemptCount, msec);
      start = clock();
    }
  }
  mpz_clear(testNum);
  mpz_clear(pow);
  mpz_clear(phi);
  mpz_clear(f);
  mpz_clear(n_1);

}

int main(int argc, char ** argv) {

  int i, exp;

  uint64_t* modulus = malloc(sizeof(uint64_t) * 4);
  modulus[0] = MOD_ELEM_0;
  modulus[1] = MOD_ELEM_1;
  modulus[2] = MOD_ELEM_2;
  modulus[3] = MOD_ELEM_3;

  printf("Loading prime factor base\n");

  FILE *fp;

  struct thread_data td;
  td.primes = initVec();
  td.smooth = initVec();
  mpz_t mQrTest, mPTest, primePow1, primePow2, aBase;
  mpz_init(td.mP);
  mpz_init(td.primeProd);
  mpz_init(mQrTest);
  mpz_init(mPTest);
  mpz_init(aBase);
  mpz_init(primePow1);
  mpz_init(primePow2);
  mpz_init(td.pollardTest);
  mpz_init(td.pollardTestPrime);
  mpz_set_ui(td.pollardTest, 1);
  mpz_set_ui(td.primeProd, 1);
  mpz_import(td.mP, 4, -1, sizeof(uint64_t), 0, 0, modulus);
  mpz_set_ui(td.pollardTestPrime, 2);
  mpz_set(mQrTest, td.mP);
  mpz_add_ui(mQrTest, mQrTest, 1);
  mpz_divexact_ui(mQrTest, mQrTest, 4);

  fp = fopen("./primes.txt", "r");
  int c = 0;
  uint64_t a;
  td.smallNonQuad = 1;
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
      mpz_set_ui(mPTest, a);
      mpz_powm(mPTest, mPTest, mQrTest, td.mP);
      mpz_powm_ui(mPTest, mPTest, 2, td.mP);
      if (a < 6500) {
        int b = a;
        while (b * a < 6500) {
          b *= a;
        }
        mpz_mul_ui(td.pollardTest, td.pollardTest, a);
      }
      if (mpz_cmp_d(mPTest, a) == 0) {
        insert(td.primes, a);

        /*mpz_set_ui(primePow1, a);
        mpz_set(primePow2, primePow1);
        while (mpz_cmp(primePow1, td.mP) < 0) {
          mpz_swap(primePow2, primePow1);
          mpz_pow_ui(primePow1, primePow1, 2);
        }
        mpz_mul_ui(td.pollardTest, td.pollardTest, a);
        printf("%d\n", a);*/
      } else if (td.smallNonQuad < UINT64_MAX) {
        td.smallNonQuad *= a;
      }
    }
  }

  /*for (i = 0; i < 100000; ++i) {
    mpz_set_ui(aBase, td.primes->elems[i]);
    // Bug with mpz_sizeinbase. It segfaults if the base is greater than 2333.
    exp = i < 341 ? mpz_sizeinbase(td.mP, td.primes->elems[i]) : 15;
    mpz_powm_ui(aBase, aBase, exp, td.mP);
    mpz_mul(td.pollardTest, td.pollardTest, aBase);
  }*/

  fclose(fp);

  printf("Done loading prime factor base: %llu primes loaded.\n", td.primes->size);

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

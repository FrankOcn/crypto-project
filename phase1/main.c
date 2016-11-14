#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>
#include <pthread.h>
#include <time.h>

struct thread_data {
  mpz_t mP;
  struct vec_t* primes;
  struct vec_t* smooth;
  mpz_t primeProd;
  mpz_t pollardTest;
  mpz_t pollardTestPrime;
  mpz_t mQrTest;
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
  mpz_t x, y, diff, a;
  mpz_init(x);
  mpz_init(y);
  mpz_init(diff);
  mpz_init(a);
  mpz_set_ui(x, 2);
  mpz_set_ui(y, 2);
  mpz_set_ui(d, 1);
  mpz_set_ui(a, 1);
  for (i = 0; i < limit && mpz_cmp_ui(d, 1) == 0; ++i) {
    mpz_powm_ui(x, x, 2, n);
    mpz_sub_ui(x, x, 1);
    mpz_powm_ui(y, x, 2, n);
    mpz_sub_ui(y, y, 1);
    mpz_sub(diff, x, y);
    mpz_abs(diff, diff);

    mpz_mul(a, a, diff);
  }
  mpz_gcd(d, a, n);
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

  int found = 0, j;
  clock_t start, diff;

  int msec = 0;

  gmp_randstate_t rand_state;
  struct thread_data * td;
  td = (struct thread_data*)ptr;
  mpz_t testNum, pow, phi, base, modResult, gcd, softTest, aBase, poll, pollardTestPrime;

  mpz_t limit, factorNum;

  mpz_t a, a2, b, d;

  mpz_t x, y, diffm;
  mpz_init(x);
  mpz_init(y);
  mpz_init(diffm);
  mpz_init(a);
  mpz_init(b);
  mpz_init(d);

  mpz_t quotient, remainder, num;

  gmp_randinit_mt(rand_state);
  gmp_randseed_ui(rand_state, pthread_self());
  mpz_init(testNum);
  mpz_init(softTest);
  mpz_init(factorNum);
  mpz_init(modResult);
  mpz_init(aBase);
  mpz_init(pow);
  mpz_init(poll);
  mpz_init(phi);
  mpz_init(gcd);
  mpz_init(limit);
  mpz_init(pollardTestPrime);
  mpz_init(base);
  mpz_init(quotient);
  mpz_init(remainder);
  mpz_set_str(base, "1483470583978676987274572113759546747043713044116589831659689", 10);
  mpz_set_ui(pollardTestPrime, 2);
  mpz_set(phi, td->mP);
  mpz_sub_ui(phi, phi, 1);
  mpz_divexact_ui(phi, phi, 2);

  int largest = td->primes->elems[td->primes->size - 1];
  mpz_set_ui(limit, largest);
  mpz_pow_ui(limit, limit, 2);


  uint64_t smallPrimeProd = 1;

  char filename[128];

  sprintf(filename, "./smooth_thread_%u", pthread_self());

  FILE *smoothFile;

  smoothFile = fopen(filename, "w+");

  struct vec_t *factors = initVec();
  struct vec_t *powers = initVec();
  uint64_t *primes = td->primes->elems;
  int numberOfPrimes = td->primes->size;
  uint64_t prime;
  int attemptCount = 0;

  mpz_t f, n_1;
  mpz_init(f);
  mpz_init(n_1);
  int power;

  start = clock();


  mpz_urandomm(pow, rand_state, td->mP);
  mpz_powm(testNum, base, pow, td->mP);

  while (1) {

    int failed = 0, iterations = 0;
    found = 0;
    factors->size = 0;
    powers->size = 0;
    //mpz_mul(testNum, testNum, base);
    mpz_powm_ui(testNum, testNum, 3, td->mP);

    if (mpz_divisible_ui_p(testNum, 5) != 0) {
      attemptCount += 1;
      if (attemptCount % 10000 == 0) {
        diff = clock() - start;
        msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Five filter: Attempts: %lu, Time: %d\n", attemptCount, msec);
        start = clock();
      }
      continue;
    }


    mpz_set(softTest, testNum);

    // Factor out powers of 2.
    while (mpz_even_p(softTest) != 0) {
      mpz_tdiv_q_2exp(softTest, softTest, 1);
    }
    mpz_sub_ui(aBase, softTest, 1);

    mpz_powm(a, pollardTestPrime, aBase, softTest);

    if (mpz_cmp_ui(a, 1) == 0) {
      attemptCount += 1;
      if (attemptCount % 10000 == 0) {
        diff = clock() - start;
        msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Probably prime?: %lu, Time: %d\n", attemptCount, msec);
        start = clock();
      }
      continue;
    }

    for (i = 0; i < 100; ++i) {
      prime = td->primes->elems[i];
      while (mpz_divisible_ui_p(softTest, prime) != 0) {
        mpz_divexact_ui(softTest, softTest, prime);
      }
    }

    while (mpz_cmp(softTest, limit) > 0 && miller_rabin(softTest, 10, rand_state) != 1) {
      // Pollard's p-1
      iterations += 1;
      mpz_powm(modResult, pollardTestPrime, td->pollardTest, softTest);
      mpz_sub_ui(modResult, modResult, 1);
      mpz_gcd(gcd, softTest, modResult);

      if (attemptCount % 10000 == 0) {
        gmp_printf("p-1 factor: %Zd.\n", gcd);
      }
      if (mpz_cmp_ui(gcd, largest) > 0) {
        if (attemptCount % 10000 == 0) {
            gmp_printf("too large!\n", gcd);
          }
        failed = 1;
        break;
      }
      // We're not adjusting any parameters. If we find nothing, throw it away.
      if (mpz_cmp_ui(gcd, 1) == 0 || mpz_cmp(gcd, softTest) == 0) {
        if (attemptCount % 10000 == 0) {
          printf("could not factor.\n");
        }
        failed = 1;
        break;
      }

      mpz_divexact(softTest, softTest, gcd);
      if (attemptCount % 10000 == 0) {
        gmp_printf("next: %Zd.\n", softTest);
      }
    }
    /*mpz_sqrt(a, softTest);
    for (i = 0; i < 10000; ++i) {
      mpz_add_ui(a, a, 1);
      mpz_pow_ui(b, a, 2);
      mpz_sub(b, b, testNum);
      if (mpz_perfect_square_p(b)) {
        found = 1;
        break;
      }
    }

    if (found == 0) {
      attemptCount += 1;
      if (attemptCount % 10000 == 0) {
        diff = clock() - start;
        msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Cannot break Time: %d, attempts: %d\n", msec, attemptCount);
        start = clock();
      }
      continue;
    }*/

    /*mpz_gcd(gcd, softTest, td->primeProd);

    if (mpz_cmp_ui(gcd, 1) == 0) {
      attemptCount += 1;
      if (attemptCount % 10000 == 0) {
        diff = clock() - start;
        msec = diff * 1000 / CLOCKS_PER_SEC;
        gmp_printf("%Zd\n", softTest);
        printf("No factors Attempts: %lu, Time: %d, Iterations: %d\n", attemptCount, msec, iterations);
        start = clock();
      }
      continue;
    }*/

    if (mpz_cmp(softTest, limit) > 0 || failed == 1) {
      attemptCount += 1;
      if (attemptCount % 10000 == 0) {
        diff = clock() - start;
        msec = diff * 1000 / CLOCKS_PER_SEC;
        gmp_printf("%Zd\n", softTest);
        printf("Factor failure Attempts: %lu, Time: %d, Iterations: %d\n", attemptCount, msec, iterations);
        start = clock();
      }
      continue;
    }

    //gmp_printf("Factoring: %Zd.\n", softTest);
    //gmp_printf("Limit: %Zd.\n", limit);
    mpz_set(factorNum, softTest);
    //
    for (i = 1; i < numberOfPrimes; i++) {
      if (attemptCount % 10000 == 0) {
        gmp_printf("Factoring: %Zd.\n", testNum);
      }
      power = 0;
      prime = primes[i];
      while (mpz_divisible_ui_p(factorNum, 0) != 0) {
        mpz_divexact_ui(factorNum, factorNum, prime);
        power += 1;
      }
      if (power != 0) {
        found = 1;
        //printf("factor: %d^%d\n", prime, power);
        //gmp_printf("testNum: %Zd\n", testNum);
        insert(factors, prime);
        insert(powers, power);
      }

    }

    //mpz_gcd(gcd, td->primeProd, factorNum);

    //gmp_printf("Fully factored: %Zd\n", testNum);
    if (mpz_cmp_ui(factorNum, 1) == 0) {
      char powStr[100];

      mpz_get_str(powStr, 10, pow);

      fprintf(smoothFile, "%s:", powStr);

      for (i = 0; i < factors->size; ++i) {
        fprintf(smoothFile, "(%lu, %lu)", factors->elems[i], powers->elems[i]);
      }

      fprintf(smoothFile, "\n");
      fflush(smoothFile);
    }

    attemptCount += 1;
    if (attemptCount % 10000 == 0) {
      diff = clock() - start;
      msec = diff * 1000 / CLOCKS_PER_SEC;
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

  printf("Loading prime factor base\n");

  FILE *fp;

  struct thread_data td;
  td.primes = initVec();
  td.smooth = initVec();
  mpz_t mPTest, primePow1, primePow2, aBase;
  mpz_init(td.mP);
  mpz_init(td.primeProd);
  mpz_init(td.mQrTest);
  mpz_init(mPTest);
  mpz_init(aBase);
  mpz_init(primePow1);
  mpz_init(primePow2);
  mpz_init(td.pollardTest);
  mpz_init(td.pollardTestPrime);
  mpz_set_ui(td.pollardTest, 1);
  mpz_set_ui(td.primeProd, 1);
  mpz_set_str(td.mP, "3217014639198601771090467299986349868436393574029172456674199", 10);
  mpz_set_ui(td.pollardTestPrime, 2);
  mpz_set(td.mQrTest, td.mP);
  mpz_add_ui(td.mQrTest, td.mQrTest, 1);
  mpz_divexact_ui(td.mQrTest, td.mQrTest, 4);

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
      //mpz_set_ui(mPTest, a);
      //mpz_powm(mPTest, mPTest, mQrTest, td.mP);
      //mpz_powm_ui(mPTest, mPTest, 2, td.mP);
      if (a < 2000) {
        mpz_set_ui(primePow1, a);
        mpz_set(primePow2, primePow1);
        while (mpz_cmp_ui(primePow1, 2000) < 0) {
          mpz_swap(primePow2, primePow1);
          mpz_mul_ui(primePow1, primePow1, a);
        }
        mpz_mul(td.pollardTest, td.pollardTest, primePow2);
      }
      //if (a < 7000) {
        /*mpz_set_ui(primePow1, a);
        mpz_set(primePow2, primePow1);
        while (mpz_cmp_ui(primePow1, 3000) < 0) {
          mpz_swap(primePow2, primePow1);
          mpz_mul_ui(primePow1, primePow1, a);
        }*/
        //mpz_mul_ui(td.primeProd, td.primeProd, a);
      //}
      //if (mpz_cmp_d(mPTest, a) == 0) {
        insert(td.primes, a);

        /*mpz_set_ui(primePow1, a);
        mpz_set(primePow2, primePow1);
        while (mpz_cmp(primePow1, td.mP) < 0) {
          mpz_swap(primePow2, primePow1);
          mpz_pow_ui(primePow1, primePow1, 2);
        }
        mpz_mul_ui(td.pollardTest, td.pollardTest, a);
        printf("%d\n", a);*/
      //}
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
  mpz_clear(mPTest);
  clear(td.primes);
}

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
  mpz_t primeProd;
  mpz_t pollardTest;
  mpz_t pollardTestPrime;
  mpz_t mQrTest;
};

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
    mpz_powm_ui(y, y, 4, n);
    mpz_sub(diff, x, y);
    mpz_abs(diff, diff);

    mpz_gcd(d, diff, n);
  }
  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(diff);
  return i == limit ? -1 : 1;
}


void *find_smooth_function( void *ptr ) {
  int i;

  int found = 0, j;
  clock_t start, diff;

  int msec = 0;

  gmp_randstate_t rand_state;
  struct thread_data * td;
  td = (struct thread_data*)ptr;

  mpz_t quotient, remainder, num;

  gmp_randstate_t rand_state;
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


  while (1) {

    mpz_urandomm(pow, rand_state, td->mP);
    mpz_powm(testNum, base, pow, td->mP);
    int failed = 0, iterations = 0;
    found = 0;
    factors->size = 0;
    powers->size = 0;
    mpz_add_ui(pow, pow, 3);
    mpz_mod(pow, pow, td->mP);
    mpz_powm_ui(testNum, testNum, 3, td->mP);

    while (failed == 0 && mpz_sizeinbase(softTest, 10) > 50 && miller_rabin(softTest, 10, rand_state) != 1) {
      // Pollard's p-1
      mpz_powm(modResult, pollardTestPrime, td->pollardTest, softTest);
      mpz_sub_ui(modResult, modResult, 1);
      mpz_gcd(gcd, softTest, modResult);


      if (mpz_cmp_ui(gcd, largest) > 0) {

        failed = 1;
      }
      // We're not adjusting any parameters. If we find nothing, throw it away.
      if (mpz_cmp_ui(gcd, 1) == 0 || mpz_cmp(gcd, softTest) == 0) {
        failed = 1;
      }
      if (failed != 1) {
        mpz_divexact(softTest, softTest, gcd);
        iterations += 1;
      }
    }
    if (mpz_sizeinbase(softTest, 10) > 50) {
      attemptCount += 1;
      if (attemptCount % 100 == 0) {
        diff = clock() - start;
        msec = diff * 1000 / CLOCKS_PER_SEC;
        gmp_printf("P-1 factored to %Zd\n", softTest);
        printf("Factor failure size: %d Attempts: %lu, Time: %d, Iterations: %d\n", mpz_sizeinbase(softTest, 10), attemptCount, msec, iterations);
        start = clock();
      }
      continue;
    }
    failed = 0;

    iterations = 0;
    while (failed != 1 && mpz_sizeinbase(softTest, 10) > 24) {
      mpz_set_ui(x, 2);
      mpz_set_ui(y, 2);
      mpz_set_ui(d, 1);
      mpz_set_ui(gcd, 1);
      found = 0;
      for (i = 0; i < 6500 && mpz_cmp_ui(d, 1) == 0; ++i) {
        mpz_powm_ui(x, x, 2, softTest);
        mpz_add_ui(x, x, 1);
        mpz_powm_ui(y, y, 2, softTest);
        mpz_add_ui(y, y, 1);
        mpz_powm_ui(y, y, 2, softTest);
        mpz_add_ui(y, y, 1);
        mpz_sub(diffm, x, y);
        mpz_abs(diffm, diffm);
        mpz_gcd(d, diffm, softTest);
        if (mpz_cmp_ui(d, 1) != 0) {
          found = 1;
          mpz_divexact(softTest, softTest, d);
          break;
        }
      }
      if (found == 1) {
        iterations += 1;
      }
      if (found == 0 || mpz_cmp_ui(softTest, largest) > 0 && miller_rabin(softTest, 10, rand_state) == 1) {
        failed = 1;
      }
    }

    if (mpz_sizeinbase(softTest, 10) > 24 || failed == 1) {
      attemptCount += 1;
      if (attemptCount % 100 == 0) {
        diff = clock() - start;
        msec = diff * 1000 / CLOCKS_PER_SEC;
        gmp_printf("rho factored to %Zd\n", softTest);
        printf("Factor failure size: %d Attempts: %lu, Time: %d, Iterations: %d\n", mpz_sizeinbase(softTest, 10), attemptCount, msec, iterations);
        start = clock();
      }
      continue;
    }

    //gmp_printf("Factoring: %Zd.\n", softTest);
    //gmp_printf("Limit: %Zd.\n", limit);
    mpz_set(factorNum, softTest);
    //
    gmp_printf("Attempt: %d, Factoring: %Zd.\n", attemptCount, testNum);
    for (i = 1; i < numberOfPrimes; i++) {
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
      gmp_printf("Fully factored: %Zd\n", testNum);
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
  mpz_t mPTest, primePow1, primePow2, aBase;
  mpz_init(mPTest);
  mpz_init(aBase);
  mpz_init(primePow1);
  mpz_init(primePow2);

  mpz_init(td.mP);
  mpz_init(td.primeProd);
  mpz_init(td.mQrTest);
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
      if (a < 6500) {
        mpz_set_ui(primePow1, a);
        mpz_set(primePow2, primePow1);
        while (mpz_cmp_ui(primePow1, 6500) < 0) {
          mpz_swap(primePow2, primePow1);
          mpz_mul_ui(primePow1, primePow1, a);
        }
        mpz_mul(td.pollardTest, td.pollardTest, primePow2);
      }
      insert(td.primes, a);
    }
  }
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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>
#include <pthread.h>

#define MOD_ELEM_0 135291589527
#define MOD_ELEM_1 18446744073705357312U
#define MOD_ELEM_2 9223371487099224063
#define MOD_ELEM_3 512

struct thread_data {
  mpz_t mP;
  mpz_t mNr;
  mpz_t mQr;
  gmp_randstate_t rand_state;
};

void *find_smooth_function( void *ptr ) {
  struct thread_data * td;
  td = (struct thread_data*)ptr;
  mpz_t testNum, pow, phi;
  mpz_init(testNum);
  mpz_init(pow);
  mpz_init(phi);
  mpz_set(phi, td->mP);
  mpz_sub_ui(phi, phi, 1);
  mpz_divexact_ui(phi, phi, 2);
  mpz_urandomb(pow, td->rand_state, 200);
  mpz_set_ui(testNum, 5);
  mpz_powm(testNum, testNum, pow, td->mP);

  char filename[128];

  sprintf(filename, "./smooth_thread_%d", pthread_self());

  FILE *smoothFile;

  smoothFile = fopen(filename, "w+");

  while (1) {

    mpz_mul_ui(testNum, testNum, 5);

    mpz_mod(testNum, testNum, td->mP);

    mpz_add_ui(pow, pow, 1);

    mpz_mod(pow, pow, phi);

    mpz_t mTest, mGcd;
    mpz_init(mTest);
    mpz_init(mGcd);
    mpz_set(mTest, testNum);
    mpz_gcd(mGcd, mTest, td->mNr);


    if (mpz_cmp_ui(mGcd, 1) != 0) {
      //printf("Has a factor in mNr. Continuing.\n");
      mpz_clear(mTest);
      mpz_clear(mGcd);
      continue;
    }

    mpz_gcd(mGcd, mTest, td->mQr);
    int hasNonQFactor = 0;
    while (hasNonQFactor == 0 && mpz_cmp(mTest, mGcd) != 0) {
      if (mpz_cmp_ui(mGcd, 1) == 0) {
        //printf("Has a factor not in mQr. Continuing.\n");
        hasNonQFactor = 1;
        continue;
      }
      mpz_divexact(mTest, mTest, mGcd);
      mpz_gcd(mGcd, mTest, td->mQr);
    }
    mpz_clear(mTest);
    mpz_clear(mGcd);

    if (hasNonQFactor == 1) continue;

    fprintf(smoothFile, "%Zd\n", pow);

  }
  mpz_clear(testNum);
  mpz_clear(pow);
  mpz_clear(phi);
}

struct vec_t {
  uint32_t * elems;
  int size;
  int capacity;
};

void insert(struct vec_t* vec, uint32_t elem) {
  int i;
  if (vec->size == vec->capacity) {
    vec->capacity = vec->capacity * 2;
    uint32_t* newElems = (uint32_t *)malloc(vec->capacity* sizeof(uint32_t));
    for (i = 0; i < vec->size; ++i) {
      newElems[i] = vec->elems[i];
    }
    free(vec->elems);
    vec->elems = newElems;
  }
  vec->elems[vec->size] = elem;
  vec->size++;
}

struct vec_t* initVec() {
  struct vec_t* newVec = (struct vec_t*)malloc(sizeof(struct vec_t));
  newVec->capacity = 64;
  newVec->size = 0;
  newVec->elems = (uint32_t *)malloc(newVec->capacity * 2 * sizeof(uint32_t));
  return newVec;
}

int main(int argc, char ** argv) {

  int i;

  if (argc != 2) return 0;

  uint64_t* modulus = malloc(sizeof(uint64_t) * 4);
  modulus[0] = MOD_ELEM_0;
  modulus[1] = MOD_ELEM_1;
  modulus[2] = MOD_ELEM_2;
  modulus[3] = MOD_ELEM_3;

  printf("Loading prime factor base\n");

  FILE *fp;

  struct thread_data td;
  mpz_t mQrTest, mPTest;
  mpz_init(td.mQr);
  mpz_init(td.mNr);
  mpz_init(td.mP);
  mpz_init(mQrTest);
  mpz_init(mPTest);
  mpz_set_ui(td.mQr, 1);
  mpz_set_ui(td.mNr, 1);
  mpz_import(td.mP, 4, -1, sizeof(uint64_t), 0, 0, modulus);
  mpz_set(mQrTest, td.mP);
  mpz_add_ui(mQrTest, mQrTest, 1);
  mpz_divexact_ui(mQrTest, mQrTest, 4);

  fp = fopen("./primes.txt", "r");
  struct vec_t * primeVector = initVec();
  int c = 0;
  uint32_t a;
  while (c != -1) {
    a = 0;
    do {
      c = fgetc(fp);
      if (c >= '0' && c <= '9') {
        a *= 10;
        uint32_t n = c - '0';
        a += n;
      }
    } while (c >= '0' && c <= '9');
    if (a != 0) {
      mpz_set_ui(mPTest, a);
      mpz_powm(mPTest, mPTest, mQrTest, td.mP);
      mpz_powm_ui(mPTest, mPTest, 2, td.mP);
      if (mpz_cmp_d(mPTest, a) == 0) {
        mpz_mul_ui(td.mQr, td.mQr, a);
        insert(primeVector, a);
      } else {
        mpz_mul_ui(td.mNr, td.mNr, a);
      }
    }
  }

  mpz_clear(mQrTest);
  mpz_clear(mPTest);
  fclose(fp);

  printf("Done loading prime factor base\n");

  int numThreads = 0;

  for (i = 0; argv[1][i] != '\0'; ++i) {

    if (argv[1][i] >= '0' && argv[1][i] <= '9') {
      numThreads *= 10;
      int n = argv[1][i] - '0';
      numThreads += n;
    }
  }

  gmp_randinit_default(td.rand_state);

  printf("Starting with %d threads.\n", numThreads);

  pthread_t* threads = malloc(numThreads * sizeof(pthread_t));

  for (i = 0; i < numThreads; ++i) {
    pthread_create(&threads[0], NULL, find_smooth_function, &td);
  }

  for (i = 0; i < numThreads; ++i) {
    void* status;
    pthread_join(threads[i], &status);
  }

}

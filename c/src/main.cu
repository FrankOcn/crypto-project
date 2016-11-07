#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>
#include <cudam.cuh>

#define CHECK(call) \
{ \
const cudaError_t error = call; \
if (error != cudaSuccess) \
{ \
printf("Error: %s:%d, ", __FILE__, __LINE__); \
printf("code:%d, reason: %s\n", error, cudaGetErrorString(error)); \
exit(1); \
} \
}

#define MOD_ELEM_0 135291589527
#define MOD_ELEM_1 18446744073705357312
#define MOD_ELEM_2 9223371487099224063
#define MOD_ELEM_3 512

struct vec_t {
  uint32_t * elems;
  int size;
  int capacity;
};

void insert(struct vec_t* vec, uint32_t elem) {
  if (vec->size == vec->capacity) {
    vec->capacity = vec->capacity * 2;
    uint32_t* newElems = (uint32_t *)malloc(vec->capacity* sizeof(uint32_t));
    for (int i = 0; i < vec->size; ++i) {
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

__global__ void checkFactor(uint32_t* primes, uint32_t* results, cudam_t comp) {
  uint32_t prime = primes[(threadIdx.x) + blockIdx.x * blockDim.x];
  uint32_t* result = results + (threadIdx.x + blockIdx.x * blockDim.x);

  cudam_t c = cudam_init();
  cudam_t p = cudam_init();

  p[0] = prime;
  while (cudam_is_zero(c) == 1 && cudam_cmp(p, comp) != 1) {
    cudam_set(comp, c);
    cudam_mod(c, c, p);
    if (cudam_is_zero(c) == 1) {
      cudam_mult_ui(p, p, prime);
      *result = *result + 1;
    }
  }
  free(c);
  free(p);
}

int main(int argc, char ** argv) {

  cudam_t modulus = cudam_init();
  modulus[0] = MOD_ELEM_0;
  modulus[1] = MOD_ELEM_1;
  modulus[2] = MOD_ELEM_2;
  modulus[3] = MOD_ELEM_3;

  printf("Loading prime factor base\n");

  FILE *fp;

  mpz_t mQr, mNr, mP, mQrTest, mPTest;
  mpz_init(mQr);
  mpz_init(mNr);
  mpz_init(mP);
  mpz_init(mQrTest);
  mpz_init(mPTest);
  mpz_set_ui(mQr, 1);
  mpz_set_ui(mNr, 1);
  mpz_import(mP, 4, -1, sizeof(uint64_t), 0, 0, modulus);
  mpz_set(mQrTest, mP);
  mpz_add_ui(mQrTest, mQrTest, 1);
  mpz_divexact_ui(mQrTest, mQrTest, 4);

  fp = fopen("./primes.txt", "r");
  vec_t * primeVector = initVec();
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
      mpz_powm(mPTest, mPTest, mQrTest, mP);
      mpz_powm_ui(mPTest, mPTest, 2, mP);
      if (mpz_cmp_d(mPTest, a) == 0) {
        mpz_mul_ui(mQr, mQr, a);
        insert(primeVector, a);
      } else {
        mpz_mul_ui(mNr, mNr, a);
      }
    }
  }

  mpz_clear(mQrTest);
  mpz_clear(mPTest);
  fclose(fp);

  dim3 block (1024);
  dim3 grid (((primeVector->size)+block.x - 1)/block.x);

  printf("Done loading prime factor base\n");

  printf("Moving factor base into GPU\n");

  uint32_t * dPrimes;

  cudaMalloc((uint32_t**)&dPrimes, primeVector->size * sizeof(uint32_t));
  cudaMemcpy(dPrimes, primeVector->elems, primeVector->size * sizeof(uint32_t) , cudaMemcpyHostToDevice);

  printf("Creating results array\n");

  uint32_t* results = (uint32_t*) malloc(primeVector->size * sizeof(uint32_t));

  uint32_t * dResults;

  cudaMalloc((uint32_t**)&dResults, primeVector->size * sizeof(uint32_t));

  mpz_t testNum, pow;
  mpz_init(testNum);
  mpz_init(pow);
  mpz_set_ui(pow, 90);
  mpz_set_ui(testNum, 5);
  mpz_powm(testNum, testNum, pow, mP);

  uint64_t* testNumGPU;

  cudaMalloc((uint64_t**)&testNumGPU, 4* sizeof(uint64_t));

  while (1) {

    mpz_powm_ui(testNum, testNum, 2, mP);

    mpz_mul_ui(pow, pow, 2);
    mpz_mod(pow, pow, mP);

    gmp_printf("pow: %Zd\n", pow);

    mpz_t mTest, mGcd;
    mpz_init(mTest);
    mpz_init(mGcd);
    mpz_set(mTest, testNum);
    mpz_gcd(mGcd, mTest, mNr);


    if (mpz_cmp_ui(mGcd, 1) != 0) {
      //printf("Has a factor in mNr. Continuing.\n");
      mpz_clear(mTest);
      mpz_clear(mGcd);
      continue;
    }

    mpz_gcd(mGcd, mTest, mQr);
    int hasNonQFactor = 0;
    while (hasNonQFactor == 0 && mpz_cmp(mTest, mGcd) != 0) {
      if (mpz_cmp_ui(mGcd, 1) == 0) {
        //printf("Has a factor not in mQr. Continuing.\n");
        hasNonQFactor = 1;
        continue;
      }
      mpz_divexact(mTest, mTest, mGcd);
      mpz_gcd(mGcd, mTest, mQr);
    }

    mpz_clear(mTest);
    mpz_clear(mGcd);

    if (hasNonQFactor == 1) continue;


    gmp_printf("FOUND SMOOTH: %Zd\n", pow);

    return 0;

    /*// Reset results array
    for (int i = 0; i < primeVector->size; ++i) {
      results[i] = 0;
    }

    //printf("Move result array into GPU\n");

    cudaMemcpy(dResults, results, primeVector->size * sizeof(uint32_t), cudaMemcpyHostToDevice);

    //printf("testNum: %llu, %llu, %llu, %llu\n", testNum[0], testNum[1], testNum[2], testNum[3]);

    cudaMemcpy(testNumGPU, testNum, 4 * sizeof(uint64_t), cudaMemcpyHostToDevice);

    //printf("about to run\n");

    checkFactor <<< block, grid >>>(dPrimes, dResults, testNumGPU);
    CHECK(cudaDeviceSynchronize());
    cudaMemcpy(results, dResults, primeVector->size * sizeof(uint32_t), cudaMemcpyDeviceToHost);

    //printf("results copied back\n");

    cudam_zero(candidate);
    candidate[0] = 1;

    for (int i = 0; i < primeVector->size; ++i) {
      if (results[i] != 0) {
        printf("result: %d prime: %d\n", results[i], primeVector->elems[i]);
        for (int j = 0; j < results[i]; ++j) {
          cudam_mult_ui(candidate, candidate, primeVector->elems[i]);
        }
      }
    }

    for (int i = 0; i < 4; ++i) {
      printf("cand: %llu test: %llu\n", candidate[i], testNum[i]);
    }
    if (cudam_cmp(candidate, testNum) == 0) {
      printf("FOUND OMG FOUND ONE OMG OMG OMG %d\n", pow);

      for (int i = 0; i < 4; ++i) {
        printf("cand: %llu test: %llu\n", candidate[i], testNum[i]);
      }
      return 0;
      mpz_clear(mQr);
      mpz_clear(mNr);
      mpz_clear(mP);
    }*/
  }

}

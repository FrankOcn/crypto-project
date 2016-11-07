#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>

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

#define HIGH_BIT 0X8000000000000000

typedef uint64_t* uint256_t;

__host__ __device__ uint256_t initInt256() {
  uint256_t newInt = (uint256_t) malloc(sizeof(uint64_t) * 4);

  for (int i = 0; i < 4; ++i) {
    newInt[i] = 0;
  }
  return newInt;
}

// Returns -1 if a < b, 0 if a = b, and 1 if a > b
__host__ __device__ int cmpInt256(uint256_t a, uint256_t b) {
  for (int i = 3; i >= 0; --i) {
    if (a[i] < b[i]) {
      return -1;
    } else if (a[i] > b[i]) {
      return 1;
    }
  }
  return  0;
}

// Calculates c = a + b
__host__ __device__ void addInt256(uint256_t a, uint256_t b, uint256_t c) {
  int carry = 0;
  for (int i = 0; i < 4; ++i) {
    uint64_t* aCell = &(a[i]);
    uint64_t* bCell = &(b[i]);
    uint64_t* cCell = &(c[i]);
    if (*aCell > UINT64_MAX - *bCell) {
      *cCell = carry + *aCell - (UINT64_MAX - *bCell) - 1;
      carry = 1;
    } else if (*bCell > UINT64_MAX - *aCell) {
      *cCell = carry + *bCell - (UINT64_MAX - *aCell) - 1;
      carry = 1;
    } else {
      if (carry == 1 && (*aCell == UINT64_MAX || *bCell == UINT64_MAX)) {
        *cCell = 0;
      } else {
        *cCell = carry + *aCell + *bCell;
        carry = 0;
      }
    }
  }
  // carry will be 1 here if overflow, but our biggest number is 202 bits. This
  // isn't very likely.
}

// Calculates c = a + b
__host__ __device__ void addInt256U(uint256_t a, uint64_t b, uint256_t c) {
  int carry;
  uint64_t* aCell = &(a[0]);
  uint64_t* cCell = &(c[0]);
  if (*aCell > UINT64_MAX - b) {
    *cCell = *aCell - (UINT64_MAX - b) - 1;
    carry = 1;
  } else if (b > UINT64_MAX - *aCell) {
    *cCell = b - (UINT64_MAX - *aCell) - 1;
    carry = 1;
  } else {
    if (*aCell == UINT64_MAX || b == UINT64_MAX) {
      *cCell = 0;
      carry = 1;
    } else {
      *cCell = carry + *aCell + b;
      carry = 0;
    }
  }
  for (int i = 1; carry == 1 && i < 4; ++i) {
    uint64_t* cCell = &(c[i]);
    if (*cCell == UINT64_MAX) {
      *cCell = 0;
    } else {
      *cCell += carry;
      carry = 0;
    }
  }
  // carry will be 1 here if overflow, but our biggest number is 202 bits. This
  // isn't very likely.
}

// Shifts elements in place. Only works up to a shift of 64. Does not carry to
// a new cell
__host__ __device__ void sllInt256(uint256_t a, int shift) {
  if (shift <= 0) return;
  int cellShift = shift/64;
  int bitShift = shift%64;
  if (cellShift > 0) {
    for (int i = 3; i >= 0; --i) {
      if (i-cellShift < 0) {
        a[i] = 0;
      } else {
        a[i] = a[i-cellShift];
      }
    }
  }
  uint64_t high = UINT64_MAX << (64-bitShift);
  uint64_t carry = 0;
  for (int i = 0; i < 4; ++i) {
    uint64_t nextCarry = (a[i] & high) >> (64-shift);
    a[i] <<= shift;
    a[i] ^= carry;
    carry = nextCarry;
  }
}
// Copy a into b
__host__ __device__ void copyInt256(uint256_t a, uint256_t b) {
  for (int i = 0; i < 4; ++i) {
    b[i] = a[i];
  }
}

__host__ __device__ void clearInt256(uint256_t a) {
  for (int i = 0; i < 4; ++i) {
    a[i] = 0;
  }
}

// c = a * b
__host__ __device__ void multInt256U(uint256_t a, uint64_t b, uint256_t c) {
  if (b <= 0) return;
  uint256_t a_1 = initInt256();
  copyInt256(a, a_1);
  clearInt256(c);
  uint64_t bit = HIGH_BIT;
  while (bit != 0) {
    sllInt256(c, 1);
    if (b & bit) {
      addInt256(a_1, c, c);
    }
    bit >>= 1;
  }
  free(a_1);
}

__host__ uint256_t strToInt256(char * str) {
  uint256_t a = initInt256();
  uint256_t c = initInt256();
  for (int i = 0; str[i] != '\0'; ++i) {
    copyInt256(a, c);
    multInt256U(c, 10, a);
    uint64_t n = str[i] - '0';
    addInt256U(a, n, a);
  }
  return a;
}

// Calculates c = a - b
__host__ __device__ void subInt256(uint256_t a, uint256_t b, uint256_t c) {
  int borrow = 0;
  for (int i = 0; i < 4; ++i) {
    uint64_t* aCell = &(a[i]);
    uint64_t* bCell = &(b[i]);
    uint64_t* cCell = &(c[i]);
    if (*aCell >= *bCell + borrow) {
      *cCell = *aCell - *bCell - borrow;
      borrow = 0;
    } else {
      *cCell = UINT64_MAX - (*bCell + borrow - *aCell) + 1;
      borrow = 1;
    }
  }
}

// r = a%m
__host__ __device__ void modInt256(uint256_t a, uint256_t m, uint256_t r) {
  uint256_t a_1 = initInt256();
  copyInt256(a, a_1);
  clearInt256(r);
  uint64_t bit;
  int currentCell;
  for (int currentBit = 255; currentBit >= 0; --currentBit) {
    sllInt256(r, 1);
    if ((currentBit+1) % 64 == 0) {
      currentCell = (currentBit/64);
      bit = HIGH_BIT;
    }
    if (a_1[currentCell] & bit) {
      r[0] ^= 1;
    }
    if (cmpInt256(r, m) != -1) {
      subInt256(r, m, r);
    }
    bit >>= 1;
  }
  free(a_1);
}

__host__ __device__ uint32_t modInt256U(uint256_t a, uint32_t m) {
  uint32_t ovScale = ((UINT64_MAX % m) + 1) % m;
  uint32_t rem = 0;
  uint32_t overflow = 0;
  for (int i = 3; i >= 0; --i) {
    rem = a[i] % m;
    rem += (ovScale * overflow)%m;
    rem %=m;
    overflow = rem;
  }
  return rem;
}

__host__ __device__ void multInt256Mod(uint256_t a, uint256_t b, uint256_t c, uint256_t p) {
  int currentBit = 256;
  int currentCell;
  uint64_t bit;
  uint256_t a_1 = initInt256();
  uint256_t b_1 = initInt256();
  copyInt256(a, a_1);
  copyInt256(b, b_1);
  clearInt256(c);
  while (currentBit > 0) {
    if (currentBit % 64 == 0) {
      currentCell = (currentBit/64) - 1;
      bit = HIGH_BIT;
    }
    sllInt256(c, 1);
    if (b_1[currentCell] & bit) {
      addInt256(a_1, c, c);
    }
    modInt256(c, p, c);
    bit >>= 1;
    currentBit--;
  }
  free(a_1);
  free(b_1);
}

__host__ __device__ void multInt256(uint256_t a, uint256_t b, uint256_t c) {
  int currentBit = 256;
  int currentCell;
  uint64_t bit = HIGH_BIT;
  uint256_t a_1 = initInt256();
  copyInt256(a, a_1);
  clearInt256(c);
  while (currentBit > 0) {
    if (currentBit % 64 == 0) {
      currentCell = (currentBit/64) - 1;
      bit = HIGH_BIT;
    }
    sllInt256(c, 1);
    if (b[currentCell] & bit) {
      addInt256(a_1, c, c);
    }
    bit >>= 1;
    currentBit--;
  }
  free(a_1);
}

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

__host__ __device__ int isZeroInt256(uint256_t a) {
  if (a[0] == 0 && a[1] == 0 && a[2] == 0 && a[3] == 0) return 1;
  return -1;
}

__host__ __device__ int isOneInt256(uint256_t a) {
  if (a[0] == 1 && a[1] == 0 && a[2] == 0 && a[3] == 0) return 1;
  return -1;
}

__global__ void checkFactor(uint32_t* primes, uint32_t* results, uint256_t comp) {
  uint32_t prime = primes[(threadIdx.x) + blockIdx.x * blockDim.x];
  uint32_t* result = results + (threadIdx.x + blockIdx.x * blockDim.x);

  uint256_t c = initInt256();
  uint256_t p = initInt256();

  p[0] = prime;
  while (isZeroInt256(c) == 1 && cmpInt256(p, comp) != 1) {
    copyInt256(comp, c);
    modInt256(c, p, c);
    if (isZeroInt256(c) == 1) {
      multInt256U(p, prime, p);
      *result = *result + 1;
    }
  }
  free(c);
  free(p);
}

int main(int argc, char ** argv) {
  /*
  FILE *fp;

  fp = fopen("./primes.txt", "w+");
  struct vec_t* v = initVec();
  int composite;
  insert(v, 2);
  insert(v, 3);
  insert(v, 5);
  insert(v, 7);
  for (int i = 11; v->size < 1048576; i += 2) {
    composite = 0;
    for (int j = 0; j < v->size; ++j) {
      if (i % v->elems[j] == 0) {
        composite = 1;
        break;
      }
    }
    if (composite == 0) {
      insert(v, i);
      fprintf(fp, "%d\n", i);
      printf("%d\n", v->size);
    }
  }
  fclose(fp);
  return 0;
  */

  uint256_t modulus = initInt256();
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


  gmp_printf("mQr: %Zd\n", mQr);
  gmp_printf("mNr: %Zd\n", mNr);

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
  uint256_t five = initInt256();
  five[0] = 5;
  uint256_t testNum = initInt256();
  testNum[0] = 1;


  // Let's move it along a bit.
  int pow = 1;
  for (; pow < 87; pow++) {
    multInt256Mod(testNum, five, testNum, modulus);
    printf("power%d: %llu,%llu,%llu,%llu\n", pow, testNum[0], testNum[1], testNum[2], testNum[3]);
  }

  uint256_t max = initInt256();
  max[4] = 1;

  uint64_t* testNumGPU;

  cudaMalloc((uint64_t**)&testNumGPU, 4* sizeof(uint64_t));

  uint256_t candidate = initInt256();
  uint256_t fiveTest = initInt256();

  while (1) {

    multInt256Mod(testNum, testNum, testNum, modulus);

    pow += 1;

    printf("pow: %d\n", pow);

    modInt256(testNum, five, fiveTest);
    mpz_t mTest, mGcd;
    mpz_init(mTest);
    mpz_init(mGcd);
    mpz_import(mTest, 4, -1, sizeof(uint64_t), 0, 0, testNum);

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

    // Reset results array
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

    clearInt256(candidate);
    candidate[0] = 1;

    for (int i = 0; i < primeVector->size; ++i) {
      if (results[i] != 0) {
        printf("result: %d prime: %d\n", results[i], primeVector->elems[i]);
        for (int j = 0; j < results[i]; ++j) {
          multInt256U(candidate, primeVector->elems[i], candidate);
        }
      }
    }

    for (int i = 0; i < 4; ++i) {
      printf("cand: %llu test: %llu\n", candidate[i], testNum[i]);
    }
    if (cmpInt256(candidate, testNum) == 0) {
      printf("FOUND OMG FOUND ONE OMG OMG OMG %d\n", pow);

      for (int i = 0; i < 4; ++i) {
        printf("cand: %llu test: %llu\n", candidate[i], testNum[i]);
      }
      return 0;
      mpz_clear(mQr);
      mpz_clear(mNr);
      mpz_clear(mP);
    }
  }

}

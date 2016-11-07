#include <stdint.h>
#include <cudam.cuh>

#define HIGH_BIT 0X8000000000000000

typedef uint64_t* cudam_t;

// Allocates space for the cudam_t type. It's really just an array of
// four uint64_t elements
__host__ __device__ cudam_t cudam_init() {
  cudam_t newInt = (cudam_t) malloc(sizeof(uint64_t) * 4);

  for (int i = 0; i < 4; ++i) {
    newInt[i] = 0;
  }
  return newInt;
}

// Returns -1 if a < b, 0 if a = b, and 1 if a > b
__host__ __device__ int cudam_cmp(cudam_t a, cudam_t b) {
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
__host__ __device__ void cudam_add(cudam_t c, cudam_t a, cudam_t b) {
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
__host__ __device__ void cudam_add_ui(cudam_t c, cudam_t a, uint64_t b) {
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
__host__ __device__ void cudam_sll(cudam_t a, int shift) {
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
__host__ __device__ void cudam_set(cudam_t a, cudam_t b) {
  for (int i = 0; i < 4; ++i) {
    b[i] = a[i];
  }
}

__host__ __device__ void cudam_zero(cudam_t a) {
  for (int i = 0; i < 4; ++i) {
    a[i] = 0;
  }
}

// c = a * b
__host__ __device__ void cudam_mult_ui(cudam_t c, cudam_t a, uint64_t b) {
  if (b <= 0) return;
  cudam_t a_1 = cudam_init();
  cudam_set(a, a_1);
  cudam_zero(c);
  uint64_t bit = HIGH_BIT;
  while (bit != 0) {
    cudam_sll(c, 1);
    if (b & bit) {
      cudam_add(c, a_1, c);
    }
    bit >>= 1;
  }
  free(a_1);
}

__host__ cudam_t cudam_from_str(char * str) {
  cudam_t a = cudam_init();
  cudam_t c = cudam_init();
  for (int i = 0; str[i] != '\0'; ++i) {
    cudam_set(a, c);
    cudam_mult_ui(a, c, 10);
    uint64_t n = str[i] - '0';
    cudam_add_ui(a, a, n);
  }
  return a;
}

// Calculates c = a - b
__host__ __device__ void cudam_sub(cudam_t c, cudam_t a, cudam_t b) {
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
__host__ __device__ void cudam_mod(cudam_t r, cudam_t a, cudam_t m) {
  cudam_t a_1 = cudam_init();
  cudam_set(a, a_1);
  cudam_zero(r);
  uint64_t bit;
  int currentCell;
  for (int currentBit = 255; currentBit >= 0; --currentBit) {
    cudam_sll(r, 1);
    if ((currentBit+1) % 64 == 0) {
      currentCell = (currentBit/64);
      bit = HIGH_BIT;
    }
    if (a_1[currentCell] & bit) {
      r[0] ^= 1;
    }
    if (cudam_cmp(r, m) != -1) {
      cudam_sub(r, r, m);
    }
    bit >>= 1;
  }
  free(a_1);
}

__host__ __device__ void cudam_multm(cudam_t c,cudam_t a, cudam_t b, cudam_t p) {
  int currentBit = 256;
  int currentCell;
  uint64_t bit;
  cudam_t a_1 = cudam_init();
  cudam_t b_1 = cudam_init();
  cudam_set(a, a_1);
  cudam_set(b, b_1);
  cudam_zero(c);
  while (currentBit > 0) {
    if (currentBit % 64 == 0) {
      currentCell = (currentBit/64) - 1;
      bit = HIGH_BIT;
    }
    cudam_sll(c, 1);
    if (b_1[currentCell] & bit) {
      cudam_add(c, a_1, c);
    }
    cudam_mod(c, c, p);
    bit >>= 1;
    currentBit--;
  }
  free(a_1);
  free(b_1);
}

__host__ __device__ void cudam_mult(cudam_t c, cudam_t a, cudam_t b) {
  int currentBit = 256;
  int currentCell;
  uint64_t bit = HIGH_BIT;
  cudam_t a_1 = cudam_init();
  cudam_set(a, a_1);
  cudam_zero(c);
  while (currentBit > 0) {
    if (currentBit % 64 == 0) {
      currentCell = (currentBit/64) - 1;
      bit = HIGH_BIT;
    }
    cudam_sll(c, 1);
    if (b[currentCell] & bit) {
      cudam_add(c, a_1, c);
    }
    bit >>= 1;
    currentBit--;
  }
  free(a_1);
}

__host__ __device__ int cudam_is_zero(cudam_t a) {
  if (a[0] == 0 && a[1] == 0 && a[2] == 0 && a[3] == 0) return 1;
  return -1;
}

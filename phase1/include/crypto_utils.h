#ifndef __CRYPTO_UTILS_H__
#define __CRYPTO_UTILS_H__
#include "vec_mpz.h"

int miller_rabin(mpz_t, size_t, gmp_randstate_t rand_state);
int pollard_rho(mpz_t, mpz_t, size_t);

#endif

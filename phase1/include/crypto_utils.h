#ifndef __CRYPTO_UTILS_H__
#define __CRYPTO_UTILS_H__
#include "vec_mpz.h"

int miller_rabin(mpz_t, size_t, gmp_randstate_t rand_state, mpz_t* tmp);
int pollard_rho(mpz_t, mpz_t, size_t, mpz_t* tmp);
void eea_bounded_mpz( mpz_t beta, mpz_t p, mpz_t k, mpz_t l, mpz_t* tmp );

#endif

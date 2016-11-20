#include "crypto_utils.h"
#include "vec_mpz.h"


int miller_rabin(mpz_t n, size_t limit, gmp_randstate_t rand_state, mpz_t* tmp ) {
  int i, j, r;
  // There's a bug if n <= 4, so putting the answers here for good measure.
  if (mpz_cmp_ui(n, 1) == 0 || mpz_cmp_ui(n, 4) == 0) {
    return -1;
  }
  if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) {
    return 1;
  }
  mpz_set(tmp[0], n);
  mpz_set(tmp[2], n);
  mpz_set(tmp[3], n);
  mpz_sub_ui(tmp[0], tmp[0], 1);
  mpz_sub_ui(tmp[2], tmp[2], 1);
  mpz_sub_ui(tmp[3], tmp[3], 4);

  while (mpz_even_p(tmp[0]) != 0) {
    r += 1;
    mpz_tdiv_q_2exp(tmp[0], tmp[0], 1);
  }

  if (mpz_cmp_ui(tmp[0], 1) == 0) {
    return -1;
  }
  int continueOuter;
  for (i = 0; i < limit; ++i) {
    continueOuter = 0;
    mpz_urandomm(tmp[4], rand_state, tmp[3]);
    mpz_add_ui(tmp[4], tmp[4], 2);
    mpz_powm(tmp[1], tmp[4], tmp[0], n);
    if (mpz_cmp_ui(tmp[1], 1) == 0 || mpz_cmp(tmp[1], tmp[2]) == 0) {
      continue;
    }
    for (j = 0; j < r - 1; ++j) {
      mpz_powm_ui(tmp[1], tmp[1], 2, n);
      if (mpz_cmp_ui(tmp[1], 1) == 0) {
        return -1;
      }
      if (mpz_cmp(tmp[1], tmp[2]) == 1) {
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

int pollard_rho(mpz_t d, mpz_t n, size_t limit, mpz_t* tmp ) {
  int i;
  mpz_set_ui(tmp[0], 2);
  mpz_set_ui(tmp[1], 2);
  mpz_set_ui(d, 1);
  for (i = 0; i < limit && mpz_cmp_ui(d, 1) == 0; ++i) {
    mpz_powm_ui(tmp[0], tmp[0], 2, n);
    mpz_powm_ui(tmp[1], tmp[1], 4, n);
    mpz_sub(tmp[2], tmp[0], tmp[1]);
    mpz_abs(tmp[2], tmp[2]);

    mpz_gcd(d, tmp[2], n);
  }
  return i == limit ? -1 : 1;
}


// Find k and l such that k == l * beta (mod p) (where beta is really a^c)
// tmp == { yLast, x, y, quot, rem, m, n }
void eea_bounded_mpz( mpz_t beta, mpz_t p, mpz_t k, mpz_t l, mpz_t* tmp )
{
    // k = p;
    mpz_set( k, p );

    // tmp[1] = 1;
    // tmp[2] = 0;
    mpz_set_ui( tmp[1], 1 );
    mpz_set_ui( tmp[2], 0 );

    // l = 0;
    // tmp[0] = 1;
    mpz_set_ui( l, 0 );
    mpz_set_ui( tmp[0], 1 );

    // while ( mpz_cmp( beta, SQRT_P ) > 0 )
    while ( mpz_sizeinbase( beta, 2 ) > 101 )
    {
        // tmp[3] = k / beta;
        // tmp[4] = k % beta;
        mpz_fdiv_qr( tmp[3], tmp[4], k, beta );

        // tmp[5] = l - tmp[3] * tmp[1];
        mpz_swap( tmp[5], l );
        mpz_submul( tmp[5], tmp[3], tmp[1] );

        // tmp[6] = tmp[0] - tmp[3] * tmp[2];
        mpz_swap( tmp[6], tmp[0] );
        mpz_submul( tmp[6], tmp[3], tmp[2] );

        // l = tmp[1];
        // tmp[0] = tmp[2];
        mpz_swap( l, tmp[1] );
        mpz_swap( tmp[0], tmp[2] );

        // tmp[1] = tmp[5];
        // tmp[2] = n;
        mpz_swap( tmp[1], tmp[5] );
        mpz_swap( tmp[2], tmp[6] );

        // k = beta;
        // beta = tmp[4];
        mpz_swap( k, beta );
        mpz_swap( beta, tmp[4] );
    }
}

#include "crypto_utils.h"
#include "vec_mpz.h"


int miller_rabin(mpz_t n, size_t limit, gmp_randstate_t rand_state) {
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
  for (i = 0; i < limit; ++i) {
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

int pollard_rho(mpz_t d, mpz_t n, size_t limit) {
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

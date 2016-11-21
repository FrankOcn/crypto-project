#include "gmp.h"
#include "stdint.h"
// From the paper at http://cr.yp.to/papers/sf-20000807.pdf
// Given a positive integer b and an positive odd integer u, print a nonnegative
// integer v < 2^b such that 1 + uv /equiv 0 (mod 2^b)
//
// This is used as a subroutine in alg23
void alg21(mpz_t v, uint64_t b, mpz_t u) {
  if (b == 1) {
    mpz_set_ui(v, b);
    return;
  }
  // Set c = ceil(b/2)
  uint64_t c = b/2;
  c += b%2;

  mpz_t v_0, u_0, u_1, z, t_c, t_b, tmp;
  mpz_inits(v_0, u_0, u_1, z, t_c, t_b, tmp, 0);
  mpz_ui_pow_ui(t_c, 2, c);
  mpz_ui_pow_ui(t_b, 2, b);
  // Find v_0 < s^c such that 1 + uv_0 = 0 (mod 2^c)
  alg21(v_0, c, u);

  // Set u_0 = u (mod 2^c)
  mpz_mod(u_0, u, t_c);

  // Set u_1 = floor(u/2^c) (mod 2^c)
  mpz_set(u_1, u);
  mpz_fdiv_q(u_1, u_1, t_c);
  mpz_mod(u_1, u_1, t_b);

  // Set z = ((1 + u_0*v_0)/2^c + u_1*v_0)(mod 2^c)
  mpz_mul(z, u_0, v_0);
  mpz_add_ui(z, z, 1);
  mpz_divexact(z, z, t_c);
  mpz_mul(tmp, u_1, v_0);
  mpz_add(z, z, tmp);

  // Set v = v_0 + (2^c)*z*v_0 (mod 2^b)
  mpz_mul(v, z, v_0);
  mpz_mul(v, v, t_c);
  mpz_add(v, v, v_0);
  mpz_mod(v, v, t_b);

  mpz_clears(v_0, u_0, u_1, z, t_c, t_b, tmp, 0);
}

// From the paper at http://cr.yp.to/papers/sf-20000807.pdf
// Given a positive integer b, an odd positive integer u < 2^c, and a
// nonnegative integer x < 2^{c+b}, find a nonnegative integer r < 2^{c+1}
// such that 2^b*r = x (mod u)
//
// This is used as a subroutine in the remainder tree generation
void alg23(mpz_t r, uint64_t b, mpz_t x, mpz_t u) {
  mpz_t v, x_0, x_1, q, t_b;
  mpz_inits(v, x_0, x_1, q, t_b, 0);
  mpz_ui_pow_ui(t_b, 2, b);
  // Find v < 2^b such that 1 + uv = 0 (mod 2^b) using algo21
  alg21(v, b, u);
  // Set x_0 = x (mod 2^b) and x_1 = floor(x/2^b)
  mpz_fdiv_qr(x_1, x_0, x, t_b);
  // Set q = v*x_0 (mod 2^b)
  mpz_mul(q, v, x_0);
  mpz_mod(q, q, t_b);

  // Set r = x_1 + (x_0 + u*q)/2^b
  mpz_mul(r, q, u);
  mpz_add(r, r, x_0);
  mpz_divexact(r, r, t_b);
  mpz_add(r, r, x_1);
  mpz_clears(v, x_0, x_1, q, t_b, 0);

}

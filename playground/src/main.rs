extern crate gmp;

use gmp::mpz::Mpz;

fn main() {
    let p = Mpz::from_str_radix("1608507319599300885545233649993174934218196787014586228337099", 10).unwrap();
    let mut f: Mpz = 3u64.into();
    let two: Mpz = 2u64.into();
    while &p%&f != Mpz::zero() {
        f = f + &two;
    }
    println!("f: {}", f);
}

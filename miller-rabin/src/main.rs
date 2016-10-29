/// My first attempt at generating a strong prime.
///

extern crate gmp;
use gmp::mpz::Mpz;
use gmp::rand::RandState;
use std::time::{SystemTime, UNIX_EPOCH};

struct PrimeGenerator {
    rng: RandState
}

impl PrimeGenerator {
    fn new() -> PrimeGenerator {
        // Seed our RNG
        let t = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs();
        let mut rs = RandState::new();
        rs.seed_ui(t);
        PrimeGenerator { rng: rs }
    }
    fn mr(&mut self, n: &Mpz, k: u64) -> bool {

        // Determine d*2^r
        let mut d = n - 1;
        let mut r = 0;

        while &d & Mpz::one() == Mpz::zero() {
            d = d >> 1;
            r += 1;
        }

        'outer: for _ in 0..k {
            let a: Mpz = self.rng.urandom(&(n - 4)) + 2;
            let mut x = a.powm(&d, n);
            if x == Mpz::one() || x == n - 1 { continue 'outer; }

            for _ in 0..r {
                x = x.powm(&(2u64.into()), n);
                if x == Mpz::one() { return false; }
                if x == n - 1 { continue 'outer; }
            }
            return false;
        }
        true
    }
}

fn main() {
    let mut pg = PrimeGenerator::new();
    let p = Mpz::from_str_radix("1608507319599300885545233649993174934218196787014586228337099", 10).unwrap();
    if pg.mr(&p, 1000000000) {
        println!("{}", p);
    } else {
        println!("Not prime!");
    }
}

extern crate gmp;
use gmp::mpz::Mpz;
use std::env;

fn g(x: &Mpz, n: &Mpz) -> Mpz {
    (x.pow(2)-Mpz::one()) % n
}

fn rho(n: &Mpz) -> Mpz {
    let mut x: Mpz = Mpz::from_str_radix("2", 10).unwrap();
    let mut y: Mpz = Mpz::from_str_radix("2", 10).unwrap();
    let mut d = Mpz::one();
    while d == Mpz::one() {
        x = g(&x, n);
        y = g(&g(&y, n), n);
        d = (&x - &y).abs().gcd(n);
    }
    d
}

fn main() {

    let args = env::args().collect::<Vec<_>>();
    match args.len() {
        2 => {
            let d = Mpz::from_str_radix(&args[1], 10).unwrap();

            let f = rho(&d);

            println!("{}, {}", &f, d/&f);

        },
        _ => {
            println!("Usage: pollard-rho d")
        }
    }
}

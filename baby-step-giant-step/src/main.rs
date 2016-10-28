extern crate gmp;
extern crate rust_mpfr;

use rust_mpfr::mpfr::Mpfr;
use std::collections::HashMap;
use gmp::mpz::Mpz;
use std::env;

fn main() {

    let args = env::args().collect::<Vec<_>>();
    match args.len() {
        4 => {
            let pf = Mpfr::new2_from_str(300, args[1].as_str(), 10).unwrap();
            let m: u64 = (&pf.sqrt().floor()).into();

            let p = Mpz::from_str_radix(&args[1], 10).unwrap();
            let mut y = Mpz::from_str_radix(&args[3], 10).unwrap();
            let g = Mpz::from_str_radix(&args[2], 10).unwrap();
            let a_i = g.invert(&p).unwrap().powm(&m.into(), &p);
            let mut a_j = Mpz::one();
            let mut map = HashMap::new();
            println!("{}", m);
            for j in 0..m { map.insert(a_j.clone(), j); a_j = (a_j * &g).modulus(&p)}
            for i in 0..m {
                match map.get(&y) {
                    Some(j) => {
                        // If we have a match, return
                        let result = i * m + j;
                        println!("{}", result);
                        return;
                    },
                    None => {
                        y = (y * &a_i).modulus(&p);
                    }
                }
            }
        },
        _ => {
            println!("Usage: baby-giant <p> <g> <qr>")
        }
    }
}

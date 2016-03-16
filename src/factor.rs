extern crate num;
extern crate primal;

use std::cmp::min;
use std::collections::BTreeMap;

use self::num::integer::gcd;

use self::primal::is_prime;


const TRIAL_FACTOR_LIMIT: u32 = 100;


/// Factors an integer into its prime factors.
///
pub fn factor(n: u32) -> BTreeMap<u32, u32> {
    assert!(n != 0);

    let mut factorization = BTreeMap::new();

    // remove small factors
    let mut rest = n;
    let mut trial_factor = 2_u32;
    while trial_factor < TRIAL_FACTOR_LIMIT && rest >= trial_factor * trial_factor {
        let mut exponent = 0;
        while rest % trial_factor == 0 {
            rest /= trial_factor;
            exponent += 1;
        }
        if exponent != 0 {
            factorization.insert(trial_factor, exponent);
        }
        if trial_factor == 2 {
            trial_factor = 3;
        } else {
            trial_factor += 2;
        }
    }

    if rest < trial_factor * trial_factor {
        // rest is 1 or prime
        if rest > 1 {
            factorization.insert(rest, 1);
        }
    } else {
        // use Pollard's rho algorithm to find large factors
        let mut unfactored = vec![rest];
        while !unfactored.is_empty() {
            let u = unfactored.pop().unwrap();

            if is_prime(u as u64) {
                *factorization.entry(u).or_insert(0) += 1;
            } else {
                let f = find_large_factor(u);
                unfactored.push(f);
                unfactored.push(u / f);
            }
        }
    }
    factorization
}


/// Finds a factor using Pollard's rho algorithm.
///
/// The expected runtime is O(sqrt(f)) where f is the smallest factor of n.
/// This implies an expected runtime on the order of the fourth root of n.
/// It is the fastest algorithm currently known for up to about 100 bits.
///
/// n must be composite.
/// The returned factor may be composite.
///
fn find_large_factor(n: u32) -> u32 {

    /// Generates a pseudo-random sequence of numbers below m.
    ///
    /// c determines the sequence, should not be 0 or -2
    /// x is the previous value
    ///
    /// All values must fit in 32 bits to avoid overflow.
    ///
    fn next_random(m: u32, c: u32, x: u32) -> u32 {
        ((x as u64 * x as u64 + c as u64) % m as u64) as u32
    }

    let mut f: u32 = n;
    let mut c: u32 = 1;
    while f == n {
        let mut limit_power = 0;

        let mut k = 0;
        let mut y = 1;

        'search: loop {
            limit_power += 1;
            let limit = 1_u32 << limit_power;
            let product_limit = 1_u32 << (limit_power / 2);

            let x = y;
            while k < limit {
                let saved_y = y;

                // multiply together differences in the random sequence
                // doing gcd once for several numbers improves efficiency
                let mut product = 1;
                let mut j = 0;
                let j_limit = min(product_limit, limit - k);
                while j < j_limit {
                    y = next_random(n, c, y);
                    product = {
                        let n = n as u64;
                        let x = x as u64;
                        let y = y as u64;
                        let product = product as u64;
                        let x_minus_y = if x >= y {
                            x - y
                        } else {
                            n - y + x
                        };

                        ((product * x_minus_y) % n) as u32
                    };
                    j += 1;
                }
                f = gcd::<u32>(product, n);

                if f == 1 {
                    // no common factor found, move on
                    k += j_limit;
                } else {
                    // restart and find the factor
                    y = saved_y;
                    loop {
                        k == 1;
                        y = next_random(n, c, y);
                        let x_minus_y = if x >= y {
                            x - y
                        } else {
                            n - y + x
                        };
                        f = gcd::<u32>(x_minus_y, n);
                        if f != 1 {
                            break 'search;
                        }
                    }
                }
            }
        }

        // if we find only the trivial factor f == n
        // we retry with another random sequence
        c += 1;
    }
    f
}


/// Finds the largest prime factor of an integer.
///
pub fn largest_prime_factor(n: u32) -> Option<u32> {
    factor(n).keys().last().cloned()
}


#[cfg(test)]
mod tests {
    extern crate test;

    use std::collections::BTreeMap;
    use self::test::Bencher;

    use super::*;

    #[test]
    #[should_panic]
    fn test_zero() {
        factor(0);
    }

    #[test]
    fn test_one() {
        let expected = BTreeMap::new();
        let actual = factor(1);
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_two() {
        let mut expected = BTreeMap::new();
        expected.insert(2, 1);
        let actual = factor(2);
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_four() {
        let mut expected = BTreeMap::new();
        expected.insert(2, 2);
        let actual = factor(4);
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_six() {
        let mut expected = BTreeMap::new();
        expected.insert(2, 1);
        expected.insert(3, 1);
        let actual = factor(6);
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_power() {
        let mut expected = BTreeMap::new();
        expected.insert(2, 31);
        let actual = factor(1 << 31);
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_32bit_overflow() {
        let p1 = 65521_u32;
        let p2 = 65551_u32;

        let mut expected = BTreeMap::new();
        expected.insert(p1, 1);
        expected.insert(p2, 1);

        let actual = factor(p1 * p2);

        assert_eq!(expected, actual);
    }

    #[bench]
    fn bench_factor_low_range(b: &mut Bencher) {
        b.iter(|| factor_range(1000, 1000));
    }

    #[bench]
    fn bench_factor_high_range(b: &mut Bencher) {
        b.iter(|| factor_range(1_000_000_000, 1000));
    }

    fn factor_range(start: u32, size: u32) {
        for n in start..start + size {
            factor(n);
        }
    }

    #[bench]
    fn bench_factor_low_composite(b: &mut Bencher) {
        b.iter(|| factor(97 * 103));
    }

    #[bench]
    fn bench_factor_high_composite(b: &mut Bencher) {
        b.iter(|| factor(65521 * 65551));
    }

    #[bench]
    fn bench_factor_low_prime(b: &mut Bencher) {
        b.iter(|| factor(997));
    }

    #[bench]
    fn bench_factor_high_prime(b: &mut Bencher) {
        b.iter(|| factor(4294967291));
    }
}

use siphasher::sip::SipHasher13;
use std::hash::Hasher;

/// Derives 2r independent round keys from the 32-byte seed,
/// with domain separation per key index.
fn round_keys(seed: [u8; 32], rounds: usize) -> Vec<((u64, u64), (u64, u64))> {
    let master = SipHasher13::new_with_keys(
        u64::from_le_bytes(seed[0..8].try_into().unwrap()),
        u64::from_le_bytes(seed[8..16].try_into().unwrap()),
    );
    (0..(2 * rounds))
        .map(|i| {
            let mut h = SipHasher13::new_with_keys(
                u64::from_le_bytes(seed[16..24].try_into().unwrap()),
                u64::from_le_bytes(seed[24..32].try_into().unwrap()),
            );
            h.write_u64(master.finish() ^ 0x9E3779B97F4A7C15u64.wrapping_mul(i as u64 + 1));
            let a = h.finish();
            h.write_u8(0xFF); // domain-separate
            let b = h.finish();
            h.write_u8(0xFF);
            let c = h.finish();
            h.write_u8(0xFF);
            let d = h.finish();
            ((a, b), (c, d))
        })
        .collect()
}

/// Round function F: keyed SipHash of the half-word, mod p later.
fn round_f(x: u64, key: (u64, u64)) -> u64 {
    let mut h = SipHasher13::new_with_keys(key.0, key.1);
    h.write_u64(x);
    return h.finish();
}

pub fn permute<T, U>(n: T, x: T, s: U, r: usize) -> T 
where 
    T: Copy + PartialOrd 
    + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Mul<Output = T> 
    + std::ops::Rem<Output = T> + std::ops::Div<Output = T>
    + Into<u64> + TryFrom<u64>,
    <T as TryFrom<u64>>::Error: std::fmt::Debug,
    U: Into<[u8; 32]>,
{
    assert!(n.into() > 0, "n must be positive.");
    let keys = round_keys(s.into(), r);
    let p: T = T::try_from(n.into().isqrt()).unwrap();
    return keys.into_iter().fold(x, |acc, k| feistel_round_fwd(n, p, acc, k));
}

pub fn unpermute<T, U>(n: T, x: T, s: U, r: usize) -> T
where
    T: Copy + PartialOrd
    + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Mul<Output = T> 
    + std::ops::Rem<Output = T> + std::ops::Div<Output = T>
    + Into<u64> + TryFrom<u64>,
    <T as TryFrom<u64>>::Error: std::fmt::Debug,
    U: Into<[u8; 32]>,
{
    assert!(n.into() > 0, "n must be positive.");
    let keys = round_keys(s.into(), r);
    let p: T = T::try_from(n.into().isqrt()).unwrap();
    keys.into_iter().rev().fold(x, |acc, k| feistel_round_inv(n, p, acc, k))
}

fn feistel_round_fwd<T>(n: T, p: T, x: T, keys: ((u64, u64), (u64, u64))) -> T
where
    T: Copy + PartialOrd
    + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Mul<Output = T> 
    + std::ops::Rem<Output = T> + std::ops::Div<Output = T>
    + Into<u64> + TryFrom<u64>,
    <T as TryFrom<u64>>::Error: std::fmt::Debug,
{
    let p2 = p * p;
    if x < p2 {
        let x1: u64 = (x / p).into();
        let x2: u64 = (x % p).into();
        let y1 = (x1 + round_f(x2, keys.0)) % p.into();
        let y2 = (x2 + round_f(y1, keys.1)) % p.into();
        return p * T::try_from(y1).unwrap() + T::try_from(y2).unwrap() + n - p2;
    } 
    else {
        return x - p2;
    }
}

fn feistel_round_inv<T>(n: T, p: T, x: T, keys: ((u64, u64), (u64, u64))) -> T
where
    T: Copy + PartialOrd
    + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Mul<Output = T> 
    + std::ops::Rem<Output = T> + std::ops::Div<Output = T>
    + Into<u64> + TryFrom<u64>,
    <T as TryFrom<u64>>::Error: std::fmt::Debug,
{
    let p2 = p * p;
    if x < n - p2 {
        return x + p2;
    } else {
        // came from core: x = p*y2 + y1 + r  (where r = n - p2)
        let off: u64 = (x - (n - p2)).into();
        let p64: u64 = p.into();
        let y1 = off / p64;
        let y2 = off % p64;
        let x2 = (y2 + p64 - round_f(y1, keys.1) % p64) % p64;
        let x1 = (y1 + p64 - round_f(x2, keys.0) % p64) % p64;
        return T::try_from(p64 * x1 + x2).unwrap();
    }
}

#[test]
fn roundtrip() {
    let n: u32 = 1000; // p = 31, r = 1000 - 961 = 39
    for x in 0..n {
        assert_eq!(
            unpermute(n, permute(n, x, [7u8; 32], 8), [7u8; 32], 8),
            x
        );
    }
}
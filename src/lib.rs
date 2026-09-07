use rand::prelude::*;
use rand::rngs::StdRng;

pub fn permute<T, U>(n: T, x: T, s: U, r: usize) -> T 
where 
    T: Copy + PartialOrd 
    + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Mul<Output = T> 
    + std::ops::Rem<Output = T> + std::ops::Div<Output = T>
    + Into<u64> + TryFrom<u64>,
    <T as TryFrom<u64>>::Error: std::fmt::Debug,
    U: Into<[u8; 32]>,
{
    assert!(n.into()>0, "n must be positive.");
    let mut rng = StdRng::from_seed(s.into());
    let p: T = T::try_from(n.into().isqrt()).expect("Error converting back to T type.");
    let mut res = x;

    for _ in 0..r {
        let r1 = rng.next_u64();
        let r2 = rng.next_u64();
        res = feistel_round(n, p, res, r1, r2);
    }
    return res;
}

pub fn permute_batch<T, U>(n: T, x: &[T], s: U, r: usize) -> Vec<T>
where 
    T: Copy + PartialOrd 
    + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Mul<Output = T> 
    + std::ops::Rem<Output = T> + std::ops::Div<Output = T>
    + Into<u64> + TryFrom<u64>,
    <T as TryFrom<u64>>::Error: std::fmt::Debug,
    U: Into<[u8; 32]>,
{
    assert!(n.into()>0, "n must be positive.");
    let mut rng = StdRng::from_seed(s.into());
    let p: T = T::try_from(n.into().isqrt()).expect("Error converting back to T type.");
    let mut res: Vec<T> = x.to_vec();

    for _ in 0..r {
        let r1 = rng.next_u64();
        let r2 = rng.next_u64();
        res = res.iter().map(|xx| feistel_round(n, p, *xx, r1, r2)).collect();
    }
    return res;
}

fn feistel_round<T>(n: T, p: T, x: T, r1: u64, r2: u64) -> T
where 
    T: Copy + PartialOrd 
    + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Mul<Output = T> 
    + std::ops::Rem<Output = T> + std::ops::Div<Output = T>
    + Into<u64> + TryFrom<u64>,
    <T as TryFrom<u64>>::Error: std::fmt::Debug,
{
    // Flush the numbers from 0 to p^2-1
    // Then puts the remaining n-p^2 numbers under the stack
    let p2 = p*p;

    if x < p2 
    {
        let x1: T = x / p; 
        let x2: T = x % p;
        
        //Line and column flush
        let rr = StdRng::seed_from_u64(r1.wrapping_add(58353983u64.wrapping_mul(x2.into()))).next_u64();
        let y1: T = x2;
        let y2: T = T::try_from((x1.into() + rr) % p.into()).expect("Error converting back to T type.");

        let rr = StdRng::seed_from_u64(r2.wrapping_add(81960493u64.wrapping_mul(y2.into()))).next_u64();
    
        let z1: T = y2;
        let z2: T = T::try_from((y1.into() + rr) % p.into()).expect("Error converting back to T type.");
    
        return p*z1 + z2 + n - p2;
    }
    else {
        return x - p2;
    }
}
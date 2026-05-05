use std::ops::{Add, Mul, Neg, Sub};

mod butterfly;

/// Compile-time specified 31-bit modulus.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Modulus<const M: u32>;

impl<const M: u32> Modulus<M> {
    /// 31-bit modulus used for modular arithmetic
    const MODULUS: i64 = M as i64;

    /// Magic number for Plantard multiplication
    const MAGIC_D: i64 = {
        let m = M as u64;
        // 1 * 1 = 3 * 3 = 1 (mod 4)
        let mut inv_m = m & 3;
        // n inv_n = 1 (mod 2^k) => (n inv_n - 1)^2 = 0 (mod 2^{2k})
        // => n inv_n (2 - n inv_n) = 1 (mod 2^{2k})
        let mut i = u64::BITS.ilog2() - 1;
        while i > 0 {
            i -= 1;
            inv_m = inv_m.wrapping_mul(2_u64.wrapping_sub(m.wrapping_mul(inv_m)));
        }
        assert!(m.wrapping_mul(inv_m) == 1);

        inv_m as i64
    };

    /// Magic number for Plantard multiplication
    const MAGIC_A: i64 = {
        let lz = M.leading_zeros();
        assert!(lz < 31, "Modulus `M` should be at least 2.");

        1 << lz
            .checked_sub(1)
            .expect("Modulus `M` should be at most 31-bit integer.")
    };

    /// Product of operands of `imul` must be greater than or equal to this.
    pub const MIM_PRODUCT: i64 = {
        (1 - Self::MAGIC_A)
            .checked_mul(M as i64)
            .unwrap()
            .checked_shl(32)
            .unwrap()
    };

    /// Performs signed Plantard multiplication.
    ///
    /// # Constraints
    ///
    /// Let `A` to be `1 << (M.leading_zeros() - 1)`
    ///
    /// - `a * b >= -(2^A - 1) * M * 2^32`.
    /// - `a * b` must not overflow.
    #[inline(always)]
    pub const fn imul(a: i64, b: i64) -> i64 {
        let c = a.wrapping_mul(b).wrapping_mul(Self::MAGIC_D);
        ((c >> 32) + Self::MAGIC_A).wrapping_mul(Self::MODULUS) >> 32
    }


    /// Performs `a.pow(exp)`.
    ///
    /// # Constraints
    ///
    /// Let `A` to be `1 << (M.leading_zeros() - 1)`.
    ///
    /// - `a >= -(2^A - 1) * 2^32`.
    /// - `a^2` must not overflow.
    pub const fn pow(mut a: i64, mut exp: u32) -> i64 {
        let mut res = const { Self::i2p(1) };

        while exp > 0 {
            if exp & 1 == 1 {
                res = Self::imul(a, res);
            }

            exp >>= 1;
            a = Self::imul(a, a);
        }

        res
    }

    /// Converts an integer to its Plantard representation,
    /// i.e. computes `(-i * 2^64) mod M`.
    ///
    /// # Preconditions
    ///
    /// Let `A` to be `1 << (M.leading_zeros() - 1)`.
    ///
    /// - `i >= -(2^A - 1) * 2^32`.
    /// - `i * M` must not overflow.
    /// - Any value representable in `i32` satisfies these conditions.
    pub const fn i2p(i: i64) -> i64 {
        //`(2^128 % M) * MAGIC_D`
        let init: i64 = const {
            let pow_2_64 = (1 << 64) % M as u128;
            let pow_2_128 = pow_2_64 * pow_2_64 % M as u128;
            (pow_2_128 as i64).wrapping_mul(Self::MAGIC_D)
        };

        let c = i.wrapping_mul(init);
        ((c >> 32) + Self::MAGIC_A).wrapping_mul(Self::MODULUS) >> 32
    }

    /// Converts a value in Plantard representation to a standard integer,
    /// i.e. computes `(-p / 2^64) mod M`.
    ///
    /// # Preconditions
    ///
    /// Let `A` to be `1 << (M.leading_zeros() - 1)`.
    ///
    /// - `p >= -(2^A - 1) * M * 2^32`.
    /// - Any value representable in `i32` satisfies these conditions.
    pub const fn p2i(p: i64) -> i64 {
        let c = p.wrapping_mul(/* 1 times */ Self::MAGIC_D);
        ((c >> 32) + Self::MAGIC_A).wrapping_mul(Self::MODULUS) >> 32
    }

    pub const fn p2i2p(p: i64) -> i64 {
        let i = Self::p2i(p);
        Self::i2p(i)
    }

    /// Performs `p % self` in Plantard form.
    ///
    /// # Precondition
    ///
    /// Let `A` to be `1 << (M.leading_zeros() - 1)`.
    ///
    /// - `i >= -(2^A - 1) * 2^32`.
    /// - `i * M` must not overflow.
    /// - Any value representable in `i32` satisfies these conditions.
    pub const fn reduce(p: i64) -> i64 {
        // \bar{1}
        let one = const { Self::i2p(1).wrapping_mul(Self::MAGIC_D) };

        let c = p.wrapping_mul(one);
        ((c >> 32) + Self::MAGIC_A).wrapping_mul(Self::MODULUS) >> 32
    }
}

#[cfg(test)]
mod test_modulus {
    use super::Modulus;

    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn imul(a: i32, b: i32) {
            type M = Modulus::<998_244_353>;

            let prod = {
                let a = M::i2p(a as i64);
                let b = M::i2p(b as i64);

                let prod = M::imul(a, b);
                M::p2i(prod)
            };

            let naive = a as i64 * b as i64 % M::MODULUS;
            let test = if naive.is_negative() {
                [naive, naive + M::MODULUS]
            } else {
                [naive, naive - M::MODULUS]
            };

            assert!(test.contains(&prod), "{prod} will be in {test:?}")
        }
    }
}

pub trait NTTFriendlyPrime {
    /// A primitive root
    const PRIMITIVE_ROOT: i64;
    /// `P = A 2^L + 1 (A: odd)`
    const A: u32;
    /// `P = A 2^L + 1 (A: odd)`
    const L: u32;
}

impl NTTFriendlyPrime for Modulus<998_244_353> {
    const PRIMITIVE_ROOT: i64 = 3;
    const A: u32 = 998_244_353 >> Self::L;
    const L: u32 = (998_244_353 as u32 - 1).trailing_zeros();
}

#[derive(Debug, Clone)]
pub struct Polynomial<const M: u32>
where
    Modulus<M>: NTTFriendlyPrime,
{
    seq: Vec<i64>,
    /// seq[i].abs() < SC * M for all i
    scaling_factor: u32,
    /// Upper bound of degree.
    degree: usize,
}

impl<const M: u32> Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    pub fn new() -> Self {
        todo!()
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            seq: Vec::with_capacity(capacity),
            scaling_factor: 1,
            degree: 0,
        }
    }

    pub fn ones(n: usize) -> Self {
        todo!()
    }

    pub fn zeros(n: usize) -> Self {
        todo!()
    }

    /// Gets `i`-th coefficient.
    ///
    /// Consider to use [`eval`], [`sum`] or [`prod`] if possible for performance.
    ///
    /// [`eval`]: Self::eval
    /// [`sum`]: Self::sum
    /// [`prod`]: Self::prod
    pub fn get(&self, i: usize) -> Option<i32> {
        match self.seq.get(i) {
            Some(v) => Some(Modulus::<M>::p2i(*v) as i32),
            None => None,
        }
    }

    pub fn set(&mut self, i: usize, v: i32) -> bool {
        if let Some(u) = self.seq.get_mut(i) {
            *u = Modulus::<M>::i2p(v as i64);
            self.degree = self.degree.max(i);

            true
        } else {
            false
        }
    }

    pub fn eval(&self, x: i32) -> i32 {
        let max_sc_imul = const {
            Modulus::<M>::MIM_PRODUCT
                .unsigned_abs()
                .div_euclid(M as u64 * M as u64)
        };

        if (self.scaling_factor as u64) < max_sc_imul {
            let x = Modulus::<M>::i2p(x as i64);
            let result = self.seq[..=self.degree]
                .iter()
                .rev()
                .fold(0, |acc, v| Modulus::<M>::imul(acc, x) + v);
            Modulus::<M>::p2i(result) as i32
        } else {
            let result = self.seq[..=self.degree].iter().rev().fold(0, |acc, v| {
                // M < 2^31
                (acc * x as i64 + Modulus::<M>::p2i(*v)) % M as i64
            });
            result as i32
        }
    }

    pub fn sum(&self) -> i32 {
        let max_sc_p2i = const { Modulus::<M>::MIM_PRODUCT.unsigned_abs() / M as u64 };
        let chunk_size = {
            let max = max_sc_p2i / self.scaling_factor as u64;
            max.min(usize::MAX as u64) as usize
        };

        // |sum| < M 2^L < M^2
        let sum = self.seq[..=self.degree]
            .chunks(chunk_size)
            .fold(0, |sum, chunk| sum + Modulus::<M>::p2i(chunk.iter().sum()));
        (sum % M as i64) as i32
    }

    pub fn prod(&self) -> i32 {
        let max_sc_imul = const {
            Modulus::<M>::MIM_PRODUCT
                .unsigned_abs()
                .div_euclid(M as u64 * M as u64)
        };

        if self.scaling_factor as u64 <= max_sc_imul {
            let prod = self.seq[..=self.degree]
                .iter()
                .fold(const { Modulus::<M>::i2p(1) }, |acc, v| {
                    Modulus::<M>::imul(acc, *v)
                });
            Modulus::<M>::p2i(prod) as i32
        } else {
            let prod = self.seq[..=self.degree].iter().fold(1, |acc, v| {
                // |v| < INTERVAL * M
                acc * Modulus::<M>::p2i(*v) % M as i64
            });
            prod as i32
        }
    }
}

impl<const M: u32> Add for Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    type Output = Self;

    fn add(mut self, mut rhs: Self) -> Self::Output {
        // save allocation cost
        if self.seq.len() < rhs.seq.len() {
            std::mem::swap(&mut self, &mut rhs);
        }

        self.degree = self.degree.max(rhs.degree);

        if let Some(sum) = self.scaling_factor.checked_add(rhs.scaling_factor) {
            // no reduction
            Iterator::zip(self.seq.iter_mut(), rhs.seq).for_each(|(l, r)| *l += r);
            self.scaling_factor = sum;
        } else if M >> 29 == 0
        /* i.e. α >= 2 */
        {
            // add then reduce
            Iterator::zip(self.seq.iter_mut(), rhs.seq)
                .for_each(|(l, r)| *l = Modulus::<M>::p2i2p(*l + r));
            self.scaling_factor = 1;
        } else {
            // reduce
            let larger = std::cmp::max_by_key(&mut self, &mut rhs, |poly| poly.scaling_factor);
            Iterator::for_each(larger.seq.iter_mut(), |v| *v = Modulus::<M>::p2i2p(*v));
            larger.scaling_factor = 1;
            // then add
            Iterator::zip(self.seq.iter_mut(), rhs.seq).for_each(|(l, r)| *l += r);
            self.scaling_factor += rhs.scaling_factor;
        }

        self
    }
}

impl<const M: u32> Sub for Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl<const M: u32> Mul for Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let d = self.degree + rhs.degree;
        let n = (d + 1).next_power_of_two();

        let mut lhs = self.seq;
        lhs.resize(n, 0);

        let mut rhs = rhs.seq;
        rhs.resize(n, 0);

        Self::wrapping_mul_assign(&mut lhs, &mut rhs);

        Self {
            seq: lhs,
            scaling_factor: 1,
            degree: d,
        }
    }
}

impl<const M: u32> Neg for Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.seq.iter_mut().for_each(|v| *v = -*v);
        self
    }
}

impl<const M: u32> From<Vec<i32>> for Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    fn from(value: Vec<i32>) -> Self {
        let seq: Vec<_> = value
            .into_iter()
            .map(|i| Modulus::<M>::i2p(i as i64))
            .collect();
        let degree = seq.len().checked_sub(1).unwrap_or(0);

        Self {
            seq,
            scaling_factor: 1,
            degree,
        }
    }
}

impl<const M: u32> Extend<i32> for Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    fn extend<T: IntoIterator<Item = i32>>(&mut self, iter: T) {
        self.seq
            .extend(iter.into_iter().map(|i| Modulus::<M>::i2p(i as i64)));
        self.degree = self.seq.len().saturating_sub(1)
    }
}

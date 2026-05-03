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

    /// Performs signed Plantard multiplication.
    ///
    /// # Constraints
    ///
    /// TODO
    #[inline(always)]
    pub const fn imul(a: i64, b: i64) -> i64 {
        let c = a.wrapping_mul(b).wrapping_mul(Self::MAGIC_D);
        ((c >> 32) + Self::MAGIC_A).wrapping_mul(Self::MODULUS) >> 32
    }

    /// Performs `a.pow(exp)`.
    ///
    /// # Constraints
    ///
    /// `a * a` must not overflow.
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
    /// - `i >= -((2^A - 1) * M * 2^32)`, where `A = 1 << (M.leading_zeros() - 1)`.
    /// - `i < (2^64 - M * 2^(32 + A)) / (M - 1)`.
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
    /// - `p >= -(2^A - 1) * M * 2^32`, where `A = 1 << (M.leading_zeros() - 1)`.
    /// - Any value representable in `i32` satisfies this condition.
    pub const fn p2i(p: i64) -> i64 {
        let c = p.wrapping_mul(/* 1 times */ Self::MAGIC_D);
        ((c >> 32) + Self::MAGIC_A).wrapping_mul(Self::MODULUS) >> 32
    }

    /// Performs `p % self` in Plantard form.
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
    const PRIMITIVE_ROOT: i64;

    /// P = A 2^L + 1 (A: odd)
    const A: u32;
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
}

/// Low-level apis for multiplication of polynomials (convolution).
impl<const M: u32> Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    /// Performs an in-place Cooley–Tukey butterfly without normalization.
    ///
    /// This function may leave `seq` in a non-reduced state.
    ///
    /// Returns `true` if the resulting `seq` happens to be fully [`reduced`].
    ///
    /// [`reduced`]: Modulus::reduce
    ///
    /// # Preconditions
    ///
    /// - values in `seq` must be in Plantard form.
    /// - `seq` must be reduced.
    /// - `seq.len()` must be a power of two.
    /// - `seq.len() <= (1 << L)`, where `L = (M - 1).trailing_zeros()`.
    ///
    /// # Time complexity
    ///
    /// Θ(N log N), where N = `seq.len()`.
    pub fn butterfly(seq: &mut [i64]) -> bool {
        assert!(
            seq.len().is_power_of_two(),
            "`seq.len()` must be a power of two."
        );
        assert!(
            seq.len() >> Modulus::<M>::L <= 1,
            "Modulus `M` does not support NTT for this sequence length (too large)."
        );

        let lut: [i64; 32] = const {
            let mut lut = [0; _];

            let mut r = Modulus::<M>::pow(
                Modulus::<M>::i2p(Modulus::<M>::PRIMITIVE_ROOT),
                Modulus::<M>::A,
            );
            let mut ir = Modulus::<M>::pow(r, M.checked_sub(2).unwrap());

            let mut i = 2;
            let l = Modulus::<M>::L as usize;
            while i <= l {
                lut[l - i] = r;

                let mut j = l - i;
                while j + 2 < l {
                    j += 1;
                    lut[j] = Modulus::<M>::imul(lut[j], ir);
                }

                r = Modulus::<M>::imul(r, r);
                ir = Modulus::<M>::imul(ir, ir);
                i += 1;
            }

            lut
        };

        let mut w = seq.len() >> 1;
        let mut step = 0_u64;
        let interval = const {
            let lb = (1 - Modulus::<M>::MAGIC_A) /* times M */ << 32;
            lb.unsigned_abs() / M as u64
        };
        while w > 0 {
            // r = \bar{1}
            {
                let (pre, suf) = seq.split_at_mut(w);
                let (suf, _) = suf.split_at_mut(w);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = *s;
                    *s = *p - x;
                    *p = *p + x;
                });
            }

            let mut r = lut[0];
            for (i, pair) in seq.chunks_exact_mut(w << 1).enumerate().skip(1) {
                let (pre, suf) = pair.split_at_mut(w);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = Modulus::<M>::imul(*s, r);
                    *s = *p - x;
                    *p = *p + x;
                });

                r = Modulus::<M>::imul(r, lut[i.trailing_ones() as usize])
            }

            w >>= 1;

            step += 1;
            if step.is_multiple_of(interval) {
                seq.iter_mut().for_each(|v| *v = Modulus::<M>::reduce(*v));
            }
        }

        !step.is_multiple_of(interval)
    }

    /// Performs an in-place inverse Cooley–Tukey butterfly without normalization.
    ///
    /// This function may leave `seq` in a non-reduced state.
    ///
    /// Returns `true` if the resulting `seq` happens to be fully reduced.
    ///
    /// This is the inverse operation of [`butterfly`].
    ///
    /// # Preconditions
    ///
    /// - values in `seq` must be in Plantard form.
    /// - `seq` must be reduced.
    /// - `seq.len()` must be a power of two.
    /// - `seq.len() <= (1 << L)`, where `L = (M - 1).trailing_zeros()`.
    ///
    /// # Time complexity
    ///
    /// Θ(N log N), where N = `seq.len()`.
    pub fn butterfly_inv(seq: &mut [i64]) -> bool {
        assert!(
            seq.len().is_power_of_two(),
            "`seq.len()` must be a power of two."
        );
        assert!(
            seq.len() >> Modulus::<M>::L <= 1,
            "Modulus `M` does not support NTT for this sequence length (too large)."
        );

        let lut: [i64; 32] = const {
            let mut lut = [0; _];

            let mut ir = Modulus::<M>::pow(
                Modulus::<M>::i2p(Modulus::<M>::PRIMITIVE_ROOT),
                Modulus::<M>::A,
            );
            let mut r = Modulus::<M>::pow(ir, M.checked_sub(2).unwrap());

            let mut i = 2;
            let l = Modulus::<M>::L as usize;
            while i <= l {
                lut[l - i] = r;

                let mut j = l - i;
                while j + 2 < l {
                    j += 1;
                    lut[j] = Modulus::<M>::imul(lut[j], ir);
                }

                r = Modulus::<M>::imul(r, r);
                ir = Modulus::<M>::imul(ir, ir);
                i += 1;
            }

            lut
        };

        let mut w = 1;
        let mut step = 0_u64;
        let interval = const {
            let lb = (1 - Modulus::<M>::MAGIC_A) /* times M */ << 32;
            lb.unsigned_abs() / M as u64
        };
        while w < seq.len() {
            // r = \bar{1}
            {
                let (pre, suf) = seq.split_at_mut(w);
                let (suf, _) = suf.split_at_mut(w);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = *s;
                    *s = *p - x;
                    *p = *p + x;
                });
            }

            let mut r = lut[0];
            for (i, pair) in seq.chunks_exact_mut(w << 1).enumerate().skip(1) {
                let (pre, suf) = pair.split_at_mut(w);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = *s;
                    *s = Modulus::<M>::imul(*p - x, r);
                    *p = *p + x;
                });

                r = Modulus::<M>::imul(r, lut[i.trailing_ones() as usize])
            }

            w <<= 1;

            step += 1;
            if step.is_multiple_of(interval) {
                seq.iter_mut().for_each(|v| *v = Modulus::<M>::reduce(*v));
            }
        }

        !step.is_multiple_of(interval)
    }

    /// Performs an in-place convolution using Cooley–Tukey butterflies,
    /// storing the result in `lhs`.
    ///
    /// This function applies [`butterfly`] to `rhs` and then [`reduce`]s it,
    /// leaving it modified.
    ///
    /// [`reduce`]: Modulus::reduce
    ///
    /// # Preconditions
    ///
    /// - `lhs` and `rhs` must be in Plantard form.
    /// - `lhs` and `rhs` must be reduced.
    /// - `lhs.len() == rhs.len()`.
    /// - `lhs.len()` must be a power of two.
    ///
    /// # Time complexity
    ///
    /// Θ(N log N), where N = `seq.len()`.
    pub fn wrapping_mul_assign(lhs: &mut [i64], rhs: &mut [i64]) {
        assert_eq!(lhs.len(), rhs.len(), "lengths of operands must match");
        assert!(
            lhs.len().is_power_of_two(),
            "length of operands must be a power of two"
        );

        let frac_1_n = {
            // `1 / 2^i (mod M)`
            let lut = const {
                let mut lut = [0; 32];
                lut[0] = Modulus::<M>::i2p(1);
                lut[1] = Modulus::<M>::i2p((M + 1).div_ceil(2) as i64);

                let mut i = 2;
                while i < 32 {
                    lut[i] = Modulus::<M>::imul(lut[i - 1], lut[1]);
                    i += 1;
                }

                lut
            };

            let exp = lhs.len().trailing_zeros();
            lut[exp as usize]
        };

        Self::butterfly(lhs);
        Self::butterfly(rhs);
        lhs.iter_mut().zip(rhs.iter_mut()).for_each(|(l, r)| {
            // Since `rhs` is reduced, precondition of `imul` is never violated
            *r = Modulus::<M>::reduce(*r);
            *l = Modulus::<M>::imul(*l, *r);
        });
        Self::butterfly_inv(lhs);
        // normalize and reduce the result
        lhs.iter_mut()
            .for_each(|v| *v = Modulus::<M>::imul(*v, frac_1_n));
    }
}

#[test]
fn butterfly() {
    for n in (0..23).map(|d| 1 << d) {
        const MOD: u32 = 998_244_353;

        let mut seq = Vec::from_iter(0..n as i64);
        seq.extend_from_slice(&vec![0; n]);
        let test: Vec<_> = seq
            .iter()
            .map(|v| (v * n as i64 * 2).rem_euclid(MOD as i64))
            .collect();

        type P = Polynomial<MOD>;
        P::butterfly(&mut seq);
        P::butterfly_inv(&mut seq);
        seq.iter_mut().for_each(|v| *v = v.rem_euclid(MOD as i64));

        assert_eq!(seq, test)
    }
}

impl<const M: u32> Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    pub fn new(n: usize) -> Self {
        todo!()
    }

    pub fn ones(n: usize) -> Self {
        todo!()
    }

    pub fn zeros(n: usize) -> Self {
        todo!()
    }

    pub fn get(&mut self, i: usize) -> Option<i32> {
        match self.seq.get(i) {
            Some(v) => Some(Modulus::<M>::p2i(*v) as i32),
            None => None,
        }
    }

    pub fn set(&mut self, i: usize, v: i32) -> bool {
        if let Some(u) = self.seq.get_mut(i) {
            *u = Modulus::<M>::i2p(v as i64);
            true
        } else {
            false
        }
    }

    pub fn eval(&self, x: i32) -> i32 {
        let x = Modulus::<M>::i2p(x as i64);
        let result = self
            .seq
            .iter()
            .rev()
            // assume `v.abs() < 2^31`
            .fold(0, |acc, v| Modulus::<M>::imul(acc, x) + v);
        Modulus::<M>::p2i(result) as i32
    }

    pub fn sum(&self) -> i32 {
        // `seq.len() < M => sum.abs() < M 2^32`
        let sum: i64 = self.seq.iter().sum();
        Modulus::<M>::p2i(sum) as i32
    }

    pub fn prod(&self) -> i32 {
        let prod = self
            .seq
            .iter()
            .fold(1, |acc, v| Modulus::<M>::imul(acc, *v));
        todo!()
    }
}

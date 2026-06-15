pub struct Modulus<const M: u32>;

impl<const M: u32> Modulus<M> {
    const _CHECK: () = {
        assert!(M % 2 == 1, "Modulus `M` must be an odd integer");
        assert!(M >> 31 == 0, "Modulus `M` must be less than 2^31");
    };

    pub const INV_MODULUS: u64 = {
        let m = M as u64;
        // 1 * 1 = 3 * 3 = 1 (mod 4)
        let mut inv_m = m & 3;
        // Newton's method
        // n inv_n = 1 (mod 2^k) => (n inv_n - 1)^2 = 0 (mod 2^{2k})
        // => n inv_n (2 - n inv_n) = 1 (mod 2^{2k})
        let mut i = u64::BITS.ilog2() - 1;
        while i > 0 {
            i -= 1;
            inv_m = inv_m.wrapping_mul(2_u64.wrapping_sub(m.wrapping_mul(inv_m)));
        }
        assert!(m.wrapping_mul(inv_m) == 1);

        inv_m
    };

    /// Computes `a * b (mod M)` in Plantard representation.
    ///
    /// # Preconditions
    ///
    /// `a * b + M * 2^32 < 2^64` OR `a < M` OR `b < M`
    #[inline(always)]
    pub const fn mul(a: u32, b: u32) -> u32 {
        // FIXME: const hack
        debug_assert!(
            (a as u64 * b as u64)
                .checked_add((M as u64) << 32)
                .is_some(),
            "`a * b + M * 2^32` must fit in `u64`"
        );

        let x = (a as u64)
            .wrapping_mul(b as u64)
            .wrapping_mul(Self::INV_MODULUS)
            >> 32;
        let x = (x + 1) as u32 as u64 * M as u64;
        (x >> 32) as u32
    }

    /// Computes `p.pow(exp) (mod M)` in Plantard representation.
    ///
    /// # Time complexity
    ///
    /// O(log `exp`)
    pub const fn pow(mut p: u32, mut exp: u32) -> u32 {
        p %= M;

        let mut result = const { Self::i2p(1) };
        while exp > 0 {
            if exp & 1 == 1 {
                result = Self::mul(result, p);
            }
            // `p < M`
            p = Self::mul(p, p);
            exp >>= 1;
        }

        result
    }

    /// Converts an integer into Plantard representation.
    pub const fn i2p(i: u32) -> u32 {
        // `2^128 (mod M)`
        let pow_2_128 = const {
            let pow_2_64 = (M as u64).wrapping_neg() % M as u64;
            ((pow_2_64 * pow_2_64) % M as u64) as u32
        };

        Self::mul(i, pow_2_128)
    }

    /// Converts a value from Plantard representation back to a standard integer.
    pub const fn p2i(p: u32) -> u32 {
        Self::mul(p, 1)
    }
}

#[cfg(test)]
mod plantard_mul {
    use proptest::prelude::*;

    use super::Modulus;

    const P: u32 = 998_244_353;
    type M = Modulus<P>;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn conversion(i: u32) {
            let p = M::i2p(i);
            assert_eq!(M::p2i(p), i % P)
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 10))]
        #[test]
        fn pow(i: u32) {
            let p = M::i2p(i);
            let mut pow = M::i2p(1);
            for exp in 0..1 << 5 {
                assert_eq!(M::pow(p, exp), pow);
                pow = M::mul(pow, p);
            }
        }
    }
}

pub trait NTTFriendlyPrime {
    const PRIMITIVE_ROOT: u32;
    /// `P = (A << D) + 1`
    const A: u32;
    /// `P = (A << D) + 1`
    const D: u32;
}

impl NTTFriendlyPrime for Modulus<998_244_353> {
    const PRIMITIVE_ROOT: u32 = 3;
    const A: u32 = 998_244_353 >> Self::D;
    const D: u32 = (998_244_353_u32 - 1).trailing_zeros();
}

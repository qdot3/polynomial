#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mint<const M: u32>(u32);

impl<const M: u32> Mint<M> {
    const _CHECK_INVARIANT: () = {
        assert!(M % 2 == 1, "Modulus `M` must be an odd integer.");
        assert!(M > 2);
    };

    /// Modular inverse of `M` modulo `2^32`.
    const INV_M: u32 = {
        // 1 * 1 = 3 * 3 = 1 (mod 4)
        let mut inv_m = M & 3;
        // Newton's method
        // n inv_n = 1 (mod 2^k) => (n inv_n - 1)^2 = 0 (mod 2^{2k})
        // => n inv_n (2 - n inv_n) = 1 (mod 2^{2k})
        let mut i = u32::BITS.ilog2() - 1;
        while i > 0 {
            i -= 1;
            inv_m = inv_m.wrapping_mul(2_u32.wrapping_sub(M.wrapping_mul(inv_m)));
        }
        assert!(M.wrapping_mul(inv_m) == 1);

        inv_m
    };

    pub const fn new(x: u32) -> Self {
        let init = const {
            let pow_2_64 = (M as u64).wrapping_neg() % (M as u64);
            pow_2_64 as u32
        };

        // `x * init < (M << 32)`
        Self::mul(Self(x), Self(init))
    }

    /// Returns remainder.
    pub const fn get(self) -> u32 {
        self.mul(Self(1)).0
    }

    /// Performs modular multiplication.
    pub const fn mul(self, rhs: Self) -> Self {
        let (p_hi, p_lo) = {
            let prod = (self.0 as u64) * (rhs.0 as u64);
            debug_assert!(prod < (M as u64) << 32, "");

            ((prod >> 32) as u32, prod as u32)
        };

        let t_hi = {
            let t = p_lo.wrapping_mul(Self::INV_M) as u64 * (M as u64);
            (t >> 32) as u32
        };

        let (mut rem, borrow) = p_hi.overflowing_sub(t_hi);
        rem = rem.wrapping_add(if borrow { M } else { 0 });

        Self(rem)
    }

    pub const fn pow(mut self, mut exp: u32) -> Self {
        let mut res = Self::new(1);
        while exp > 0 {
            if exp & 1 == 1 {
                res = res.mul(self);
            }
            exp >>= 1;
            self = self.mul(self);
        }
        res
    }

    pub const fn inv(self) -> Option<Self> {
        if self.0 > 0 {
            Some(self.pow(M - 2))
        } else {
            None
        }
    }

    pub const fn add(self, rhs: Self) -> Self {
        self.sub(Self(M - rhs.0))
    }

    pub const fn sub(self, rhs: Self) -> Self {
        let (mut diff, borrow) = self.0.overflowing_sub(rhs.0);
        diff = diff.wrapping_add(if borrow { M } else { 0 });

        Self(diff)
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::Mint;

    const P: u32 = 998_244_353;
    type M = Mint<P>;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn mul(a in 1..P, b in 1..P) {
            let prod = M::new(a).mul(M::new(b)).get();
            let test = (a as u64 * b as u64 % P as u64) as u32;
            assert_eq!(prod, test);
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn inv(a in 1..P) {
            let v = M::new(a);
            let inv_v = v.inv().unwrap();

            assert_eq!(v.mul(inv_v).get(), 1);
        }
    }
}

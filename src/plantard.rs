#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mint<const M: u32>(u32);

impl<const M: u32> Mint<M> {
    const _CHECK: () = {
        assert!(M % 2 == 1, "Modulus `M` must be an odd integer");
        assert!(M >> 31 == 0, "Modulus `M` must be less than 2^31");
        assert!(M > 2);
    };

    const INV_M: u64 = {
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

    pub const fn new(x: u32) -> Self {
        // `2^128 (mod M)`
        let pow_2_128 = const {
            let pow_2_64 = (M as u64).wrapping_neg() % M as u64;
            ((pow_2_64 * pow_2_64) % M as u64) as u32
        };

        Self::mul(Self(x), Self(pow_2_128))
    }

    pub const fn get(self) -> u32 {
        self.mul(Self(1)).0
    }

    #[inline(always)]
    pub const fn mul(self, rhs: Self) -> Self {
        // FIXME: const hack
        debug_assert!(
            (self.0 as u64 * rhs.0 as u64)
                .checked_add((M as u64) << 32)
                .is_some(),
            "`prod + M * 2^32` must fit in `u64`"
        );

        let q = (self.0 as u64)
            .wrapping_mul(rhs.0 as u64)
            .wrapping_mul(Self::INV_M)
            >> 32;
        let rem = ((q + 1) as u32 as u64 * M as u64) >> 32;

        Self(rem as u32)
    }

    pub const fn pow(mut self, mut exp: u32) -> Self {
        let mut result = Self::new(1);
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(self);
            }
            // `p < M`
            self = self.mul(self);
            exp >>= 1;
        }

        result
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

use crate::{Modulus, NTTFriendlyPrime};

pub struct Butterfly<const M: u32>
where
    Modulus<M>: NTTFriendlyPrime;

impl<const M: u32> Butterfly<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    const _CHECK: () = {
        assert!(M >> 31 == 0, "Modulus `M` must be less than 2^31");
    };

    /// w <- w * LUT[i.trailing_ones()]
    const LUT: [u32; 30] = {
        let mut lut = [0; _];

        let mut r = Modulus::<M>::pow(
            Modulus::<M>::i2p(Modulus::<M>::PRIMITIVE_ROOT),
            Modulus::<M>::A,
        );
        let mut ir = Modulus::<M>::pow(r, M.checked_sub(2).unwrap());

        let mut i = 2;
        let l = Modulus::<M>::D as usize;
        while i <= l {
            lut[l - i] = r;

            let mut j = l - i;
            while j + 2 < l {
                j += 1;
                lut[j] = Modulus::<M>::mul(lut[j] as u64, ir as u64);
            }

            r = Modulus::<M>::mul(r as u64, r as u64);
            ir = Modulus::<M>::mul(ir as u64, ir as u64);
            i += 1;
        }

        lut
    };

    const LUT_INV: [u32; 30] = {
        let mut lut = Self::LUT;

        let mut i = 0;
        while i < lut.len() {
            // 0 -> 0
            lut[i] = Modulus::<M>::pow(lut[i], M.checked_sub(2).unwrap());
            i += 1;
        }

        lut
    };

    /// Performs in-place radix-2 Cooley–Tukey NTT.
    ///
    /// The input is in natural order, and the output is in bit-reversed order.
    ///
    /// If `REDUCE` is `true`, all output coefficients are reduced modulo `M`.
    ///
    /// # Preconditions
    ///
    /// - `seq.len().is_power_of_two()`
    /// - `seq.len() <= 1 << Modulus::<M>::D`
    /// - `seq[i] < M` for all `i`
    pub fn op<const REDUCE: bool>(seq: &mut [u32]) {
        assert!(seq.len().is_power_of_two());
        assert!(
            (seq.len() - 1) >> Modulus::<M>::D == 0,
            "`seq.len()` is too large"
        );
        debug_assert!(seq.iter().all(|v| *v < M));

        // maximum number of butterfly stages that can be accumulated
        // before a modular reduction is required.
        let max_scaling_factor = const {
            // (M-1) + n M < 2^32
            let n = (M.wrapping_neg() + 1) / M;
            // n >= 2 <=> M < (2^32 + 1)/3
            assert!(n > 1);
            n
        };

        let mut w = seq.len();
        let mut scaling_factor = 1;
        while w >= 2 {
            let mut r = const { Modulus::<M>::i2p(1) };
            for (i, pair) in seq.chunks_exact_mut(w).enumerate() {
                let (pre, suf) = pair.split_at_mut(w / 2);

                // every stage contributes at most one additional `M` to each coefficient.
                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = Modulus::<M>::mul(*s as u32 as u64, r as u64);
                    debug_assert!(x < M);

                    *s = *p + (M - x);
                    *p = *p + x;
                });

                // advance to the next twiddle factor.
                r = Modulus::<M>::mul(
                    r as u64,
                    Self::LUT[i.trailing_ones() as usize % Self::LUT.len()] as u64,
                )
            }
            w >>= 1;
            scaling_factor += 1;

            // reduce coefficients before the next doubling would overflow.
            if scaling_factor == max_scaling_factor {
                scaling_factor = 1;
                seq.iter_mut().for_each(|v| *v %= M);
            }
        }

        let reduced = scaling_factor == 1;
        if REDUCE && !reduced {
            seq.iter_mut().for_each(|v| *v %= M);
        }
    }

    /// Performs in-place radix-2 Cooley–Tukey INTT.
    ///
    /// The input is in bit-reversed order, and the output is in natural order.
    ///
    /// If `REDUCE` is `true`, all output coefficients are reduced modulo `M`.
    ///
    /// # Preconditions
    ///
    /// - `seq.len().is_power_of_two()`
    /// - `seq.len() <= 1 << Modulus::<M>::D`
    /// - `seq[i] < M` for all `i`
    pub fn inv<const REDUCE: bool>(seq: &mut [u32]) {
        assert!(seq.len().is_power_of_two());
        assert!(
            (seq.len() - 1) >> Modulus::<M>::D == 0,
            "`seq.len()` is too large"
        );
        debug_assert!(seq.iter().all(|v| *v < M));

        // upper bound on each coefficient to avoid overflow
        let upper_bound = const {
            assert!(M.leading_zeros() >= 1);
            M << M.leading_zeros()
        };

        let mut w = 2;
        let mut offset = M;
        while w <= seq.len() {
            let mut r = const { Modulus::<M>::i2p(1) };
            for (i, pair) in seq.chunks_exact_mut(w).enumerate() {
                let (pre, suf) = pair.split_at_mut(w / 2);

                // each stage can at most double a coefficient.
                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = *s;
                    *s = Modulus::<M>::mul((*p + offset - x) as u64, r as u64);
                    *p = *p + x;
                });

                // advance to the next twiddle factor.
                r = Modulus::<M>::mul(
                    r as u64,
                    Self::LUT_INV[i.trailing_ones() as usize % Self::LUT_INV.len()] as u64,
                )
            }
            w <<= 1;
            // `offset` now equals the current coefficient bound.
            offset <<= 1;

            // reduce coefficients before the next doubling would overflow.
            if offset == upper_bound {
                offset = M;
                seq.iter_mut().for_each(|v| *v %= M);
            }
        }

        let reduced = offset == M;
        if REDUCE && !reduced {
            seq.iter_mut().for_each(|v| *v %= M);
        }
    }

    /// Performs an in-place circular convolution.
    ///
    /// On return:
    ///
    /// - `lhs` contains the convolution result modulo `M`.
    /// - `rhs` contains the result of [`Butterfly::op::<true>()`](Self::op).
    ///
    /// # Preconditions
    ///
    /// - `lhs.len() == rhs.len()`
    /// - `lhs` and `rhs` satisfy the preconditions of [`Butterfly::op`]
    pub fn circular_convolution(lhs: &mut [u32], rhs: &mut [u32]) {
        assert_eq!(lhs.len(), rhs.len());

        let frac_1_n = {
            // LUT of 2^{-i} (mod M) in Plantard representation.
            let lut = const {
                let mut lut = [0; 32];

                let mut i = 0;
                let mut pow2 = Modulus::<M>::i2p(1);
                let two = Modulus::<M>::i2p(2) as u64;
                while i < lut.len() {
                    lut[i as usize] = Modulus::<M>::pow(pow2, M.checked_sub(2).unwrap());

                    pow2 = Modulus::<M>::mul(pow2 as u64, two);
                    i += 1;
                }

                lut
            };

            let exp = lhs.len().trailing_zeros();
            lut[exp as usize]
        };

        Self::op::<false>(lhs);
        Self::op::<true>(rhs);

        lhs.iter_mut().zip(rhs.iter()).for_each(|(l, r)| {
            // Since `r < M < 2^31`, precondition is always satisfied
            *l = Modulus::<M>::mul(*l as u64, *r as u64);
        });

        Self::inv::<false>(lhs);
        lhs.iter_mut().for_each(|l| {
            // Since `frac_1_n < M < 2^31`, precondition is always satisfied
            *l = Modulus::<M>::mul(*l as u64, frac_1_n as u64);
        });
    }
}

#[cfg(test)]
mod butterfly {
    use proptest::prelude::*;

    use super::{Butterfly, Modulus};

    const N: usize = 1 << 8;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 10))]
        #[test]
        fn butterfly(src in proptest::collection::vec(0..998_244_353_u32, N)) {
            let mut tar = src.clone();
            Butterfly::<998_244_353>::op::<true>(&mut tar);
            Butterfly::<998_244_353>::inv::<true>(&mut tar);

            assert!((0..N).all(|i| (src[i] as u64 * N as u64 % 998_244_353) as u32 == tar[i]));
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 10))]
        #[test]
        fn convolution(
            mut lhs in proptest::collection::vec(0..998_244_353_u32, N),
            mut rhs in proptest::collection::vec(0..998_244_353_u32, N),
        ) {
            lhs[N / 2..].fill(0);
            rhs[N / 2..].fill(0);

            type M = Modulus::<998_244_353>;

            lhs.iter_mut().for_each(|v| *v = M::i2p(*v as u32) );
            rhs.iter_mut().for_each(|v| *v = M::i2p(*v as u32) );

            let mut naive = vec![0; N];
            for i in 0..N / 2 {
                for j in 0.. N / 2 {
                    naive[i + j] = (
                        naive[i + j] + M::mul(lhs[i] as u64, rhs[j] as u64)
                    ) % 998_244_353;
                }
            }

            Butterfly::<998_244_353>::circular_convolution(&mut lhs, &mut rhs);

            assert!((0..N).all(|i| lhs[i] == naive[i]));
        }
    }
}

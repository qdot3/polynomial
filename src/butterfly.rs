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
    const LUT: [u32; 32] = {
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
                lut[j] = Modulus::<M>::mul(lut[j], ir);
            }

            r = Modulus::<M>::mul(r, r);
            ir = Modulus::<M>::mul(ir, ir);
            i += 1;
        }

        lut
    };

    const LUT_INV: [u32; 32] = {
        let mut lut = Self::LUT;

        let mut i = 0;
        while i < lut.len() {
            // 0 -> 0
            lut[i] = Modulus::<M>::pow(lut[i], M.checked_sub(2).unwrap());
            i += 1;
        }

        lut
    };

    /// Performs in-place radix-2 Cooley–Tukey NTT with reduction.
    ///
    /// The input is in natural order, and the output is in bit-reversed order.
    /// After completion, all elements satisfy `seq[i] < M`.
    ///
    /// # Preconditions
    ///
    /// - `seq.len().is_power_of_two()`
    /// - `seq.len() <= 1 << Modulus::<M>::D`
    /// - `seq[i] < M` for all `i`
    #[inline(always)]
    pub fn op(seq: &mut [u32]) {
        assert!(seq.len().is_power_of_two());
        assert!(
            (seq.len() - 1) >> Modulus::<M>::D == 0,
            "`seq.len()` is too large"
        );
        debug_assert!(seq.iter().all(|v| *v < M));

        let mut w = seq.len();
        while w >= 2 {
            let mut r = const { Modulus::<M>::i2p(1) };
            for (i, pair) in seq.chunks_exact_mut(w).enumerate() {
                let (pre, suf) = pair.split_at_mut(w / 2);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = Modulus::<M>::mul(*s, r);
                    debug_assert!(x < M);

                    // s < M
                    *s = p.wrapping_sub(x);
                    *s = s.wrapping_add(if s.cast_signed().is_negative() { M } else { 0 });
                    // p < M
                    *p = *p + x;
                    *p = p.wrapping_sub(if *p >= M { M } else { 0 });
                });

                // advance to the next twiddle factor.
                r = Modulus::<M>::mul(r, Self::LUT[i.trailing_ones() as usize % Self::LUT.len()])
            }
            w >>= 1;
        }
    }

    /// A target-feature-gated wrapper of [`Butterfly::op`].
    #[target_feature(enable = "avx2")]
    pub fn op_avx2(seq: &mut [u32]) {
        Self::op(seq);
    }

    /// Performs in-place radix-2 Cooley–Tukey INTT.
    ///
    /// The input is in bit-reversed order, and the output is in natural order.
    /// After completion, all elements satisfy `seq[i] < M`.
    ///
    /// # Preconditions
    ///
    /// - `seq.len().is_power_of_two()`
    /// - `seq.len() <= 1 << Modulus::<M>::D`
    /// - `seq[i] < M` for all `i`
    #[inline(always)]
    pub fn inv(seq: &mut [u32]) {
        assert!(seq.len().is_power_of_two());
        assert!(
            (seq.len() - 1) >> Modulus::<M>::D == 0,
            "`seq.len()` is too large"
        );
        debug_assert!(seq.iter().all(|v| *v < M));

        let mut w = 2;
        while w <= seq.len() {
            let mut r = const { Modulus::<M>::i2p(1) };
            for (i, pair) in seq.chunks_exact_mut(w).enumerate() {
                let (pre, suf) = pair.split_at_mut(w / 2);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = *s;
                    // s < M
                    *s = Modulus::<M>::mul(*p + M - x, r);
                    // p < M
                    *p = *p + x;
                    *p = p.wrapping_sub(if *p >= M { M } else { 0 });
                });

                // advance to the next twiddle factor.
                r = Modulus::<M>::mul(
                    r,
                    Self::LUT_INV[i.trailing_ones() as usize % Self::LUT_INV.len()],
                )
            }
            w <<= 1;
        }
    }

    /// A target-feature-gated wrapper of [`Butterfly::inv`].
    #[target_feature(enable = "avx2")]
    pub fn inv_axv2<const REDUCE: bool>(seq: &mut [u32]) {
        Self::inv(seq)
    }

    /// Performs in-place circular convolution.
    ///
    /// On return:
    ///
    /// - `lhs` contains the convolution result.
    /// - `rhs` contains the result of [`Butterfly::op(rhs)`](Self::op).
    ///
    /// # Preconditions
    ///
    /// - `lhs.len() == rhs.len()`
    /// - `lhs` and `rhs` satisfy the preconditions of [`Butterfly::op`]
    #[inline(always)]
    pub fn circular_convolution(lhs: &mut [u32], rhs: &mut [u32]) {
        assert_eq!(lhs.len(), rhs.len());

        let frac_1_n = {
            // LUT of 2^{-i} (mod M) in Plantard representation.
            let lut = const {
                let mut lut = [0; 32];

                let mut i = 0;
                let mut pow2 = Modulus::<M>::i2p(1);
                let two = Modulus::<M>::i2p(2);
                while i < lut.len() {
                    lut[i as usize] = Modulus::<M>::pow(pow2, M.checked_sub(2).unwrap());

                    pow2 = Modulus::<M>::mul(pow2, two);
                    i += 1;
                }

                lut
            };

            let exp = lhs.len().trailing_zeros();
            lut[exp as usize % lut.len()]
        };

        Self::op(lhs);
        Self::op(rhs);

        lhs.iter_mut().zip(rhs.iter()).for_each(|(l, r)| {
            // Since `l, r < M`, precondition is always satisfied
            *l = Modulus::<M>::mul(*l, *r);
        });

        Self::inv(lhs);
        lhs.iter_mut().for_each(|l| {
            // Since `l, frac_1_n < M`, precondition is always satisfied
            *l = Modulus::<M>::mul(*l, frac_1_n);
        });
    }

    /// A target-feature-gated wrapper of [`Butterfly::circular_convolution`].
    #[target_feature(enable = "avx2")]
    pub fn circular_convolution_avx2(lhs: &mut [u32], rhs: &mut [u32]) {
        Self::circular_convolution(lhs, rhs);
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
            Butterfly::<998_244_353>::op(&mut tar);
            Butterfly::<998_244_353>::inv(&mut tar);

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
                        naive[i + j] + M::mul(lhs[i] , rhs[j] )
                    ) % 998_244_353;
                }
            }

            Butterfly::<998_244_353>::circular_convolution(&mut lhs, &mut rhs);

            assert!((0..N).all(|i| lhs[i] == naive[i]));
        }
    }
}

use crate::{Mint, NTTFriendlyPrime, Prime};

pub struct Butterfly<const M: u32>
where
    Prime<M>: NTTFriendlyPrime;

impl<const M: u32> Butterfly<M>
where
    Prime<M>: NTTFriendlyPrime,
{
    const _CHECK: () = {
        assert!(M >> 31 == 0, "Modulus `M` must be less than 2^31");
        assert!(
            Self::LUT.len() > Prime::<M>::D as usize,
            "out-of-bounds error"
        );
        assert!(
            Self::LUT_INV.len() > Prime::<M>::D as usize,
            "out-of-bounds error"
        );
    };

    /// w <- w * LUT[i.trailing_ones()]
    const LUT: [Mint<M>; 32] = {
        let mut lut = [Mint::new(1); _];

        let mut r = Mint::new(Prime::PRIMITIVE_ROOT).pow(Prime::A);
        let mut ir = r.inv().expect("");

        let mut i = 2;
        let l = Prime::D as usize;
        while i <= l {
            lut[l - i] = r;

            let mut j = l - i;
            while j + 2 < l {
                j += 1;
                lut[j] = lut[j].mul(ir);
            }

            r = r.mul(r);
            ir = ir.mul(ir);
            i += 1;
        }

        lut
    };

    const LUT_INV: [Mint<M>; 32] = {
        let mut lut = Self::LUT;

        let mut i = 0;
        while i < lut.len() {
            // 0 -> 0
            lut[i] = lut[i].inv().unwrap();
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
    ///
    /// # Time complexity
    ///
    /// Θ(N log N)
    #[inline(always)]
    pub fn op(seq: &mut [Mint<M>]) {
        assert!(seq.len().is_power_of_two());
        let size = seq.len().trailing_zeros();
        assert!(size <= Prime::D);

        let mut w = 1_usize << size;
        while w >= 2 {
            let mut r = Mint::new(1);
            for (i, pair) in seq.chunks_exact_mut(w).enumerate() {
                let (pre, suf) = pair.split_at_mut(w / 2);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = s.mul(r);

                    *s = p.sub(x);
                    *p = p.add(x);
                });

                // advance to the next twiddle factor.
                r = r.mul(
                    // SAFETY: `i.trailing_ones() <= Modulus::D < LUT.len()`
                    *unsafe { Self::LUT.get_unchecked(i.trailing_ones() as usize) },
                );
            }
            w >>= 1;
        }
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
    ///
    /// # Time complexity
    ///
    /// Θ(N log N)
    #[inline(always)]
    pub fn inv(seq: &mut [Mint<M>]) {
        assert!(seq.len().is_power_of_two());
        assert!((seq.len() - 1) >> Prime::D == 0, "`seq.len()` is too large");

        let mut w = 2;
        while w <= seq.len() {
            let mut r = Mint::new(1);
            for (i, pair) in seq.chunks_exact_mut(w).enumerate() {
                let (pre, suf) = pair.split_at_mut(w / 2);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = *s;

                    *s = p.sub(x).mul(r);
                    *p = p.add(x);
                });

                // advance to the next twiddle factor.
                r = r.mul(
                    // SAFETY: `i.trailing_ones() <= Modulus::D < LUT_INV.len()`
                    *unsafe { Self::LUT_INV.get_unchecked(i.trailing_ones() as usize) },
                )
            }
            w <<= 1;
        }
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
    /// - `seq` must be represented in Plantard form
    /// - `lhs.len() == rhs.len()`
    /// - `lhs` and `rhs` satisfy the preconditions of [`Butterfly::op`].
    ///
    /// # Time complexity
    ///
    /// Θ(N log N)
    #[inline(always)]
    pub fn circular_convolution(lhs: &mut [Mint<M>], rhs: &mut [Mint<M>]) {
        assert_eq!(lhs.len(), rhs.len());

        let frac_1_n = {
            // LUT of 2^{-i} (mod M) in Plantard representation.
            let lut = const {
                let mut lut = [Mint::new(0); 32];

                let mut i = 0;
                while i < lut.len() {
                    lut[i as usize] = Mint::new(2).pow(i as u32).inv().expect("must not `0`");
                    i += 1;
                }

                lut
            };

            let exp = lhs.len().trailing_zeros();
            lut[exp as usize % lut.len()]
        };

        Self::op(lhs);
        Self::op(rhs);
        lhs.iter_mut()
            .zip(rhs.iter())
            .for_each(|(l, r)| *l = l.mul(*r).mul(frac_1_n));
        Self::inv(lhs);
    }

    /// A target-feature-gated wrapper of [`Butterfly::circular_convolution`].
    #[target_feature(enable = "avx2")]
    pub fn circular_convolution_avx2(lhs: &mut [Mint<M>], rhs: &mut [Mint<M>]) {
        Self::circular_convolution(lhs, rhs);
    }
}

#[cfg(test)]
mod butterfly {
    use proptest::prelude::*;

    use super::{Butterfly, Mint};

    const N: usize = 1 << 8;

    const P: u32 = 998_244_353;
    type M = Mint<P>;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 12))]
        #[test]
        fn butterfly(src in proptest::collection::vec(0..998_244_353_u32, N)) {
            let src: Vec<_> = src.into_iter().map(|v| M::new(v)).collect();
            let mut tar = src.clone();
            Butterfly::<998_244_353>::op(&mut tar);
            Butterfly::<998_244_353>::inv(&mut tar);

            assert!((0..N).all(|i| src[i].mul(M::new(N as u32)) == tar[i]));
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

            let mut lhs: Vec<_> = lhs.into_iter().map(|v| M::new(v)).collect();
            let mut rhs: Vec<_> = rhs.into_iter().map(|v| M::new(v)).collect();

            let mut naive = vec![M::new(0); N];
            for i in 0..N / 2 {
                for j in 0.. N / 2 {
                    naive[i + j] = lhs[i].mul(rhs[j]).add(naive[i + j])
                }
            }

            Butterfly::<998_244_353>::circular_convolution(&mut lhs, &mut rhs);

            assert!((0..N).all(|i| lhs[i] == naive[i]));
        }
    }
}

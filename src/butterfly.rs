use crate::{Modulus, NTTFriendlyPrime, Polynomial};

/// Low-level apis for multiplication of polynomials (convolution).
impl<const M: u32> Polynomial<M>
where
    Modulus<M>: NTTFriendlyPrime,
{
    /// w <- w * LUT[i.trailing_ones()]
    const LUT: [i64; 32] = {
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

    const LUT_INV: [i64; 32] = {
        let mut lut = Self::LUT;

        let mut i = 0;
        while i < lut.len() {
            // 0 -> 0
            lut[i] = Modulus::<M>::pow(lut[i], M.checked_sub(2).unwrap());
            i += 1;
        }

        lut
    };

    const SKIP_REDC_INTERVAL: u32 = {
        Modulus::<M>::MIM_PRODUCT
            .unsigned_abs()
            .div_euclid(M as u64 * M as u64)
            .checked_ilog2()
            .unwrap()
    };

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
        debug_assert!(
            seq.iter().all(|v| v.unsigned_abs() < M as u64),
            "`seq` must be reduced."
        );

        let mut w = seq.len() >> 1;
        let mut step = 0;
        while w > 0 {
            // r = \bar{1}
            {
                let (pre, suf) = seq[..w * 2].split_at_mut(w);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = *s;
                    *s = *p - x;
                    *p = *p + x;
                });
            }

            let mut r = Self::LUT[0];
            for (i, pair) in seq.chunks_exact_mut(w << 1).enumerate().skip(1) {
                let (pre, suf) = pair.split_at_mut(w);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = Modulus::<M>::imul(*s, r);
                    *s = *p - x;
                    *p = *p + x;
                });

                r = Modulus::<M>::imul(r, Self::LUT[i.trailing_ones() as usize])
            }

            w >>= 1;

            step += 1;
            if step == Self::SKIP_REDC_INTERVAL {
                step = 0;
                seq.iter_mut().for_each(|v| *v = Modulus::<M>::reduce(*v));
            }
        }

        step == 0
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
        debug_assert!(
            seq.iter().all(|v| v.unsigned_abs() < M as u64),
            "`seq` must be reduced."
        );

        let mut w = 1;
        let mut step = 0;
        while w < seq.len() {
            // r = \bar{1}
            {
                let (pre, suf) = seq[..2 * w].split_at_mut(w);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = *s;
                    *s = *p - x;
                    *p = *p + x;
                });
            }

            let mut r = Self::LUT_INV[0];
            for (i, pair) in seq.chunks_exact_mut(w << 1).enumerate().skip(1) {
                let (pre, suf) = pair.split_at_mut(w);

                pre.iter_mut().zip(suf.iter_mut()).for_each(|(p, s)| {
                    let x = *s;
                    *s = Modulus::<M>::imul(*p - x, r);
                    *p = *p + x;
                });

                r = Modulus::<M>::imul(r, Self::LUT_INV[i.trailing_ones() as usize])
            }

            w <<= 1;

            step += 1;
            if step == Self::SKIP_REDC_INTERVAL {
                step = 0;
                seq.iter_mut().for_each(|v| *v = Modulus::<M>::reduce(*v));
            }
        }

        step == 0
    }

    /// Performs an in-place convolution using Cooley–Tukey butterflies,
    /// storing the result in `lhs`.
    ///
    /// The result in `lhs` is normalized and reduced.
    /// This function applies [`butterfly`] to `rhs`, leaving it modified.
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
    /// Θ(N log N), where N = `lhs.len()`.
    pub fn wrapping_mul_assign(lhs: &mut [i64], rhs: &mut [i64]) {
        assert_eq!(lhs.len(), rhs.len(), "lengths of operands must match");
        assert!(
            lhs.len().is_power_of_two(),
            "length of operands must be a power of two"
        );
        debug_assert!(
            lhs.iter().all(|v| v.unsigned_abs() < M as u64),
            "`lhs` must be reduced."
        );
        debug_assert!(
            rhs.iter().all(|v| v.unsigned_abs() < M as u64),
            "`rhs` must be reduced."
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
            // Since scaling factors of `lhs` and `rhs` are the same,
            // at least one multiplication is available without reduction.
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

        let mut seq = Vec::from_iter((0..n as i64).map(|v| v % MOD as i64));
        seq.extend_from_slice(&vec![0; n]);
        let test: Vec<_> = seq
            .iter()
            .map(|v| (v * n as i64 * 2).rem_euclid(MOD as i64))
            .collect();

        type P = Polynomial<MOD>;
        if !P::butterfly(&mut seq) {
            // seq.iter_mut().for_each(|v| *v = Modulus::<MOD>::reduce(*v));
        }
        P::butterfly_inv(&mut seq);
        seq.iter_mut().for_each(|v| *v = v.rem_euclid(MOD as i64));

        assert_eq!(seq, test)
    }
}

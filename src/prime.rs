pub struct Prime<const M: u32>;

pub trait NTTFriendlyPrime {
    const PRIMITIVE_ROOT: u32;
    /// `P = (A << D) + 1`
    const A: u32;
    /// `P = (A << D) + 1`
    const D: u32;
}

macro_rules! ntt_prime_impl {
    ( $p:expr, $r:expr) => {
        impl NTTFriendlyPrime for Prime<$p> {
            const PRIMITIVE_ROOT: u32 = $r;

            const A: u32 = $p >> Self::D;
            const D: u32 = ($p - 1_u32).trailing_zeros();
        }
    };
}
ntt_prime_impl!(998_244_353, 3);



/// NTT-friendly prime numbers of the form `a * 2^22 + 1`, where `a` is odd.
pub const PRIME_22: [u32; 52] = [
    104857601, 113246209, 138412033, 155189249, 163577857, 230686721, 415236097, 666894337,
    683671553, 918552577, 935329793, 943718401, 985661441, 1161822209, 1212153857, 1321205761,
    1438646273, 1572864001, 1790967809, 1866465281, 2025848833, 2151677953, 2168455169, 2319450113,
    2344615937, 2369781761, 2403336193, 2470445057, 2629828609, 2671771649, 2680160257, 2705326081,
    2722103297, 2747269121, 2780823553, 2805989377, 2998927361, 3074424833, 3175088129, 3208642561,
    3435134977, 3451912193, 3510632449, 3552575489, 3577741313, 3602907137, 3628072961, 3686793217,
    3837788161, 3938451457, 4013948929, 4106223617,
];

/// NTT-friendly prime numbers of the form `a * 2^23 + 1`, where `a` is odd.
pub const PRIME_23: [u32; 20] = [
    377487361, 595591169, 645922817, 880803841, 897581057, 998244353, 1300234241, 1484783617,
    2088763393, 2558525441, 2810183681, 2910846977, 2994733057, 3112173569, 3313500161, 3414163457,
    3615490049, 3665821697, 3749707777, 4253024257,
];

/// NTT-friendly prime numbers of the form `a * 2^24 + 1`, where `a` is odd.
pub const PRIME_24: [u32; 9] = [
    754974721, 1224736769, 2130706433, 2533359617, 2634022913, 2868903937, 3238002689, 3942645761,
    4076863489,
];

/// NTT-friendly prime numbers of the form `a * 2^25 + 1`, where `a` is odd.
pub const PRIME_25: [u32; 6] = [
    167772161, 1107296257, 1711276033, 2113929217, 2717908993, 4194304001,
];

/// NTT-friendly prime numbers of the form `a * 2^26 + 1`, where `a` is odd.
pub const PRIME_26: [u32; 4] = [469762049, 1811939329, 2483027969, 2885681153];

/// NTT-friendly prime numbers of the form `a * 2^27 + 1`, where `a` is odd.
pub const PRIME_27: [u32; 3] = [2013265921, 2281701377, 3892314113];

/// NTT-friendly prime numbers of the form `a * 2^28 + 1`, where `a` is odd.
pub const PRIME_28: u32 = 3489660929;

/// NTT-friendly prime numbers of the form `a * 2^30 + 1`, where `a` is odd.
pub const PRIME_30: u32 = 3221225473;

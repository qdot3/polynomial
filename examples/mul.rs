/// Verified at <https://judge.yosupo.jp/problem/convolution_mod>
use std::io::{stdin, stdout, BufWriter, Write};

use output::IntBuffer;
use polynomial::{Butterfly, Modulus};
use reader::FastBufReader;

fn main() {
    let mut input = FastBufReader::<{ 1 << 18 }, _>::new(stdin().lock());
    let mut output = BufWriter::with_capacity(1 << 20, stdout().lock());
    let mut buf = IntBuffer::new();

    let n: usize = input.parse_next_token().unwrap();
    let m: usize = input.parse_next_token().unwrap();

    const MOD: u32 = 998_244_353;
    type M = Modulus<MOD>;

    let len = (n + m - 1).next_power_of_two();
    let mut a = {
        let mut a = Vec::with_capacity(len);
        for _ in 0..n {
            let value: u32 = input.parse_next_token().unwrap();
            a.push(M::i2p(value));
        }
        a.resize(len, 0);
        a
    };
    let mut b = {
        let mut b = Vec::with_capacity(len);
        for _ in 0..m {
            let value: u32 = input.parse_next_token().unwrap();
            b.push(M::i2p(value));
        }
        b.resize(len, 0);
        b
    };

    if is_x86_feature_detected!("avx2") {
        // SAFETY: guaranteed above
        unsafe { Butterfly::circular_convolution_avx2(&mut a, &mut b) };
    } else {
        Butterfly::circular_convolution(&mut a, &mut b);
    }

    for i in 0..n + m - 1 {
        let v = M::p2i(a[i]);

        let _ = output.write(buf.format(v).as_bytes());
        if i + 1 < n + m - 1 {
            let _ = output.write(b" ");
        }
    }
    output.write(b"\n").unwrap();
}

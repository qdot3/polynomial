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

    a.truncate(n + m - 1);
    a.iter_mut().for_each(|a| *a = M::p2i(*a));
    output.write(buf.format(a[0]).as_bytes()).unwrap();
    for a in a.into_iter().skip(1) {
        output.write(b" ").unwrap();
        output.write(buf.format(a).as_bytes()).unwrap();
    }
    output.write(b"\n").unwrap();
}

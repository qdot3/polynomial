/// Verified at <https://judge.yosupo.jp/problem/convolution_mod>
use std::io::{stdin, stdout, BufWriter, Write};

use output::IntBuffer;
use polynomial::Polynomial;
use reader::FastBufReader;

fn main() {
    let mut input = FastBufReader::<{ 1 << 18 }, _>::new(stdin().lock());

    let n: usize = input.parse_next_token().unwrap();
    let m: usize = input.parse_next_token().unwrap();
    let a: Vec<i32> = input.parse_next_token_vec(n).unwrap();
    let b: Vec<i32> = input.parse_next_token_vec(m).unwrap();

    let mut output = BufWriter::with_capacity(1 << 20, stdout().lock());
    let mut buf = IntBuffer::new();

    const MOD: u32 = 998_244_353;

    let a = Polynomial::<MOD>::from(a);
    let b = Polynomial::from(b);

    let c = a * b;
    for i in 0..n + m - 1 {
        let v = c.get(i).unwrap();
        let v = if v.is_negative() { v + MOD as i32 } else { v } as u32;

        let _ = output.write(buf.format(v).as_bytes());
        if i + 1 < n + m - 1 {
            let _ = output.write(b" ");
        }
    }
    output.write(b"\n").unwrap();
}

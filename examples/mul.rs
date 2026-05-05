/// Verified at <https://judge.yosupo.jp/problem/convolution_mod>
use std::io::stdin;

use input::{bind, FastInput};
use polynomial::Polynomial;

fn main() {
    let mut input: FastInput<std::io::StdinLock<'_>> = FastInput::new(stdin().lock());
    bind! { input >> n: usize, m: usize, a: [i32; n], b: [i32; m], }

    const MOD: u32 = 998_244_353;

    let a = Polynomial::<MOD>::from(a);
    let b = Polynomial::from(b);

    let c = a * b;

    let mut output = Vec::with_capacity(11 * (n + m));
    let mut buf = [0; 10];
    for i in 0..n + m - 1 {
        let i = c.get(i).unwrap();
        let mut i = if i.is_negative() { i + MOD as i32 } else { i } as u32;

        let mut j = buf.len();
        loop {
            let d = i % 10;
            i /= 10;

            j -= 1;
            buf[j] = d as u8 + b'0';

            if i == 0 {
                break;
            }
        }

        output.extend_from_slice(&buf[j..]);
        output.push(b' ');
    }
    output.pop();

    println!("{}", unsafe { String::from_utf8_unchecked(output) })
}

/// PCG64 (XSL-RR 128/64). Seeded via SplitMix64 expansion of a u64, giving a
/// reproducible permutation stream for the Monte-Carlo p-value. This is not
/// numpy's `Generator.permutation` bit-for-bit; it is a documented seeded
/// estimator (the pseudo-F statistic itself is exact).
pub struct Pcg64 {
    state: u128,
    inc: u128,
}

const MUL: u128 = 0x2360_ed05_1fc6_5da4_4385_df64_9fcc_f645;

impl Pcg64 {
    pub fn seed(seed: u64) -> Pcg64 {
        let mut sm = seed;
        let mut next = || {
            sm = sm.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = sm;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        let init_state = (u128::from(next()) << 64) | u128::from(next());
        let init_seq = (u128::from(next()) << 64) | u128::from(next());
        let mut rng = Pcg64 {
            state: 0,
            inc: (init_seq << 1) | 1,
        };
        rng.step();
        rng.state = rng.state.wrapping_add(init_state);
        rng.step();
        rng
    }

    #[inline]
    fn step(&mut self) {
        self.state = self.state.wrapping_mul(MUL).wrapping_add(self.inc);
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.step();
        let rot = (self.state >> 122) as u32;
        let xsl = ((self.state >> 64) ^ self.state) as u64;
        xsl.rotate_right(rot)
    }

    /// Uniform in `0..bound` via Lemire's multiply-shift with rejection.
    #[inline]
    pub fn below(&mut self, bound: u64) -> u64 {
        let mut x = self.next_u64();
        let mut m = u128::from(x) * u128::from(bound);
        let mut low = m as u64;
        if low < bound {
            let thresh = bound.wrapping_neg() % bound;
            while low < thresh {
                x = self.next_u64();
                m = u128::from(x) * u128::from(bound);
                low = m as u64;
            }
        }
        (m >> 64) as u64
    }

    /// In-place Fisher-Yates shuffle.
    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = self.below((i + 1) as u64) as usize;
            v.swap(i, j);
        }
    }
}

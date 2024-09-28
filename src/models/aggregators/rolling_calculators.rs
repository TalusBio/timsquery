// Rolling median calculator

pub struct RollingMedianCalculator<T: PartialOrd + Copy + Clone> {
    window_size: usize,
    data: Vec<(T, usize)>,
    index: usize,
}

impl<T: PartialOrd + Copy + Clone> RollingMedianCalculator<T> {
    pub fn new(window_size: usize) -> Self {
        Self {
            window_size,
            data: Vec::with_capacity(window_size),
            index: 0,
        }
    }

    pub fn add(&mut self, value: T) {
        self.data.push((value, self.index));
        self.index += 1;
        if self.data.len() > self.window_size {
            self.data.retain(|x| x.1 >= self.index - self.window_size);
        }
        self.reorder();
    }

    fn reorder(&mut self) {
        self.data
            .sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    }

    pub fn median(&self) -> Option<T> {
        if self.data.len() < self.window_size {
            None
        } else {
            let median_index = self.data.len() / 2;
            Some(self.data[median_index].0)
        }
    }
}

pub fn rolling_median<T: PartialOrd + Copy + Clone>(
    values: &[T],
    window_size: usize,
    pad_value: T,
) -> Vec<T> {
    let mut out = vec![pad_value; values.len()];
    let mut rolling = RollingMedianCalculator::new(window_size);
    let offset = window_size / 2;
    for (i, value) in values.iter().enumerate() {
        rolling.add(*value);
        if i >= (window_size - 1) {
            out[i - offset] = rolling.median().unwrap();
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rolling_median_calculator() {
        let mut calc = RollingMedianCalculator::new(3);
        calc.add(10.0);
        calc.add(20.0);
        calc.add(30.0);
        assert_eq!(calc.median(), Some(20.0));
        calc.add(1.0);
        calc.add(1.0);
        calc.add(1.0);
        calc.add(1.0);
        calc.add(1.0);
        calc.add(1.0);
        assert_eq!(calc.median(), Some(1.0));
    }

    #[test]
    fn test_rolling_median() {
        let input = vec![1.0, 2.0, 30.0, 4.0, 5.0, 60.0, 7.0, 8.0, 9.0];
        let out = rolling_median(&input, 3, f64::NAN);
        let expect_out = vec![f64::NAN, 2.0, 4.0, 5.0, 5.0, 7.0, 8.0, 8.0, f64::NAN];

        // assert_eq!(
        //     out,
        //     expect_out,
        // );
        //
        for i in 0..out.len() {
            if expect_out[i].is_nan() {
                assert!(out[i].is_nan());
            } else {
                assert!(
                    (out[i] - expect_out[i]).abs() < 1e-6,
                    "Expected {:?}, got {:?}",
                    expect_out,
                    out
                );
            }
        }
    }
}

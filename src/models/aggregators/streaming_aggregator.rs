use log::debug;

// Generic streaming aggregator that takes a pair of unsigned ints one with a value
// and another with a weight and in a streaming fashion adds the value to the accumulator
// to calculate the total, mean and variance.

#[derive(Debug, Clone, Copy)]
pub enum StreamingAggregatorError {
    DivisionByZero,
    NotEnoughData,
}

type Result<T> = std::result::Result<T, StreamingAggregatorError>;

pub struct StreamingMeanAggregator {
    pub aggregation: u128,
    pub weight: u128, // Pretty sure an u64 is enough ...
    pub points: u64,
}

impl StreamingMeanAggregator {
    pub fn new() -> Self {
        Self {
            aggregation: 0,
            weight: 0,
            points: 0,
        }
    }

    fn add(&mut self, value: u128, weight: u128) {
        self.aggregation += value;
        self.weight += weight;
        self.points += 1;
    }

    fn fold(&mut self, other: Self) {
        self.aggregation += other.aggregation;
        self.weight += other.weight;
    }

    fn mean(&self) -> Result<f64> {
        if self.weight == 0 {
            return Err(StreamingAggregatorError::NotEnoughData);
        }
        Ok(self.aggregation as f64 / self.weight as f64)
    }
}

/// Ref impl in javascript ...
/// https://nestedsoftware.com/2018/03/27/calculating-standard-deviation-on-streaming-data-253l.23919.html
/// https://nestedsoftware.com/2019/09/26/incremental-average-and-standard-deviation-with-sliding-window-470k.176143.html
/// Read the blog ... its amazing.
///
#[derive(Debug, Clone, Copy, Default)]
pub struct RunningStatsCalculator {
    weight: u64,
    mean_n: f64,
    d_: f64,
    // count: u64,
}

impl RunningStatsCalculator {
    pub fn new(weight: u64, mean: f64) -> Self {
        Self {
            weight,
            mean_n: mean,
            d_: 0.0,
            // count: 0,
        }
    }

    /// Add a new value to the running stats calculator.
    pub fn add(&mut self, value: f64, weight: u64) {
        // Update the mean
        // I am using a default weight of 1 for now ... not 100% sure how much
        // tha treally matters but it would fix the division by zero.
        debug_assert!(weight >= 1, "Weight must be >= 1");
        let weight_ratio = weight as f64 / self.weight as f64;
        let delta = value - self.mean_n;
        let last_mean_n = self.mean_n;
        self.mean_n += delta * weight_ratio;

        // Update the variance
        let to_add = weight as f64 * (value - self.mean_n) * (value - last_mean_n);
        self.d_ += to_add;

        // Update the weight
        self.weight += weight;
    }

    pub fn mean(&self) -> Result<f64> {
        if self.weight == 0 {
            return Err(StreamingAggregatorError::NotEnoughData);
        }
        Ok(self.mean_n)
    }

    pub fn variance(&self) -> Result<f64> {
        if self.weight == 0 {
            return Err(StreamingAggregatorError::NotEnoughData);
        }
        Ok(self.d_.abs() / self.weight as f64)
    }

    pub fn standard_deviation(&self) -> Result<f64> {
        let variance = self.variance()?;
        if variance.is_nan() {
            debug!("variance is nan, state -> {:?}", self);
        };
        if variance.is_infinite() {
            debug!("variance is infinite, state -> {:?}", self);
        };
        if variance.is_sign_negative() {
            debug!("variance is negative, state -> {:?}", self);
        };
        Ok(variance.sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_running_stats_calculator() {
        let mut calc = RunningStatsCalculator::new(10, 0.0);
        calc.add(10.0, 2);
        calc.add(10.0, 2);
        calc.add(10.0, 2);
        calc.add(10.0, 2);
        calc.add(10.0, 2);
        assert!(calc.mean().unwrap() < 5.6);
        assert!(calc.mean().unwrap() > 4.4);
        assert!(calc.variance().unwrap() > 15.);
        assert!(calc.variance().unwrap() < 25.);
        assert!(calc.standard_deviation().unwrap() > 4.5);
        assert!(calc.standard_deviation().unwrap() < 5.5);
    }

    // https://www.kaggle.com/datasets/carlmcbrideellis/data-anscombes-quartet?resource=download
    // ascombes quarted data from kaggle
    //
    // Both have real mean of 7.5 and std of 1.94
    const ASCOMBES_3: [f64; 11] = [
        7.46, 6.77, 12.74, 7.11, 7.81, 8.84, 6.08, 5.39, 8.15, 6.42, 5.73,
    ];
    const ASCOMBES_4: [f64; 11] = [
        6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.5, 5.56, 7.91, 6.89,
    ];

    #[test]
    fn test_running_stats_calculator_ascombes_3() {
        let mut calc = RunningStatsCalculator::new(1, ASCOMBES_3[0]);
        for i in 1..ASCOMBES_3.len() {
            calc.add(ASCOMBES_3[i], 1);
        }
        assert!(calc.mean().unwrap() < 7.6);
        assert!(calc.mean().unwrap() > 7.4);
        assert!(calc.standard_deviation().unwrap() > 1.92);
        assert!(calc.standard_deviation().unwrap() < 1.99);
    }

    #[test]
    fn test_running_stats_calculator_ascombes_4() {
        let mut calc = RunningStatsCalculator::new(1, ASCOMBES_4[0]);
        for i in 1..ASCOMBES_4.len() {
            calc.add(ASCOMBES_4[i], 1);
        }
        assert!(calc.mean().unwrap() < 7.6);
        assert!(calc.mean().unwrap() > 7.4);

        // Note that the tolerance here is a hair higher ... bc there
        // is an outlier value.
        assert!(
            calc.standard_deviation().unwrap() > 1.91,
            "Expected > 1.92, got {}",
            calc.standard_deviation().unwrap()
        );
        assert!(calc.standard_deviation().unwrap() < 1.99);
    }
}

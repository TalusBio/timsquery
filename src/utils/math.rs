/// Stirling's approximation for log factorial
// Code shamelessly stolen from
// https://github.com/lazear/sage/blob/445f5ab156c24f5c8f21098717692077e3b1d1ee/crates/sage/src/scoring.rs#L149C1-L157C2
//
pub fn lnfact(n: u16) -> f64 {
    if n == 0 {
        1.0
    } else {
        let n = n as f64;
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

pub fn lnfact_float(n: f64) -> f64 {
    if n < 1.0 {
        0.0
    } else {
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

/// Logarigthmic mean of a slice of values.
pub fn lnmean(vals: &[f64]) -> f64 {
    let mut sum = 0.0;
    for val in vals {
        sum += val.ln();
    }
    (sum / vals.len() as f64).exp()
}

fn cosine_similarity(a: &[f64], b: &[f64]) -> Option<f64> {
    // Check if vectors have the same length and are not empty
    if a.len() != b.len() || a.is_empty() {
        return None;
    }

    // Calculate dot product (numerator)
    let dot_product: f64 = a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum();

    // Calculate magnitudes (denominator)
    let magnitude_a: f64 = a.iter().map(|&x| x * x).sum::<f64>().sqrt();

    let magnitude_b: f64 = b.iter().map(|&x| x * x).sum::<f64>().sqrt();

    // Avoid division by zero
    if magnitude_a == 0.0 || magnitude_b == 0.0 {
        return None;
    }

    // Calculate cosine similarity
    Some(dot_product / (magnitude_a * magnitude_b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cosine_similarity() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let result = cosine_similarity(&a, &b).unwrap();
        assert!((result - 0.974631846).abs() < 1e-8);
    }

    #[test]
    fn test_identical_vectors() {
        let a = vec![1.0, 1.0, 1.0];
        let result = cosine_similarity(&a, &a).unwrap();
        assert!((result - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_empty_vectors() {
        let a: Vec<f64> = vec![];
        let b: Vec<f64> = vec![];
        assert!(cosine_similarity(&a, &b).is_none());
    }

    #[test]
    fn test_different_lengths() {
        let a = vec![1.0, 2.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!(cosine_similarity(&a, &b).is_none());
    }

    #[test]
    fn test_zero_vector() {
        let a = vec![0.0, 0.0, 0.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!(cosine_similarity(&a, &b).is_none());
    }
}

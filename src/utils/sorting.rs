fn place_at_indices<T>(original: &mut [T], indices: &mut [usize]) {
    for i in 0..indices.len() {
        while i != indices[i] {
            let new_i = indices[i];
            indices.swap(i, new_i);
            original.swap(i, new_i);
        }
    }
}

fn place_at_indices_n<T>(original: &mut [&mut [T]], indices: &mut [usize]) {
    for i in 0..indices.len() {
        while i != indices[i] {
            let new_i = indices[i];
            indices.swap(i, new_i);
            for slice in original.iter_mut() {
                slice.swap(i, new_i);
            }
        }
    }
}

pub fn argsort_by<T, F, K>(v: &[T], key: F) -> Vec<usize>
where
    F: Fn(&T) -> K,
    K: Ord,
{
    let mut indices: Vec<usize> = (0..v.len()).collect();
    indices.sort_by_key(|&i| key(&v[i]));
    indices
}

/// Sorts all passed slices by the key function.
///
/// For a similar macro that works in slices of heterogeneous types, see `sort_by_indices_multi!`.
///
/// # Example
/// ```
/// use timsquery::utils::sorting::sort_multiple_by;
/// let mut va = vec![9, 8, 7];
/// let mut vb = vec![1, 2, 3];
/// sort_multiple_by(&mut [&mut va, &mut vb], |x| *x);
/// assert_eq!(va, vec![7, 8, 9]);
/// assert_eq!(vb, vec![3, 2, 1]);
/// ```
///
pub fn sort_multiple_by<T, F, K>(vectors: &mut [&mut [T]], key: F)
where
    F: Fn(&T) -> K,
    K: Ord,
{
    let mut indices = argsort_by(vectors[0], |x| key(x));
    place_at_indices_n(vectors, &mut indices);
}

/// Macro that sorts an arbitrary number of slices by a key function applied to the
/// first slice.
///
/// # Example
/// ```
/// use timsquery::sort_by_indices_multi;
/// use timsquery::utils::sorting::argsort_by;
///
/// let mut va = vec![9, 8, 7];
/// let mut vb = vec![1, 2, 3];
/// let mut vc = vec!['a', 'b', 'c'];
/// let mut indices = argsort_by(&va, |x| *x);
/// sort_by_indices_multi!( indices, &mut va, &mut vb, &mut vc);
///
/// assert_eq!(va, vec![7, 8, 9]);
/// assert_eq!(vb, vec![3, 2, 1]);
/// assert_eq!(vc, vec!['c', 'b', 'a']);
/// ```
///
///
#[macro_export]
macro_rules! sort_by_indices_multi {
    ($indices:expr, $($data:expr),+ $(,)?) => {{
        let mut indices = $indices.to_vec();
        for idx in 0..indices.len() {
            if indices[idx] != idx {
                let mut current_idx = idx;
                loop {
                    let target_idx = indices[current_idx];
                    indices[current_idx] = current_idx;
                    if indices[target_idx] == target_idx {
                        break;
                    }
                    $(
                        $data.swap(current_idx, target_idx);
                    )+
                    current_idx = target_idx;
                }
            }
        }
    }};
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_place_at_indices() {
        let mut va = vec![1, 2, 3, 4];
        let mut indices = vec![3, 2, 1, 0];
        place_at_indices(&mut va, &mut indices);
        assert_eq!(va, vec![4, 3, 2, 1]);
    }

    #[test]
    fn test_place_at_indices_n() {
        let mut va = vec![1, 2, 3, 4];
        let mut vb = vec![5, 6, 7, 8];
        let mut indices = vec![3, 2, 1, 0];
        place_at_indices_n(&mut [&mut va, &mut vb], &mut indices);
        assert_eq!(va, vec![4, 3, 2, 1]);
        assert_eq!(vb, vec![8, 7, 6, 5]);
    }

    #[test]
    fn test_argsort_by() {
        let va = vec![1, 2, 3, 4];
        let indices = argsort_by(&va, |x| *x);
        assert_eq!(indices, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_sort_multiple_slices() {
        let mut tof_indices: Vec<u32> = vec![1, 2, 3, 4, 1, 2];
        let mut scan_numbers: Vec<usize> = vec![0, 0, 1, 1, 2, 2];
        let mut intensities: Vec<u32> = vec![10, 20, 30, 40, 50, 60];
        let mut indices = argsort_by(&tof_indices, |x| *x);
        sort_by_indices_multi!(
            &mut indices,
            &mut tof_indices,
            &mut scan_numbers,
            &mut intensities
        );
        assert_eq!(tof_indices, vec![1, 1, 2, 2, 3, 4]);
        assert_eq!(scan_numbers, vec![0, 2, 0, 2, 1, 1]);
        assert_eq!(intensities, vec![10, 50, 20, 60, 30, 40]);
    }
}

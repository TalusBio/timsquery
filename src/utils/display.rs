pub fn glimpse_vec<T: std::fmt::Debug>(v: &[T]) -> String {
    // Display short slices as [1,2,3], long slices as
    // [1,2,3 ... 123,124,125] (len = 125)
    let len = v.len();
    let mut out = String::new();
    if len > 10 {
        out.push_str("[");
        out.push_str(
            &v[..3]
                .iter()
                .map(|x| format!("{:?}, ", x))
                .collect::<String>(),
        );
        out.push_str("... ");
        out.push_str(
            &v[len - 3..]
                .iter()
                .map(|x| format!("{:?}, ", x))
                .collect::<String>(),
        );
        out.push_str(&format!("] len = {}", len));
    } else {
        out.push_str("[");
        out.push_str(&v.iter().map(|x| format!("{:?}, ", x)).collect::<String>());
        out.push_str("]");
    }
    out
}

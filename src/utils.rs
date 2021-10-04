pub fn sorted_vec_intersection_count<T: Ord>(a: &[T], b: &[T]) -> usize {
    let mut count = 0;
    let mut b_iter = b.iter();
    if let Some(mut current_b) = b_iter.next() {
        for current_a in a {
            while current_b < current_a {
                current_b = match b_iter.next() {
                    Some(current_b) => current_b,
                    None => return count,
                };
            }
            if current_a == current_b {
                count += 1;
            }
        }
    }
    count
}

use std::{cmp::Ordering, fmt::Debug};

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

pub fn compare_vecs_w_order<T: PartialEq>(a: &[T], b: &[T]) -> bool {
    if a.len() != b.len() {
        return false;
    }

    for i in 0..a.len() {
        if a[i] != b[i] {
            return false;
        }
    }

    true
}

pub fn compare_vecs_no_order<T: Ord + PartialEq + Clone>(a: &Vec<T>, b: &Vec<T>) -> bool {
    if a.len() != b.len() {
        return false;
    }

    let mut a = a.clone();
    a.sort();
    let mut b = b.clone();
    b.sort();
    for i in 0..a.len() {
        if a[i] != b[i] {
            println!("0.2");
            return false;
        }
    }

    true
}

pub fn compare_vecs_w_order_customeq<T, E: Fn(&T, &T) -> bool>(a: &[T], b: &[T], equal: &E) -> bool {
    if a.len() != b.len() {
        return false;
    }

    for i in 0..a.len() {
        if !equal(&a[i], &b[i]) {
            return false;
        }
    }

    true
}

pub fn compare_vecs_no_order_customeq<T: Clone + Debug, F: Fn(&T, &T) -> Ordering, E: Fn(&T, &T) -> bool>(a: &Vec<T>, b: &Vec<T>, compare: &F, equal: &E) -> bool {
    if a.len() != b.len() {
        println!("1");
        return false;
    }

    let mut a = a.clone();
    a.sort_by(compare);
    let mut b = b.clone();
    b.sort_by(compare);
    for i in 0..a.len() {
        if !equal(&a[i], &b[i]) {
            println!("at[{}]: {:?} vs {:?}", i, a[i], b[i]);
            return false;
        }
    }

    true
}
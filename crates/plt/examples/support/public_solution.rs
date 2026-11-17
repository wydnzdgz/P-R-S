// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Public modular elimination for defensive parameter experiments.
//! The caller supplies no secret, and changing column order selects different
//! sparse solutions. Failure of these attempts is not a hardness estimate.

use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    traits::{MatrixDimensions, MatrixGetEntry, MatrixSetEntry},
};
use std::error::Error;

pub fn public_solution(a: &MatZq, target: &MatZq) -> Result<MatZ, Box<dyn Error>> {
    let order: Vec<_> = (0..a.get_num_columns() as usize).collect();
    public_solution_in_order(a, target, &order)
}

pub fn public_solution_in_order(
    a: &MatZq,
    target: &MatZq,
    order: &[usize],
) -> Result<MatZ, Box<dyn Error>> {
    let q = i64::try_from(Z::from(a.get_mod()))?;
    let (rows, columns) = (a.get_num_rows() as usize, a.get_num_columns() as usize);
    if target.get_num_rows() != rows as i64
        || target.get_num_columns() != 1
        || target.get_mod() != a.get_mod()
        || !Z::from(q).is_prime()
    {
        return Err("expected compatible target and prime modulus".into());
    }
    let mut sorted = order.to_vec();
    sorted.sort_unstable();
    if sorted != (0..columns).collect::<Vec<_>>() {
        return Err("column order must be a permutation".into());
    }
    let mut augmented = vec![vec![0i64; columns + 1]; rows];
    for (i, row) in augmented.iter_mut().enumerate() {
        for (j, value) in row[..columns].iter_mut().enumerate() {
            let entry: Z = a.get_entry(i, order[j])?;
            *value = i64::try_from(entry)?;
        }
        let value: Z = target.get_entry(i, 0)?;
        row[columns] = i64::try_from(value)?;
    }
    let multiply =
        |x: i64, y: i64| ((i128::from(x) * i128::from(y)).rem_euclid(i128::from(q))) as i64;
    let mut pivots = Vec::new();
    for column in 0..columns {
        let row = pivots.len();
        let Some(pivot) = (row..rows).find(|&i| augmented[i][column] != 0) else {
            continue;
        };
        augmented.swap(row, pivot);
        let (mut inverse, mut base, mut exponent) = (1, augmented[row][column], q - 2);
        while exponent > 0 {
            if exponent % 2 == 1 {
                inverse = multiply(inverse, base);
            }
            base = multiply(base, base);
            exponent /= 2;
        }
        for value in &mut augmented[row][column..=columns] {
            *value = multiply(*value, inverse);
        }
        let pivot_values = augmented[row][column..=columns].to_vec();
        for (i, other) in augmented.iter_mut().enumerate() {
            if i != row {
                let coefficient = other[column];
                for (value, pivot_value) in other[column..=columns].iter_mut().zip(&pivot_values) {
                    *value = (i128::from(*value)
                        - i128::from(coefficient) * i128::from(*pivot_value))
                    .rem_euclid(i128::from(q)) as i64;
                }
            }
        }
        pivots.push(column);
        if pivots.len() == rows {
            break;
        }
    }
    if pivots.len() != rows {
        return Err("public matrix is not full row rank".into());
    }
    let mut x = MatZ::new(columns, 1);
    for (row, column) in pivots.into_iter().enumerate() {
        let value = augmented[row][columns];
        x.set_entry(
            order[column],
            0,
            if value > q / 2 { value - q } else { value },
        )?;
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn permuted_pivots_preserve_congruence_and_sparse_bound_at_large_modulus() {
        let q = 2_305_843_009_213_693_951i64;
        let a = MatZq::from((&MatZ::from_str("[[1, 0, 2, 3],[0, 1, 5, 7]]").unwrap(), q));
        let y = MatZq::from((
            &MatZ::from_str("[[2305843009213693950],[123456789]]").unwrap(),
            q,
        ));
        for order in [vec![0, 1, 2, 3], vec![3, 2, 1, 0], vec![2, 0, 3, 1]] {
            let x = public_solution_in_order(&a, &y, &order).unwrap();
            assert_eq!(&a * &x, y);
            assert!(x.norm_eucl_sqrd().unwrap() <= Z::from(q / 2).pow(2).unwrap() * 2);
        }
        assert_eq!(&a * public_solution(&a, &y).unwrap(), y);
    }

    #[test]
    fn rejects_invalid_orders_and_rank_deficiency() {
        let a = MatZq::new(2, 3, 257);
        let y = MatZq::new(2, 1, 257);
        assert!(public_solution_in_order(&a, &y, &[0, 0, 2]).is_err());
        assert!(public_solution_in_order(&a, &y, &[0, 1, 8]).is_err());
        assert!(public_solution(&a, &y).is_err());
    }

    use qfall_math::traits::Pow;
    use std::str::FromStr;
}

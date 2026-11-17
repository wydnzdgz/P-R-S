// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Vandermonde target decomposition for arbitrary validated t-of-k configurations.

use super::{PltConfig, PltError};
use qfall_math::{
    integer::Z,
    integer_mod_q::MatZq,
    traits::{
        MatrixDimensions, MatrixGetEntry, MatrixGetSubmatrix, MatrixSetEntry, MatrixSetSubmatrix,
    },
};

/// Returns v_j tensor matrix, v_j=(1,a_j,...,a_j^(t-1))^T with a_j=j+1.
/// `matrix` has n rows, one or more columns and modulus q from config.
pub fn participant_image(
    config: &PltConfig,
    participant: usize,
    matrix: &MatZq,
) -> Result<MatZq, PltError> {
    if participant >= config.participants() {
        return Err(PltError::InvalidParticipant);
    }
    if matrix.get_num_rows() != config.n() || matrix.get_mod() != config.modulus() {
        return Err(PltError::InvalidTarget);
    }
    let mut out = MatZq::new(
        config.target_rows(),
        matrix.get_num_columns(),
        config.modulus(),
    );
    let mut power = 1;
    for block in 0..config.threshold() {
        let scaled = matrix * &Z::from(power);
        out.set_submatrix(
            block as i64 * config.n(),
            0,
            &scaled,
            0,
            0,
            config.n() - 1,
            matrix.get_num_columns() - 1,
        )
        .unwrap();
        power = next_power(power, participant, config.modulus());
    }
    Ok(out)
}

/// Computes all z_j in the supplied recovery-set order using (V_T^-1 tensor I_n)y.
/// Uses blocks so n>1 and permuted sets share the same algorithm.
/// Rejects wrong cardinality, repeated/invalid indices, target shape or modulus.
pub fn decompose_target(
    config: &PltConfig,
    recovery_set: &[usize],
    target: &MatZq,
) -> Result<Vec<MatZq>, PltError> {
    config.validate_set(recovery_set)?;
    config.validate_target(target)?;
    let t = config.threshold() as i64;
    let mut v = MatZq::new(t, t, config.modulus());
    for (column, &participant) in recovery_set.iter().enumerate() {
        let mut power = 1;
        for row in 0..t {
            v.set_entry(row, column, power).unwrap();
            power = next_power(power, participant, config.modulus());
        }
    }
    let inverse = v.inverse().ok_or(PltError::InvalidRecoverySet)?;
    let blocks: Vec<_> = (0..t)
        .map(|i| {
            target
                .get_submatrix(i * config.n(), (i + 1) * config.n() - 1, 0, 0)
                .unwrap()
        })
        .collect();
    let mut components = Vec::with_capacity(config.threshold());
    for row in 0..t {
        let mut z = MatZq::new(config.n(), 1, config.modulus());
        for column in 0..t {
            let coefficient: Z = inverse.get_entry(row, column).unwrap();
            z += &blocks[column as usize] * &coefficient;
        }
        components.push(z);
    }
    Ok(components)
}

/// Modular products fit i128 because both factors are positive i64 values.
fn next_power(power: i64, participant: usize, modulus: i64) -> i64 {
    (i128::from(power) * (participant as i128 + 1) % i128::from(modulus)) as i64
}

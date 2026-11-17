// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Fallible entry points shared by the legacy PSF API and PLT.

use super::{PSFPerturbation, gen_trapdoor, randomized_nearest_plane_gadget, short_basis_gadget};
use crate::sample::randomness::fresh_seed;
use qfall_math::{
    error::MathError,
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::{MatQ, Q},
    traits::{Concatenate, MatrixDimensions, MatrixGetEntry, MatrixSetEntry, Pow},
};
use rand::{SeedableRng, rngs::StdRng};
use rand_distr::{Distribution, StandardNormal};

type Trapdoor = (MatZ, MatQ, (MatZ, MatQ));

impl PSFPerturbation {
    /// Generates a trapdoor, returning parameter, sampling or numerical errors.
    /// This does not turn allocation failures or foreign-library aborts into errors.
    pub fn try_trap_gen(&self) -> Result<(MatZq, Trapdoor), MathError> {
        let n = i64::try_from(&self.gp.n).map_err(|_| dimension_error())?;
        let k = i64::try_from(&self.gp.k).map_err(|_| dimension_error())?;
        let bar = i64::try_from(&self.gp.m_bar).map_err(|_| dimension_error())?;
        if n < 1 || k < 1 || bar < 1 || self.gp.base < 2 || self.r <= 0 || self.s <= 0 {
            return Err(dimension_error());
        }
        n.checked_mul(k)
            .and_then(|w| w.checked_add(bar))
            .ok_or_else(dimension_error)?;
        let a_bar = MatZq::sample_uniform_seeded(n, bar, &self.gp.q, fresh_seed()?);
        let tag = MatZq::identity(n, n, &self.gp.q);
        let (a, r) = gen_trapdoor(&self.gp, &a_bar, &tag)?;
        let sigma = self.s.pow(2)? * MatQ::identity(a.get_num_columns(), a.get_num_columns());
        let factor = self.try_compute_sqrt_sigma_2(&r, &sigma)?;
        let basis = short_basis_gadget(&self.gp);
        let gso = MatQ::from(&basis).gso();
        Ok((a, (r, factor, (basis, gso))))
    }

    /// Computes the perturbation factor with shape, symmetry, finiteness,
    /// positive-pivot and relative reconstruction-residual checks (tolerance 1e-8).
    /// Input follows the existing convention: target squared width is r^2*mat_sigma.
    pub fn try_compute_sqrt_sigma_2(
        &self,
        mat_r: &MatZ,
        mat_sigma: &MatQ,
    ) -> Result<MatQ, MathError> {
        let dim = mat_r
            .get_num_rows()
            .checked_add(mat_r.get_num_columns())
            .ok_or_else(dimension_error)?;
        if self.r <= 0
            || self.gp.base < 2
            || mat_sigma.get_num_rows() != dim
            || mat_sigma.get_num_columns() != dim
            || !mat_sigma.is_symmetric()
        {
            return Err(dimension_error());
        }
        let b = mat_r.concat_vertical(&MatZ::identity(
            mat_r.get_num_columns(),
            mat_r.get_num_columns(),
        ))?;
        let sigma: MatQ = (1.0 / (2.0 * Q::PI))
            * self.r.pow(2)?
            * (mat_sigma
                - (self.gp.base.pow(2)? + 1) * &b * b.transpose()
                - MatQ::identity(dim, dim));
        checked_cholesky(&sigma)
    }

    /// Samples a single conditioned preimage through the same MP12 core as samp_p,
    /// returning recoverable input/sampling errors instead of unwrapping them.
    pub fn try_samp_p(&self, a: &MatZq, td: &Trapdoor, target: &MatZq) -> Result<MatZ, MathError> {
        let (r, factor, (basis, gso)) = td;
        let rows = i64::try_from(&self.gp.n).map_err(|_| dimension_error())?;
        let gadget_columns: Z = &self.gp.n * &self.gp.k;
        let dim = a.get_num_columns();
        if self.r <= 0
            || !f64::from(&self.r).is_finite()
            || self.gp.base < 2
            || self.gp.k < 1
            || self.gp.base.pow(&self.gp.k)? < a.get_mod()
            || a.get_num_rows() != rows
            || target.get_num_rows() != rows
            || target.get_num_columns() != 1
            || a.get_mod() != target.get_mod()
            || a.get_mod() != self.gp.q
            || r.get_num_columns() != gadget_columns
            || r.get_num_rows() + r.get_num_columns() != dim
            || factor.get_num_rows() != dim
            || factor.get_num_columns() != dim
            || basis.get_num_rows() != gadget_columns
            || basis.get_num_columns() != gadget_columns
            || gso.get_num_rows() != gadget_columns
            || gso.get_num_columns() != gadget_columns
        {
            return Err(dimension_error());
        }
        for i in 0..dim {
            for j in 0..dim {
                let value: Q = factor.get_entry(i, j)?;
                if !f64::from(&value).is_finite() {
                    return Err(numerical_error());
                }
            }
        }
        let p = secure_perturbation(factor, &self.r)?;
        let v = target - a * &p;
        let z = randomized_nearest_plane_gadget(self, &v, basis, gso)?;
        let b = r.concat_vertical(&MatZ::identity(r.get_num_columns(), r.get_num_columns()))?;
        Ok(p + b * z)
    }
}

/// Replaces qfall-math's probability/Xorshift continuous branch on the MP12 path.
/// Each call seeds a CSPRNG with 256 OS bits; rounding gets a separate fresh seed.
fn secure_perturbation(factor: &MatQ, rounding: &Q) -> Result<MatZ, MathError> {
    let mut rng = StdRng::from_seed(fresh_seed()?);
    let mut standard = MatQ::new(factor.get_num_columns(), 1);
    for i in 0..factor.get_num_columns() {
        let value: f64 = StandardNormal.sample(&mut rng);
        if !value.is_finite() {
            return Err(numerical_error());
        }
        standard.set_entry(i, 0, value)?;
    }
    (factor * standard).randomized_rounding_seeded(rounding, fresh_seed()?)
}

fn dimension_error() -> MathError {
    MathError::MismatchingMatrixDimension(
        "invalid MP12 parameters, target or trapdoor shape".into(),
    )
}
fn numerical_error() -> MathError {
    MathError::ConversionError("non-finite or inaccurate MP12 numerical computation".into())
}

/// Rejects unsupported numerical inputs before constructing a rational from NaN/Inf.
fn checked_cholesky(matrix: &MatQ) -> Result<MatQ, MathError> {
    if !matrix.is_square() || !matrix.is_symmetric() {
        return Err(dimension_error());
    }
    let dim = matrix.get_num_rows() as usize;
    let mut a = vec![vec![0.0; dim]; dim];
    let mut l = vec![vec![0.0; dim]; dim];
    let mut scale = 0.0f64;
    for (i, row) in a.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            let entry: Q = matrix.get_entry(i, j)?;
            *value = f64::from(&entry);
            if !value.is_finite() {
                return Err(numerical_error());
            }
            scale = scale.max(value.abs());
        }
    }
    for i in 0..dim {
        for j in 0..=i {
            let dot: f64 = (0..j).map(|k| l[i][k] * l[j][k]).sum();
            let value = a[i][j] - dot;
            if !value.is_finite() {
                return Err(numerical_error());
            }
            l[i][j] = if i == j {
                if value <= 0.0 {
                    return Err(MathError::NonPositive(
                        "MP12 perturbation matrix is not numerically positive definite".into(),
                    ));
                }
                value.sqrt()
            } else {
                value / l[j][j]
            };
            if !l[i][j].is_finite() {
                return Err(numerical_error());
            }
        }
    }
    // Bound forward error amplification by the infinity-norm condition estimate.
    // Solve A^{-1} columnwise through L, avoiding a second generic inversion.
    let mut inverse_row_sums = vec![0.0; dim];
    for column in 0..dim {
        let mut y = vec![0.0; dim];
        for i in 0..dim {
            let rhs = if i == column { 1.0 } else { 0.0 };
            y[i] = (rhs - (0..i).map(|j| l[i][j] * y[j]).sum::<f64>()) / l[i][i];
        }
        let mut z = vec![0.0; dim];
        for i in (0..dim).rev() {
            z[i] = (y[i] - (i + 1..dim).map(|j| l[j][i] * z[j]).sum::<f64>()) / l[i][i];
            inverse_row_sums[i] += z[i].abs();
        }
    }
    let norm = a
        .iter()
        .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0, f64::max);
    let inverse_norm = inverse_row_sums.into_iter().fold(0.0, f64::max);
    let condition = norm * inverse_norm;
    if !condition.is_finite() || condition * f64::EPSILON > 1e-8 {
        return Err(numerical_error());
    }
    let mut out = MatQ::new(dim, dim);
    for i in 0..dim {
        for j in 0..dim {
            let reconstructed: f64 = (0..=i.min(j)).map(|k| l[i][k] * l[j][k]).sum();
            if !reconstructed.is_finite() || (reconstructed - a[i][j]).abs() > 1e-8 * scale {
                return Err(numerical_error());
            }
            out.set_entry(i, j, l[i][j])?;
        }
    }
    Ok(out)
}

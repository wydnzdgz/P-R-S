// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Dimension-aware adapter for the existing integer MP12 sampler.

use super::{PltConfig, PltError, SamplingProfile};
use qfall_math::{
    error::MathError,
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::{MatQ, Q},
    traits::{
        Concatenate, MatrixDimensions, MatrixGetEntry, MatrixGetSubmatrix, MatrixSetSubmatrix, Pow,
    },
};
use qfall_tools::{
    primitive::psf::{PSF, PSFPerturbation},
    sample::g_trapdoor::{gadget_classical::gen_gadget_mat, gadget_parameters::GadgetParameters},
};

type Mp12Trapdoor = <PSFPerturbation as PSF>::Trapdoor;

/// Holds a public matrix, the sampler configuration and its cached trapdoor.
pub(super) struct Mp12Key {
    a: MatZq,
    sampler: PSFPerturbation,
    td: Mp12Trapdoor,
    profile: SamplingProfile,
}

impl Mp12Key {
    pub(super) const ROUNDING: i64 = 4;

    /// Generates the requested dimensions with ternary R. Selects the smallest
    /// power of two s>=64 satisfying s^2 > 5*||[R;I]||_F^2+1 for every ternary R.
    /// This avoids carrying a fixed-size positivity assumption into larger thresholds.
    pub(super) fn generate(config: &PltConfig, rows: i64, columns: i64) -> Result<Self, PltError> {
        let mut gp = GadgetParameters::init_default(rows, config.modulus());
        let gadget_columns = &gp.n * &gp.k;
        gp.m_bar = Z::from(columns) - &gadget_columns;
        let worst_case_norm = (&gp.m_bar + 1) * &gadget_columns;
        let bound = 5 * worst_case_norm + 1;
        let mut s = Z::from(64);
        while &s * &s <= bound {
            s *= 2;
        }
        let sampler = PSFPerturbation {
            gp,
            r: Q::from(Self::ROUNDING),
            s: Q::from(s),
        };
        let (a, td) = sampler.try_trap_gen().map_err(backend_error)?;
        let profile = SamplingProfile {
            rows,
            columns,
            gadget_length: i64::try_from(&sampler.gp.k).map_err(|_| PltError::InvalidParameters)?,
            m_bar: i64::try_from(&sampler.gp.m_bar).map_err(|_| PltError::InvalidParameters)?,
            rounding_parameter: Self::ROUNDING,
            initial_width: &sampler.s * &sampler.r,
            distribution: "MP12 PlusMinusOneZero: P(-1)=1/4, P(0)=1/2, P(1)=1/4",
        };
        let b =
            td.0.concat_vertical(&MatZ::identity(
                td.0.get_num_columns(),
                td.0.get_num_columns(),
            ))
            .map_err(backend_error)?;
        let g = gen_gadget_mat(rows, &sampler.gp.k, &sampler.gp.base);
        if &a * &b != MatZq::from((&g, config.modulus())) {
            return Err(PltError::AlgebraicMismatch);
        }
        for i in 0..td.0.get_num_rows() {
            for j in 0..td.0.get_num_columns() {
                let value: Z = td.0.get_entry(i, j).map_err(backend_error)?;
                if !(Z::MINUS_ONE..=Z::ONE).contains(&value) {
                    return Err(PltError::AlgebraicMismatch);
                }
            }
        }
        Ok(Self {
            a,
            sampler,
            td,
            profile,
        })
    }

    /// Borrows the public matrix.
    pub(super) fn a(&self) -> &MatZq {
        &self.a
    }

    /// Returns non-secret algorithm parameters; never returns the secret R.
    pub(super) fn profile(&self) -> &SamplingProfile {
        &self.profile
    }

    /// Independent reconstruction of the total target squared width, including
    /// randomized rounding and gadget contributions, in the paper convention.
    pub(super) fn diagnostics(&self, sigma: &MatQ) -> Result<(f64, f64), PltError> {
        let factor = self.perturbation_factor(sigma)?;
        let dim = self.a.get_num_columns();
        let b = self
            .td
            .0
            .concat_vertical(&MatZ::identity(
                self.td.0.get_num_columns(),
                self.td.0.get_num_columns(),
            ))
            .map_err(backend_error)?;
        let base_squared = self.sampler.gp.base.pow(2).map_err(backend_error)? + 1;
        let r_squared = self.sampler.r.pow(2).map_err(backend_error)?;
        let perturbation: MatQ = (1.0 / (2.0 * Q::PI))
            * (sigma
                - &r_squared * &base_squared * &b * b.transpose()
                - &r_squared * MatQ::identity(dim, dim));
        let reconstructed: MatQ = (2.0 * Q::PI) * &factor * factor.transpose()
            + &r_squared * &base_squared * &b * b.transpose()
            + &r_squared * MatQ::identity(dim, dim);
        let difference = &reconstructed - sigma;
        let mut error = 0.0f64;
        let mut scale = 0.0f64;
        for i in 0..dim {
            for j in 0..dim {
                let delta: Q = difference.get_entry(i, j).map_err(backend_error)?;
                let expected: Q = sigma.get_entry(i, j).map_err(backend_error)?;
                error = error.max(f64::from(&delta).abs());
                scale = scale.max(f64::from(&expected).abs());
            }
        }
        let relative = error / scale;
        if !relative.is_finite() || relative > 1e-8 {
            return Err(PltError::NumericalFailure);
        }
        Ok((super::assessment::condition_inf(&perturbation)?, relative))
    }

    /// Samples independent columns at the configured spherical width r*s.
    /// Every column is checked and any failure aborts without retry.
    pub(super) fn sample_columns(&self, targets: &MatZq) -> Result<MatZ, PltError> {
        if targets.get_num_rows() != self.a.get_num_rows() || targets.get_mod() != self.a.get_mod()
        {
            return Err(PltError::InvalidTarget);
        }
        let mut result = MatZ::new(self.a.get_num_columns(), targets.get_num_columns());
        for column in 0..targets.get_num_columns() {
            let target = targets.get_column(column).map_err(backend_error)?;
            let sample = self
                .sampler
                .try_samp_p(&self.a, &self.td, &target)
                .map_err(backend_error)?;
            if &self.a * &sample != target {
                return Err(PltError::AlgebraicMismatch);
            }
            if sample.norm_eucl_sqrd().map_err(backend_error)?
                > self.profile.initial_width.pow(2).map_err(backend_error)?
                    * self.a.get_num_columns()
            {
                return Err(PltError::NormExceeded);
            }
            result
                .set_column(column, &sample, 0)
                .map_err(backend_error)?;
        }
        Ok(result)
    }

    /// Samples with paper squared-width matrix Sigma, after its positivity bound
    /// has been checked by the caller. The existing API multiplies its input by r^2,
    /// so pass Sigma/r^2. It handles the conversion to continuous variance 1/(2*pi).
    /// Never apply PSF::f_a/check_domain here: those assume the initial spherical width.
    pub(super) fn sample_with_squared_width(
        &self,
        target: &MatZq,
        sigma: &MatQ,
    ) -> Result<MatZ, PltError> {
        let factor = self.perturbation_factor(sigma)?;
        let td = (self.td.0.clone(), factor, self.td.2.clone());
        self.sampler
            .try_samp_p(&self.a, &td, target)
            .map_err(backend_error)
    }

    /// Converts a squared width matrix to the continuous perturbation factor.
    fn perturbation_factor(&self, sigma: &MatQ) -> Result<MatQ, PltError> {
        let normalized = Q::from((1, Self::ROUNDING.pow(2))) * sigma;
        self.sampler
            .try_compute_sqrt_sigma_2(&self.td.0, &normalized)
            .map_err(backend_error)
    }

    /// Returns ||[R;I]||_F^2 exactly, including the identity rows.
    pub(super) fn full_trapdoor_norm_squared(&self) -> Z {
        frobenius_squared(&self.td.0) + self.td.0.get_num_columns()
    }
}

/// Computes the exact sum of squared integer entries, without f64 conversion.
pub(super) fn frobenius_squared(matrix: &MatZ) -> Z {
    let mut sum = Z::ZERO;
    for row in 0..matrix.get_num_rows() {
        for column in 0..matrix.get_num_columns() {
            let value: Z = matrix.get_entry(row, column).unwrap();
            sum += &value * &value;
        }
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use qfall_math::traits::{Concatenate, MatrixSetEntry};

    /// Checks the r^2 and 2*pi conversion independently using a non-diagonal width.
    /// This checks the perturbation covariance, not the sampler's statistical distance.
    #[test]
    fn non_spherical_perturbation_has_the_requested_covariance() {
        let config = PltConfig::new(2, 3).unwrap();
        let key = Mp12Key::generate(&config, 1, 10).unwrap();
        let mut sigma: MatQ = 3200 * MatQ::identity(10, 10);
        sigma.set_entry(0, 1, 16).unwrap();
        sigma.set_entry(1, 0, 16).unwrap();
        let factor = key.perturbation_factor(&sigma).unwrap();
        let b = key.td.0.concat_vertical(&MatZ::identity(9, 9)).unwrap();
        let expected: MatQ = &sigma - 80 * &b * b.transpose() - 16 * MatQ::identity(10, 10);
        let actual: MatQ = (2.0 * Q::PI) * &factor * factor.transpose();
        let difference = &actual - &expected;
        for row in 0..10 {
            for column in 0..10 {
                let error: Q = difference.get_entry(row, column).unwrap();
                assert!(f64::from(&error).abs() < 1e-8);
            }
        }
    }
}

/// Preserves recoverable backend error categories without catching arbitrary panics.
pub(super) fn backend_error(error: MathError) -> PltError {
    match error {
        MathError::NonPositive(_) => PltError::CovarianceNotPositiveDefinite,
        MathError::ConversionError(_) => PltError::NumericalFailure,
        other => PltError::BackendFailure(other.to_string()),
    }
}

#[cfg(test)]
mod failure_tests {
    use super::*;
    use qfall_math::traits::{MatrixSetEntry, Pow};

    #[test]
    fn covariance_failures_are_results_not_panics() {
        let config = PltConfig::new(2, 3).unwrap();
        let key = Mp12Key::generate(&config, 1, 10).unwrap();
        assert_eq!(
            key.perturbation_factor(&MatQ::new(10, 10)),
            Err(PltError::CovarianceNotPositiveDefinite)
        );
        let mut asym: MatQ = 3200 * MatQ::identity(10, 10);
        asym.set_entry(0, 1, 1).unwrap();
        assert!(matches!(
            key.perturbation_factor(&asym),
            Err(PltError::BackendFailure(_))
        ));
        assert!(matches!(
            key.perturbation_factor(&MatQ::identity(9, 9)),
            Err(PltError::BackendFailure(_))
        ));
        let huge: MatQ = Q::from(Z::from(2).pow(2048).unwrap()) * MatQ::identity(10, 10);
        assert_eq!(
            key.perturbation_factor(&huge),
            Err(PltError::NumericalFailure)
        );
    }

    #[test]
    fn backend_and_issued_norm_errors_propagate() {
        let config = PltConfig::new(2, 3).unwrap();
        let mut key = Mp12Key::generate(&config, 1, 10).unwrap();
        assert_eq!(
            key.sample_columns(&MatZq::new(2, 1, 257)),
            Err(PltError::InvalidTarget)
        );
        assert_eq!(
            key.sample_columns(&MatZq::new(1, 1, 17)),
            Err(PltError::InvalidTarget)
        );
        let target = MatZq::identity(1, 1, 257);
        key.profile.initial_width = Q::ZERO; // Fault injection; actual sampler remains unchanged.
        assert_eq!(key.sample_columns(&target), Err(PltError::NormExceeded));
        key.sampler.r = Q::ZERO;
        assert!(matches!(
            key.sample_with_squared_width(&target, &MatQ::identity(10, 10)),
            Err(PltError::BackendFailure(_))
        ));
    }
}

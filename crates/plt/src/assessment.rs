// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Explicit numerical evidence and fail-closed parameter assessment.

use super::*;
use qfall_math::traits::MatrixGetEntry;

/// Engineering numerical checks; not a statistical-distance or security certificate.
#[derive(Debug, Serialize)]
pub struct SamplingDiagnostics {
    pub dimension: i64,
    pub gram_condition_inf: f64,
    pub covariance_condition_inf: f64,
    pub perturbation_condition_inf: f64,
    pub total_squared_width_relative_residual: f64,
    pub pullback_metric_exact: bool,
}

/// Conditions not established remain false. No current configuration is certified.
#[derive(Debug, Serialize)]
pub struct ParameterAssessment {
    /// Necessary comparison when claiming the paper's Table 2 instantiation (integer g=1).
    pub table2_issuance_spectral_upper_bound_squared: Q,
    pub table2_issuance_upper_bound_satisfied: bool,
    pub table2_m_greater_than_participants: bool,
    pub public_linear_solution_norm_bound_squared: Z,
    pub accepted_full_norm_bound_squared: Q,
    pub accepts_public_linear_solutions: bool,
    pub trapdoor_uniformity_established: bool,
    pub smoothing_conditions_established: bool,
    pub concrete_hardness_established: bool,
    pub protocol_security_approved: bool,
    pub reasons: Vec<String>,
}

impl PublicParameters {
    /// Compares beta^2 with a constructive bound for public modular elimination.
    /// A full-row-rank A admits a solution supported on at most n*t coordinates;
    /// centered lifts have squared norm at most n*t*floor(q/2)^2.
    pub fn assess_parameters(&self, sigma_squared: &Q) -> Result<ParameterAssessment, PltError> {
        validate_width(sigma_squared)?;
        let half = Z::from(self.config.modulus() / 2);
        let bound = &half * &half * self.config.target_rows();
        let t = Z::from(self.config.threshold() as i64);
        let beta = &t * &t * sigma_squared * self.config.preimage_rows();
        let accepts = beta >= bound;
        let q = Z::from(self.config.modulus());
        let upper = Q::from(&q * &q) / self.config.m();
        let initial = &self.master_profile.initial_width;
        let issuance_ok = initial * initial <= upper;
        let m_greater = self.config.m() > self.config.participants() as i64;
        let mut reasons=vec!["Trapdoor uniformity, smoothing error budget and concrete hardness are not established for this instance.".into()];
        if accepts {
            reasons.push("The accepted norm bound permits public Gaussian-elimination solutions without any trapdoor or threshold shares.".into());
        }
        if !issuance_ok {
            reasons.push("The selected issuance width violates the Table 2 necessary spectral upper bound q/sqrt(m) for the integer instantiation.".into());
        }
        if !m_greater {
            reasons.push(
                "m<=k violates the nondegenerate m-k condition of the cited Table 2 instantiation."
                    .into(),
            );
        }
        Ok(ParameterAssessment {
            table2_issuance_spectral_upper_bound_squared: upper,
            table2_issuance_upper_bound_satisfied: issuance_ok,
            table2_m_greater_than_participants: m_greater,
            public_linear_solution_norm_bound_squared: bound,
            accepted_full_norm_bound_squared: beta,
            accepts_public_linear_solutions: accepts,
            trapdoor_uniformity_established: false,
            smoothing_conditions_established: false,
            concrete_hardness_established: false,
            protocol_security_approved: false,
            reasons,
        })
    }

    /// Optional gate for callers needing actual parameter approval. Experiments keep
    /// working; no unassessed or demonstrably weak instance is silently approved.
    pub fn require_protocol_parameters(&self, sigma_squared: &Q) -> Result<(), PltError> {
        if !self
            .assess_parameters(sigma_squared)?
            .protocol_security_approved
        {
            return Err(PltError::UnapprovedParameters);
        }
        Ok(())
    }
}

impl PartialTrapdoor {
    /// Checks the actual Sigma/Gram relation and MP12 perturbation reconstruction.
    /// Does not return U or the auxiliary secret trapdoor in the report.
    pub fn sampling_diagnostics(&self, sigma_squared: &Q) -> Result<SamplingDiagnostics, PltError> {
        validate_width(sigma_squared)?;
        if sigma_squared < &self.minimum_sigma_squared {
            return Err(PltError::WidthTooSmall);
        }
        let gram = MatQ::from(&(&self.u.transpose() * &self.u));
        let covariance = sigma_squared * &self.gram_inverse;
        let exact =
            &gram * &covariance == sigma_squared * MatQ::identity(self.config.m(), self.config.m());
        let (perturbation_condition_inf, total_squared_width_relative_residual) =
            self.auxiliary.diagnostics(&covariance)?;
        Ok(SamplingDiagnostics {
            dimension: self.config.m(),
            gram_condition_inf: condition_inf(&gram)?,
            covariance_condition_inf: condition_inf(&covariance)?,
            perturbation_condition_inf,
            total_squared_width_relative_residual,
            pullback_metric_exact: exact,
        })
    }
}

/// Independent exact rational inversion followed by finite float norm evaluation.
pub(super) fn condition_inf(matrix: &MatQ) -> Result<f64, PltError> {
    let inverse = matrix
        .inverse()
        .ok_or(PltError::CovarianceNotPositiveDefinite)?;
    let condition = float_norm_inf(matrix)? * float_norm_inf(&inverse)?;
    if !condition.is_finite() {
        return Err(PltError::NumericalFailure);
    }
    Ok(condition)
}
fn float_norm_inf(matrix: &MatQ) -> Result<f64, PltError> {
    let mut maximum = 0.0f64;
    for i in 0..matrix.get_num_rows() {
        let mut sum = 0.0;
        for j in 0..matrix.get_num_columns() {
            let value: Q = matrix.get_entry(i, j).map_err(backend_error)?;
            sum += f64::from(&value).abs();
        }
        maximum = maximum.max(sum);
    }
    if !maximum.is_finite() {
        return Err(PltError::NumericalFailure);
    }
    Ok(maximum)
}

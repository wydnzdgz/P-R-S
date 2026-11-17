// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Configurable integer-lattice PLT prototype of Fig. 2 in Albrecht et al.,
//! "Partial Lattice Trapdoors: How to Split Lattice Trapdoors, Literally".
//!
//! One algorithm implementation handles every validated `PltConfig`.
//! It reuses the existing MP12 implementation and its ternary trapdoor distribution.
//! Parameter validation and functional tests do not establish cryptographic security.
//! Gaussian matrices use exp(-pi * x^T * Sigma^-1 * x); `sigma_squared`
//! denotes the paper's squared width, not ordinary variance.
//!
//! ```
//! use plt::{PltConfig, setup};
//! use qfall_math::traits::MatrixDimensions;
//! let config = PltConfig::new(3, 5).unwrap();
//! let (public, _master) = setup(&config).unwrap();
//! assert_eq!(public.a().get_num_rows(), 3);
//! ```

mod assessment;
mod mp12_adapter;
mod params;
mod vandermonde;

use self::mp12_adapter::{Mp12Key, backend_error, frobenius_squared};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::{MatQ, Q},
    traits::MatrixDimensions,
};
use serde::Serialize;
use std::{error::Error, fmt};

pub use assessment::{ParameterAssessment, SamplingDiagnostics};
pub use params::PltConfig;
pub use vandermonde::{decompose_target, participant_image};

/// Failures checked at the boundary of the prototype algorithms.
#[derive(Debug, PartialEq, Eq)]
pub enum PltError {
    /// The instance has not been approved for protocol-security claims.
    UnapprovedParameters,
    /// Invalid threshold, modulus, lattice dimensions or dimension overflow.
    InvalidParameters,
    /// Participant indices must be in 0..k.
    InvalidParticipant,
    /// A recovery set must contain exactly t distinct valid participants.
    InvalidRecoverySet,
    /// The target must be an n*t-by-1 matrix with the configured modulus.
    InvalidTarget,
    /// The supplied trapdoor belongs to another master public matrix.
    KeyMismatch,
    /// A sampled U has deficient rank over the rationals; PTrapGen aborts.
    RankDeficient,
    /// The requested squared width is below the conservative prototype bound.
    WidthTooSmall,
    /// Shares have inconsistent participants, recovery sets, targets, widths or keys.
    InvalidShares,
    /// An integer vector has the wrong dimension or is not a column vector.
    InvalidPreimage,
    /// An issued column, partial or full preimage exceeds the configured norm bound.
    NormExceeded,
    /// The supplied covariance cannot be factored with positive numerical pivots.
    CovarianceNotPositiveDefinite,
    /// The supported finite-precision range or factorization residual was exceeded.
    NumericalFailure,
    /// A modular trapdoor, auxiliary-preimage or reconstruction relation failed.
    AlgebraicMismatch,
    /// A fallible dependency operation failed; no sample is silently retried.
    BackendFailure(String),
}

impl fmt::Display for PltError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::UnapprovedParameters => {
                "parameters are experimental or fail the required security conditions"
            }
            Self::InvalidParameters => {
                "require n>0, prime q, 1<=t<=k<q and representable dimensions"
            }
            Self::InvalidParticipant => "participant index must be less than k",
            Self::InvalidRecoverySet => "recovery requires exactly t distinct participants in 0..k",
            Self::InvalidTarget => "target must be an n*t-by-1 matrix with the configured modulus",
            Self::KeyMismatch => "trapdoor and master public matrix do not match",
            Self::RankDeficient => "PTrapGen aborted: U does not have full rational column rank",
            Self::WidthTooSmall => "sigma_squared is below the conservative sampling bound",
            Self::InvalidPreimage => "integer preimage must be a 2*m*t-by-1 vector",
            Self::NormExceeded => "integer preimage exceeds the configured squared norm bound",
            Self::CovarianceNotPositiveDefinite => {
                "perturbation covariance is not numerically positive definite"
            }
            Self::NumericalFailure => "non-finite or inaccurate numerical computation",
            Self::AlgebraicMismatch => "modular algorithm invariant failed",
            Self::BackendFailure(message) => return write!(f, "PLT backend failed: {message}"),
            Self::InvalidShares => "partial preimages do not belong to the same reconstruction",
        })
    }
}

impl Error for PltError {}

/// Public record of the selected MP12 engineering instantiation, without secret matrices.
#[derive(Clone, Debug, Serialize)]
pub struct SamplingProfile {
    pub rows: i64,
    pub columns: i64,
    pub gadget_length: i64,
    pub m_bar: i64,
    pub rounding_parameter: i64,
    pub initial_width: Q,
    pub distribution: &'static str,
}

/// Master public matrix; sufficient for the online phase together with partial trapdoors.
pub struct PublicParameters {
    master_profile: SamplingProfile,
    config: PltConfig,
    a: MatZq,
}

impl PublicParameters {
    /// Returns the exact non-secret MP12 profile used to generate the master matrix.
    pub fn master_sampling_profile(&self) -> &SamplingProfile {
        &self.master_profile
    }

    /// Returns the validated configuration for this public key.
    pub fn config(&self) -> &PltConfig {
        &self.config
    }

    /// Returns A, of size n*t-by-2*m*t modulo q.
    pub fn a(&self) -> &MatZq {
        &self.a
    }
}

/// Master secret used only by PTrapGen, not by online sampling or reconstruction.
pub struct MasterTrapdoor {
    config: PltConfig,
    key: Mp12Key,
}

/// A participant's (U, C, td_C), with a cached exact inverse Gram matrix.
pub struct PartialTrapdoor {
    config: PltConfig,
    participant: usize,
    master_a: MatZq,
    u: MatZ,
    auxiliary: Mp12Key,
    gram_inverse: MatQ,
    minimum_sigma_squared: Q,
}

impl PartialTrapdoor {
    /// Returns the auxiliary key's non-secret generation parameters.
    pub fn auxiliary_sampling_profile(&self) -> &SamplingProfile {
        self.auxiliary.profile()
    }

    /// Returns the participant index, in 0..k.
    pub fn participant(&self) -> usize {
        self.participant
    }

    /// Returns the short integer matrix U, without modular reduction.
    pub fn u(&self) -> &MatZ {
        &self.u
    }

    /// Returns the auxiliary n-by-m public matrix C modulo q.
    pub fn c(&self) -> &MatZq {
        self.auxiliary.a()
    }

    /// Returns a conservative sufficient squared width for this partial trapdoor.
    /// Use the maximum across participants for a common online width.
    /// This is an implementation bound, not the paper's security parameter selection.
    pub fn minimum_sigma_squared(&self) -> &Q {
        &self.minimum_sigma_squared
    }
}

/// A sampled partial preimage together with intermediate values for correctness checks.
pub struct PartialPreimage {
    config: PltConfig,
    participant: usize,
    recovery_set: Vec<usize>,
    master_a: MatZq,
    target: MatZq,
    sigma_squared: Q,
    z: MatZq,
    d: MatZ,
    x: MatZ,
}

impl PartialPreimage {
    /// Returns the participant index.
    pub fn participant(&self) -> usize {
        self.participant
    }

    /// Returns the target component z modulo the configured q.
    pub fn z(&self) -> &MatZq {
        &self.z
    }

    /// Returns the auxiliary short integer preimage d, satisfying C*d = z mod q.
    pub fn d(&self) -> &MatZ {
        &self.d
    }

    /// Returns x = U*d over the integers, without modular reduction.
    pub fn x(&self) -> &MatZ {
        &self.x
    }
}

/// Generates the master A and its MP12 trapdoor for a validated configuration.
/// Initial sampling width is selected conservatively for the configured dimensions.
pub fn setup(config: &PltConfig) -> Result<(PublicParameters, MasterTrapdoor), PltError> {
    let key = Mp12Key::generate(config, config.target_rows(), config.preimage_rows())?;
    Ok((
        PublicParameters {
            master_profile: key.profile().clone(),
            config: config.clone(),
            a: key.a().clone(),
        },
        MasterTrapdoor {
            config: config.clone(),
            key,
        },
    ))
}

/// Generates a partial trapdoor for `participant` using the master secret.
/// Samples each column of U independently with the same master Gaussian width.
/// Returns an error for invalid indices, mismatched keys or deficient rational rank.
pub fn ptrap_gen(
    public: &PublicParameters,
    master: &MasterTrapdoor,
    participant: usize,
) -> Result<PartialTrapdoor, PltError> {
    let config = &public.config;
    if participant >= config.participants() {
        return Err(PltError::InvalidParticipant);
    }
    if public.config != master.config || public.a != *master.key.a() {
        return Err(PltError::KeyMismatch);
    }
    let auxiliary = Mp12Key::generate(config, config.n(), config.m())?;
    let image = participant_image(config, participant, auxiliary.a())?;
    let u = master.key.sample_columns(&image)?;
    validate_issued_matrix(public, participant, auxiliary.a(), &u)?;
    if u.rank() != config.m() {
        return Err(PltError::RankDeficient);
    }
    let gram_inverse = (&u.transpose() * &u)
        .inverse()
        .ok_or(PltError::RankDeficient)?;
    // lambda_max(U^T U) <= ||U||_F^2 and ||B B^T||_2 <= ||B||_F^2.
    // Thus Sigma/r^2 - 5 B B^T - I is positive definite with a wide margin.
    let minimum_sigma_squared = Q::from(
        4 * Mp12Key::ROUNDING.pow(2)
            * frobenius_squared(&u)
            * (5 * auxiliary.full_trapdoor_norm_squared() + 1),
    );
    Ok(PartialTrapdoor {
        config: config.clone(),
        participant,
        master_a: public.a.clone(),
        u,
        auxiliary,
        gram_inverse,
        minimum_sigma_squared,
    })
}

/// Samples a partial integer preimage using only public data and one partial trapdoor.
/// `recovery_set` contains exactly t distinct indices and must include this participant.
/// All parties must use the same target, recovery set and `sigma_squared`.
/// Uses Sigma = sigma_squared * (U^T U)^-1 in the existing MP12 sampler.
/// Returns an error for invalid targets/sets, mismatched keys or insufficient width.
pub fn psamp_pre(
    public: &PublicParameters,
    partial: &PartialTrapdoor,
    recovery_set: &[usize],
    target: &MatZq,
    sigma_squared: &Q,
) -> Result<PartialPreimage, PltError> {
    if public.config != partial.config || public.a != partial.master_a {
        return Err(PltError::KeyMismatch);
    }
    let components = decompose_target(&public.config, recovery_set, target)?;
    let position = recovery_set
        .iter()
        .position(|&id| id == partial.participant)
        .ok_or(PltError::InvalidRecoverySet)?;
    validate_width(sigma_squared)?;
    if sigma_squared < &partial.minimum_sigma_squared {
        return Err(PltError::WidthTooSmall);
    }
    let z = components[position].clone();
    let squared_width = sigma_squared * &partial.gram_inverse;
    let d = partial
        .auxiliary
        .sample_with_squared_width(&z, &squared_width)?;
    if partial.c() * &d != z {
        return Err(PltError::AlgebraicMismatch);
    }
    let x = &partial.u * &d;
    verify_partial_preimage(
        public,
        partial.participant,
        recovery_set,
        target,
        &x,
        sigma_squared,
    )?;
    let mut recovery_set = recovery_set.to_vec();
    recovery_set.sort_unstable();
    Ok(PartialPreimage {
        config: public.config.clone(),
        participant: partial.participant,
        recovery_set,
        master_a: public.a.clone(),
        target: target.clone(),
        sigma_squared: sigma_squared.clone(),
        z,
        d,
        x,
    })
}

/// Adds exactly t partial preimages over Z; no interpolation or modular reduction.
/// Checks count, uniqueness, configuration, target, recovery set, width and public key,
/// then verifies A*x = y mod q. This is not a malicious-share verification protocol.
pub fn rec(public: &PublicParameters, shares: &[&PartialPreimage]) -> Result<MatZ, PltError> {
    if shares.len() != public.config.threshold() {
        return Err(PltError::InvalidShares);
    }
    let first = shares[0]; // threshold is positive by configuration validation.
    let mut participants: Vec<_> = shares.iter().map(|share| share.participant).collect();
    participants.sort_unstable();
    if public.config.validate_set(&participants).is_err() || participants != first.recovery_set {
        return Err(PltError::InvalidShares);
    }
    let mut x = MatZ::new(public.config.preimage_rows(), 1);
    for share in shares {
        if share.config != public.config
            || share.master_a != public.a
            || share.target != first.target
            || share.recovery_set != first.recovery_set
            || share.sigma_squared != first.sigma_squared
        {
            return Err(PltError::InvalidShares);
        }
        verify_partial_preimage(
            public,
            share.participant,
            &share.recovery_set,
            &share.target,
            &share.x,
            &share.sigma_squared,
        )?;
        x += &share.x;
    }
    verify_preimage(public, &first.target, &x, &first.sigma_squared)?;
    Ok(x)
}

/// Verifies an external integer full preimage against a caller-selected target and width.
/// This checks algebra and shortness, not identity, authorization or threshold unforgeability.
pub fn verify_preimage(
    public: &PublicParameters,
    target: &MatZq,
    x: &MatZ,
    sigma_squared: &Q,
) -> Result<(), PltError> {
    public.config.validate_target(target)?;
    validate_width(sigma_squared)?;
    let t = Z::from(public.config.threshold() as i64);
    check_norm(&public.config, x, &(&t * &t * sigma_squared))?;
    if &public.a * x != *target {
        return Err(PltError::AlgebraicMismatch);
    }
    Ok(())
}

/// Verifies the public image and norm of an external partial preimage.
/// It does not certify membership in the secret U sublattice or authenticate the sender.
pub fn verify_partial_preimage(
    public: &PublicParameters,
    participant: usize,
    set: &[usize],
    target: &MatZq,
    x: &MatZ,
    sigma_squared: &Q,
) -> Result<(), PltError> {
    let config = &public.config;
    if participant >= config.participants() {
        return Err(PltError::InvalidParticipant);
    }
    let components = decompose_target(config, set, target)?;
    let index = set
        .iter()
        .position(|&id| id == participant)
        .ok_or(PltError::InvalidRecoverySet)?;
    validate_width(sigma_squared)?;
    check_norm(config, x, sigma_squared)?;
    if &public.a * x != participant_image(config, participant, &components[index])? {
        return Err(PltError::AlgebraicMismatch);
    }
    Ok(())
}

fn check_norm(config: &PltConfig, x: &MatZ, squared_width: &Q) -> Result<(), PltError> {
    if x.get_num_rows() != config.preimage_rows() || x.get_num_columns() != 1 {
        return Err(PltError::InvalidPreimage);
    }
    if x.norm_eucl_sqrd().map_err(backend_error)? > squared_width * config.preimage_rows() {
        return Err(PltError::NormExceeded);
    }
    Ok(())
}

fn validate_width(width: &Q) -> Result<(), PltError> {
    if width <= &Q::ZERO {
        return Err(PltError::WidthTooSmall);
    }
    if !f64::from(width).is_finite() {
        return Err(PltError::NumericalFailure);
    }
    Ok(())
}

fn validate_issued_matrix(
    public: &PublicParameters,
    id: usize,
    c: &MatZq,
    u: &MatZ,
) -> Result<(), PltError> {
    let config = &public.config;
    if u.get_num_rows() != config.preimage_rows() || u.get_num_columns() != config.m() {
        return Err(PltError::InvalidPreimage);
    }
    if u.rank() != config.m() {
        return Err(PltError::RankDeficient);
    }
    if &public.a * u != participant_image(config, id, c)? {
        return Err(PltError::AlgebraicMismatch);
    }
    Ok(())
}

#[cfg(test)]
mod acceptance_tests;

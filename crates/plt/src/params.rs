// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Validated dimensions and threshold configuration for the integer PLT algorithms.

use super::PltError;
use qfall_math::{integer::Z, integer_mod_q::MatZq, traits::MatrixDimensions};

/// Immutable configuration shared by Setup, PTrapGen, PSampPre and Rec.
/// Participants are indexed by 0..k with evaluation points 1..=k modulo prime q.
/// Constructors validate 1 <= t <= k < q, positive n and dimension arithmetic.
/// These checks ensure algebraic compatibility, not cryptographic security.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PltConfig {
    n: i64,
    q: i64,
    threshold: usize,
    participants: usize,
    m: i64,
    target_rows: i64,
    preimage_rows: i64,
}

impl PltConfig {
    /// Configures t-of-k with the original toy lattice n=1, q=257.
    /// Returns InvalidParameters for invalid thresholds or too many participants.
    pub fn new(threshold: usize, participants: usize) -> Result<Self, PltError> {
        Self::with_lattice(1, 257, threshold, participants)
    }

    /// Configures n, prime modulus q and t-of-k; derives m=n*(ceil(log2(q))+1).
    /// Rejects invalid parameters and dimensions that overflow i64.
    /// Large valid configurations may still exceed the machine's memory/precision.
    pub fn with_lattice(
        n: i64,
        q: i64,
        threshold: usize,
        participants: usize,
    ) -> Result<Self, PltError> {
        let invalid = PltError::InvalidParameters;
        if n < 1 || q < 2 || threshold == 0 || threshold > participants {
            return Err(invalid);
        }
        let k = i64::try_from(participants).map_err(|_| PltError::InvalidParameters)?;
        let t = i64::try_from(threshold).map_err(|_| PltError::InvalidParameters)?;
        if k >= q || !Z::from(q).is_prime() {
            return Err(invalid);
        }
        let gadget_length = i64::from(i64::BITS - (q - 1).leading_zeros());
        let m = n
            .checked_mul(gadget_length + 1)
            .ok_or(PltError::InvalidParameters)?;
        let target_rows = n.checked_mul(t).ok_or(PltError::InvalidParameters)?;
        let preimage_rows = m
            .checked_mul(t)
            .and_then(|v| v.checked_mul(2))
            .ok_or(PltError::InvalidParameters)?;
        Ok(Self {
            n,
            q,
            threshold,
            participants,
            m,
            target_rows,
            preimage_rows,
        })
    }

    /// Number of rows in an auxiliary public matrix C.
    pub fn n(&self) -> i64 {
        self.n
    }
    /// Prime modulus of all public matrices.
    pub fn modulus(&self) -> i64 {
        self.q
    }
    /// Number of columns in C and U.
    pub fn m(&self) -> i64 {
        self.m
    }
    /// Number of distinct participants required in a recovery set.
    pub fn threshold(&self) -> usize {
        self.threshold
    }
    /// Total number of participants.
    pub fn participants(&self) -> usize {
        self.participants
    }
    /// Number of rows in A and the public target y: n*t.
    pub fn target_rows(&self) -> i64 {
        self.target_rows
    }
    /// Number of rows in U and integer preimages: 2*m*t.
    pub fn preimage_rows(&self) -> i64 {
        self.preimage_rows
    }

    /// Lazily enumerates all t-subsets in lexicographic order, without storing
    /// all binomial(k,t) sets in memory. Intended for tests and small experiments.
    pub fn recovery_sets(&self) -> impl Iterator<Item = Vec<usize>> + use<> {
        let (t, k) = (self.threshold, self.participants);
        std::iter::successors(Some((0..t).collect::<Vec<_>>()), move |current| {
            let mut next = current.clone();
            let index = (0..t).rev().find(|&i| next[i] < k - t + i)?;
            next[index] += 1;
            for i in index + 1..t {
                next[i] = next[i - 1] + 1;
            }
            Some(next)
        })
    }

    /// Validates the cardinality, participant range and uniqueness of a set.
    pub(super) fn validate_set(&self, set: &[usize]) -> Result<(), PltError> {
        if set.len() != self.threshold || set.iter().any(|&id| id >= self.participants) {
            return Err(PltError::InvalidRecoverySet);
        }
        let mut sorted = set.to_vec();
        sorted.sort_unstable();
        if sorted.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(PltError::InvalidRecoverySet);
        }
        Ok(())
    }

    /// Checks the configured target shape and modulus before matrix arithmetic.
    pub(super) fn validate_target(&self, target: &MatZq) -> Result<(), PltError> {
        if target.get_num_rows() != self.target_rows
            || target.get_num_columns() != 1
            || target.get_mod() != self.q
        {
            return Err(PltError::InvalidTarget);
        }
        Ok(())
    }
}

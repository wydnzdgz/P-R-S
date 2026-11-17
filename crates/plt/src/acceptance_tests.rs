// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Fault injection for acceptance A1-A4; private mutations never enter the public API.

use super::*;
use qfall_math::traits::{MatrixGetEntry, MatrixSetEntry, Pow};
use std::str::FromStr;

#[test]
fn issuance_refuses_wrong_shape_rank_and_relation() {
    let config = PltConfig::new(2, 3).unwrap();
    let (public, master) = setup(&config).unwrap();
    let partial = ptrap_gen(&public, &master, 0).unwrap();
    let zero = MatZ::new(config.preimage_rows(), config.m());
    assert_eq!(
        validate_issued_matrix(&public, 0, partial.c(), &zero),
        Err(PltError::RankDeficient)
    );
    let wrong_shape = MatZ::new(config.preimage_rows() - 1, config.m());
    assert_eq!(
        validate_issued_matrix(&public, 0, partial.c(), &wrong_shape),
        Err(PltError::InvalidPreimage)
    );
    let mut wrong = partial.u.clone();
    // Negation preserves full rational rank and cannot preserve the nonzero image over q=257.
    wrong = -1 * wrong;
    assert_eq!(
        validate_issued_matrix(&public, 0, partial.c(), &wrong),
        Err(PltError::AlgebraicMismatch)
    );
}

#[test]
fn external_verifiers_reject_shortness_and_format_attacks() {
    let config = PltConfig::new(2, 3).unwrap();
    let (public, master) = setup(&config).unwrap();
    let p0 = ptrap_gen(&public, &master, 0).unwrap();
    let p1 = ptrap_gen(&public, &master, 1).unwrap();
    let width = p0.minimum_sigma_squared().max(p1.minimum_sigma_squared());
    let y = MatZq::from_str("[[12],[19]] mod 257").unwrap();
    let left = psamp_pre(&public, &p0, &[0, 1], &y, width).unwrap();
    let right = psamp_pre(&public, &p1, &[0, 1], &y, width).unwrap();
    let x = rec(&public, &[&left, &right]).unwrap();
    let mut long = x.clone();
    let value: Z = long.get_entry(0, 0).unwrap();
    long.set_entry(
        0,
        0,
        value + Z::from(2).pow(128).unwrap() * config.modulus(),
    )
    .unwrap();
    assert_eq!(public.a() * &long, y); // Algebra is unchanged; norm check must still reject it.
    assert_eq!(
        verify_preimage(&public, &y, &long, width),
        Err(PltError::NormExceeded)
    );
    assert_eq!(
        verify_preimage(&public, &y, &MatZ::new(40, 2), width),
        Err(PltError::InvalidPreimage)
    );
    assert_eq!(
        verify_preimage(&public, &y, &MatZ::new(40, 1), width),
        Err(PltError::AlgebraicMismatch)
    );
    assert_eq!(
        verify_preimage(&public, &MatZq::new(2, 1, 17), &x, width),
        Err(PltError::InvalidTarget)
    );
    assert_eq!(
        verify_preimage(&public, &y, &x, &Q::ZERO),
        Err(PltError::WidthTooSmall)
    );
    let enormous = Q::from(Z::from(2).pow(2048).unwrap());
    assert_eq!(
        verify_preimage(&public, &y, &x, &enormous),
        Err(PltError::NumericalFailure)
    );
    assert!(matches!(
        psamp_pre(&public, &p0, &[0, 1], &y, &enormous),
        Err(PltError::NumericalFailure)
    ));
    assert_eq!(
        verify_partial_preimage(&public, 0, &[0, 1], &y, &long, width),
        Err(PltError::NormExceeded)
    );
    assert_eq!(
        verify_partial_preimage(&public, 3, &[0, 1], &y, left.x(), width),
        Err(PltError::InvalidParticipant)
    );
    assert_eq!(
        verify_partial_preimage(&public, 0, &[1, 2], &y, left.x(), width),
        Err(PltError::InvalidRecoverySet)
    );
    assert_eq!(
        verify_partial_preimage(&public, 0, &[0, 1], &y, &MatZ::new(40, 1), width),
        Err(PltError::AlgebraicMismatch)
    );
}

#[test]
fn reconstruction_checks_norms_and_equations_in_addition_to_metadata() {
    let config = PltConfig::new(2, 3).unwrap();
    let (public, master) = setup(&config).unwrap();
    let p0 = ptrap_gen(&public, &master, 0).unwrap();
    let p1 = ptrap_gen(&public, &master, 1).unwrap();
    let width = p0.minimum_sigma_squared().max(p1.minimum_sigma_squared());
    let y = MatZq::from_str("[[12],[19]] mod 257").unwrap();
    let mut left = psamp_pre(&public, &p0, &[0, 1], &y, width).unwrap();
    let right = psamp_pre(&public, &p1, &[0, 1], &y, width).unwrap();
    let saved = left.x.clone();
    let value: Z = left.x.get_entry(0, 0).unwrap();
    left.x
        .set_entry(
            0,
            0,
            value + Z::from(2).pow(128).unwrap() * config.modulus(),
        )
        .unwrap();
    assert_eq!(
        public.a() * &left.x,
        participant_image(&config, 0, left.z()).unwrap()
    );
    assert_eq!(rec(&public, &[&left, &right]), Err(PltError::NormExceeded));
    left.x = MatZ::new(40, 2);
    assert_eq!(
        rec(&public, &[&left, &right]),
        Err(PltError::InvalidPreimage)
    );
    left.x = MatZ::new(40, 1);
    assert_eq!(
        rec(&public, &[&left, &right]),
        Err(PltError::AlgebraicMismatch)
    );
    left.x = saved;
    assert!(rec(&public, &[&right, &left]).is_ok());
}

#[test]
fn generation_profiles_record_the_actual_selected_instantiation() {
    for (t, k) in [(2, 3), (3, 5), (5, 10)] {
        let config = PltConfig::new(t, k).unwrap();
        let (public, master) = setup(&config).unwrap();
        let profile = public.master_sampling_profile();
        assert_eq!(profile.rows, config.target_rows());
        assert_eq!(profile.columns, config.preimage_rows());
        assert_eq!(profile.gadget_length, 9);
        assert_eq!(
            profile.m_bar,
            profile.columns - profile.rows * profile.gadget_length
        );
        assert_eq!(profile.rounding_parameter, 4);
        assert!(profile.distribution.contains("PlusMinusOneZero"));
        assert_eq!(profile.initial_width, if t == 2 { 256 } else { 512 });
        let partial = ptrap_gen(&public, &master, 0).unwrap();
        let aux = partial.auxiliary_sampling_profile();
        assert_eq!(aux.rows, config.n());
        assert_eq!(aux.columns, config.m());
        assert_eq!(aux.initial_width, 256);
    }
}

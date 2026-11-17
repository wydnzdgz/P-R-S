// This file is part of qFALL-tools, distributed under MPL-2.0.

//! B2 numerical range checks and B3 fail-closed parameter assessment.

use plt::{PltConfig, PltError, ptrap_gen, setup};
use qfall_math::{
    integer::{MatZ, Z},
    rational::{MatQ, Q},
    traits::MatrixSetEntry,
};
use qfall_tools::{
    primitive::psf::PSFPerturbation, sample::g_trapdoor::gadget_parameters::GadgetParameters,
};

#[test]
fn actual_covariance_geometry_and_numerical_reconstruction() {
    for config in [
        PltConfig::new(2, 3).unwrap(),
        PltConfig::new(3, 5).unwrap(),
        PltConfig::new(5, 10).unwrap(),
        PltConfig::with_lattice(2, 17, 2, 3).unwrap(),
    ] {
        for _ in 0..2 {
            let (public, master) = setup(&config).unwrap();
            for id in [0, config.participants() - 1] {
                let partial = ptrap_gen(&public, &master, id).unwrap();
                let diagnostics = partial
                    .sampling_diagnostics(partial.minimum_sigma_squared())
                    .unwrap();
                assert!(diagnostics.pullback_metric_exact);
                assert!(diagnostics.total_squared_width_relative_residual <= 1e-8);
                assert!(diagnostics.perturbation_condition_inf * f64::EPSILON <= 1e-8);
                assert!(diagnostics.gram_condition_inf.is_finite());
                assert!(matches!(
                    partial.sampling_diagnostics(&Q::ONE),
                    Err(PltError::WidthTooSmall)
                ));
            }
        }
    }
}

#[test]
fn numerically_ill_conditioned_positive_matrix_is_rejected() {
    let sampler = PSFPerturbation {
        gp: GadgetParameters::init_default(1, 257),
        r: Q::from(4),
        s: Q::from(64),
    };
    let r = MatZ::new(1, 9);
    // With B=[0;I], normalized Sigma=diag(2,7,...,7) gives a well-conditioned perturbation.
    let mut sigma: MatQ = 7 * MatQ::identity(10, 10);
    sigma.set_entry(0, 0, 2).unwrap();
    assert!(sampler.try_compute_sqrt_sigma_2(&r, &sigma).is_ok());
    // Exact rational positivity remains true; finite-precision amplification is unacceptable.
    sigma
        .set_entry(0, 0, Q::ONE + Q::from((1, 1_000_000_000_000i64)))
        .unwrap();
    assert!(sampler.try_compute_sqrt_sigma_2(&r, &sigma).is_err());
}

#[test]
fn unassessed_parameters_are_never_approved_for_protocol_security() {
    for (t, k) in [(2, 3), (3, 5), (5, 10)] {
        let config = PltConfig::new(t, k).unwrap();
        let (public, master) = setup(&config).unwrap();
        let partial = ptrap_gen(&public, &master, 0).unwrap();
        let assessment = public
            .assess_parameters(partial.minimum_sigma_squared())
            .unwrap();
        let half = Z::from(128);
        assert_eq!(
            assessment.public_linear_solution_norm_bound_squared,
            &half * &half * config.target_rows()
        );
        assert!(assessment.accepts_public_linear_solutions);
        assert!(!assessment.table2_issuance_upper_bound_satisfied);
        assert_eq!(assessment.table2_m_greater_than_participants, k < 10);
        assert!(!assessment.protocol_security_approved);
        assert_eq!(
            public.require_protocol_parameters(partial.minimum_sigma_squared()),
            Err(PltError::UnapprovedParameters)
        );
        // Tightening the norm alone is not evidence for the remaining conditions.
        assert!(
            !public
                .assess_parameters(&Q::ONE)
                .unwrap()
                .protocol_security_approved
        );
    }
}

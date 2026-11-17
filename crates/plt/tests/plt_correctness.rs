// This file is part of qFALL-tools, distributed under MPL-2.0.

//! The same end-to-end checks exercise configurable integer-lattice t-of-k PLT.

use plt::{
    PltConfig, PltError, decompose_target, participant_image, psamp_pre, ptrap_gen, rec, setup,
};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::Q,
    traits::{MatrixDimensions, MatrixGetEntry, MatrixSetEntry},
};
use std::str::FromStr;

/// Checks every recovery set and all four congruences with one implementation.
fn check_all_sets(config: PltConfig, trials: usize) -> usize {
    let (public, master) = setup(&config).unwrap();
    assert_eq!(public.config(), &config);
    assert_eq!(public.a().get_num_rows(), config.target_rows());
    assert_eq!(public.a().get_num_columns(), config.preimage_rows());
    let partials: Vec<_> = (0..config.participants())
        .map(|id| ptrap_gen(&public, &master, id).unwrap())
        .collect();
    drop(master);
    let width = partials
        .iter()
        .map(|p| p.minimum_sigma_squared())
        .max()
        .unwrap();
    for partial in &partials {
        assert_eq!(partial.u().get_num_rows(), config.preimage_rows());
        assert_eq!(partial.u().get_num_columns(), config.m());
        assert_eq!(partial.u().rank(), config.m());
        assert_eq!(partial.c().get_num_rows(), config.n());
        assert_eq!(partial.c().get_num_columns(), config.m());
        // 1. A U_j = v_j tensor C_j (mod q).
        assert_eq!(
            public.a() * partial.u(),
            participant_image(&config, partial.participant(), partial.c()).unwrap()
        );
    }
    let mut targets = Vec::new();
    for trial in 0..trials {
        let mut y = MatZq::sample_uniform(config.target_rows(), 1, config.modulus());
        if trial == 1 {
            y = MatZq::new(config.target_rows(), 1, config.modulus());
        }
        if trial == 2 {
            for row in 0..config.target_rows() {
                y.set_entry(row, 0, config.modulus() - 1).unwrap();
            }
        }
        targets.push(y);
    }
    let mut sets = 0;
    let mut saw_noncanonical_integer = false;
    for set in config.recovery_sets() {
        sets += 1;
        for target in &targets {
            let shares: Vec<_> = set
                .iter()
                .map(|&id| psamp_pre(&public, &partials[id], &set, target, width).unwrap())
                .collect();
            let mut sum = MatZ::new(config.preimage_rows(), 1);
            for share in &shares {
                let partial = &partials[share.participant()];
                // 2. C_j d_j = z_j (mod q).
                assert_eq!(partial.c() * share.d(), *share.z());
                assert_eq!(share.x(), &(partial.u() * share.d()));
                // 3. A x_j = v_j tensor z_j (mod q).
                assert_eq!(
                    public.a() * share.x(),
                    participant_image(&config, share.participant(), share.z()).unwrap()
                );
                assert!(share.x().norm_eucl_sqrd().unwrap() <= width * config.preimage_rows());
                sum += share.x();
            }
            let mut refs: Vec<_> = shares.iter().collect();
            refs.reverse(); // Collection order is independent of the canonical recovery set.
            let x = rec(&public, &refs).unwrap();
            assert_eq!(x, sum);
            // 4. A x = y (mod q).
            assert_eq!(public.a() * &x, *target);
            let t = Z::from(config.threshold() as i64);
            assert!(x.norm_eucl_sqrd().unwrap() <= &t * &t * width * config.preimage_rows());
            for row in 0..config.preimage_rows() {
                let value: Z = x.get_entry(row, 0).unwrap();
                saw_noncanonical_integer |= value < Z::ZERO || value >= config.modulus();
            }
        }
    }
    assert!(
        saw_noncanonical_integer,
        "integer samples must not be reduced modulo q"
    );
    sets
}

/// All 3 pairs and all 10 triples use the same algorithm and verifier.
#[test]
fn all_sets_for_two_of_three_and_three_of_five() {
    for (t, k, expected) in [(2, 3, 3), (3, 5, 10)] {
        assert_eq!(check_all_sets(PltConfig::new(t, k).unwrap(), 8), expected);
    }
}

/// Tests t=1, t=k and all 252 five-of-ten sets, with no special algorithm branches.
#[test]
fn threshold_boundaries_and_five_of_ten() {
    for (t, k, expected) in [(1, 1, 1), (1, 3, 3), (3, 3, 1), (5, 10, 252)] {
        assert_eq!(check_all_sets(PltConfig::new(t, k).unwrap(), 1), expected);
    }
}

/// A different prime modulus and n>1 exercise block decomposition and dynamic m.
#[test]
fn configurable_lattice_uses_n_row_blocks() {
    let config = PltConfig::with_lattice(2, 17, 2, 3).unwrap();
    assert_eq!(config.m(), 12); // ceil(log2(17))=5.
    assert_eq!(check_all_sets(config, 3), 3);
    assert_eq!(
        check_all_sets(PltConfig::with_lattice(1, 2, 1, 1).unwrap(), 3),
        1
    );
}

/// Independent known answers check Vandermonde order and degree beyond two.
#[test]
fn vandermonde_known_answers_and_permutations() {
    for (t, k, set, target, expected) in [
        (2, 3, vec![0, 1], "[[12],[19]] mod 257", vec![5, 7]),
        (
            3,
            5,
            vec![0, 2, 4],
            "[[23],[81],[86]] mod 257",
            vec![5, 7, 11],
        ),
    ] {
        let config = PltConfig::new(t, k).unwrap();
        let y = MatZq::from_str(target).unwrap();
        let z = decompose_target(&config, &set, &y).unwrap();
        for (component, value) in z.iter().zip(&expected) {
            let actual: Z = component.get_entry(0, 0).unwrap();
            assert_eq!(actual, *value);
        }
        let mut reversed = set.clone();
        reversed.reverse();
        let mut reversed_z = decompose_target(&config, &reversed, &y).unwrap();
        reversed_z.reverse();
        assert_eq!(z, reversed_z);
        for mut set in config.recovery_sets() {
            set.rotate_left(1);
            let components = decompose_target(&config, &set, &y).unwrap();
            let mut recovered = MatZq::new(config.target_rows(), 1, config.modulus());
            for (&id, component) in set.iter().zip(components) {
                recovered += participant_image(&config, id, &component).unwrap();
            }
            assert_eq!(recovered, y);
        }
    }
}

/// Rejects malformed configurations before allocation, plus malformed input boundaries.
#[test]
fn rejects_invalid_parameters_and_inputs() {
    for (n, q, t, k) in [
        (1, 257, 0, 3),
        (1, 257, 4, 3),
        (1, 257, 1, 257),
        (0, 257, 2, 3),
        (1, 15, 2, 3),
        (1, 1, 1, 1),
        (i64::MAX, 257, 2, 3),
    ] {
        assert_eq!(
            PltConfig::with_lattice(n, q, t, k),
            Err(PltError::InvalidParameters)
        );
    }
    for (t, k) in [(2, 3), (3, 5)] {
        let config = PltConfig::new(t, k).unwrap();
        let (public, master) = setup(&config).unwrap();
        assert!(matches!(
            ptrap_gen(&public, &master, k),
            Err(PltError::InvalidParticipant)
        ));
        let partial = ptrap_gen(&public, &master, 0).unwrap();
        let y = MatZq::new(config.target_rows(), 1, config.modulus());
        let set: Vec<_> = (0..t).collect();
        let mut duplicate = set.clone();
        duplicate[t - 1] = 0;
        let mut out_of_range = set.clone();
        out_of_range[t - 1] = k;
        for invalid in [
            vec![],
            set[..t - 1].to_vec(),
            (0..=t).collect(),
            duplicate,
            out_of_range,
        ] {
            assert_eq!(
                decompose_target(&config, &invalid, &y),
                Err(PltError::InvalidRecoverySet)
            );
        }
        for invalid in [
            MatZq::new(1, 1, 257),
            MatZq::new(config.target_rows(), 2, 257),
            MatZq::new(config.target_rows(), 1, 17),
        ] {
            assert!(matches!(
                psamp_pre(
                    &public,
                    &partial,
                    &set,
                    &invalid,
                    partial.minimum_sigma_squared()
                ),
                Err(PltError::InvalidTarget)
            ));
        }
        let absent: Vec<_> = (1..=t).collect();
        assert!(matches!(
            psamp_pre(
                &public,
                &partial,
                &absent,
                &y,
                partial.minimum_sigma_squared()
            ),
            Err(PltError::InvalidRecoverySet)
        ));
        for width in [Q::from(-1), Q::ZERO, Q::ONE] {
            assert!(matches!(
                psamp_pre(&public, &partial, &set, &y, &width),
                Err(PltError::WidthTooSmall)
            ));
        }
        let other = PltConfig::new(t, k + 1).unwrap();
        let (other_public, other_master) = setup(&other).unwrap();
        assert!(matches!(
            ptrap_gen(&public, &other_master, 0),
            Err(PltError::KeyMismatch)
        ));
        assert!(matches!(
            psamp_pre(
                &other_public,
                &partial,
                &set,
                &y,
                partial.minimum_sigma_squared()
            ),
            Err(PltError::KeyMismatch)
        ));
    }
}

/// Exactly t distinct, matching shares are required; permutations are accepted.
#[test]
fn reconstruction_checks_cardinality_and_context() {
    for (t, k) in [(2, 3), (3, 5)] {
        let config = PltConfig::new(t, k).unwrap();
        let (public, master) = setup(&config).unwrap();
        let partials: Vec<_> = (0..k)
            .map(|id| ptrap_gen(&public, &master, id).unwrap())
            .collect();
        let width = partials
            .iter()
            .map(|p| p.minimum_sigma_squared())
            .max()
            .unwrap();
        let set: Vec<_> = (0..t).collect();
        let y = MatZq::sample_uniform(config.target_rows(), 1, config.modulus());
        let mut reverse = set.clone();
        reverse.reverse();
        let shares: Vec<_> = set
            .iter()
            .map(|&id| psamp_pre(&public, &partials[id], &reverse, &y, width).unwrap())
            .collect();
        let refs: Vec<_> = shares.iter().collect();
        assert_eq!(public.a() * rec(&public, &refs).unwrap(), y);
        assert_eq!(rec(&public, &[]), Err(PltError::InvalidShares));
        assert_eq!(rec(&public, &refs[..t - 1]), Err(PltError::InvalidShares));
        let mut too_many = refs.clone();
        too_many.push(refs[0]);
        assert_eq!(rec(&public, &too_many), Err(PltError::InvalidShares));
        let mut duplicate = refs.clone();
        duplicate[t - 1] = refs[0];
        assert_eq!(rec(&public, &duplicate), Err(PltError::InvalidShares));
        let mut other_target = y.clone();
        let value: Z = other_target.get_entry(0, 0).unwrap();
        other_target.set_entry(0, 0, value + 1).unwrap();
        let alternate_set: Vec<_> = (1..=t).collect();
        let mixed_target =
            psamp_pre(&public, &partials[t - 1], &set, &other_target, width).unwrap();
        let mixed_set = psamp_pre(&public, &partials[t - 1], &alternate_set, &y, width).unwrap();
        let mixed_width = psamp_pre(&public, &partials[t - 1], &set, &y, &(width * 2)).unwrap();
        for other in [&mixed_target, &mixed_set, &mixed_width] {
            let mut mixed = refs.clone();
            mixed[t - 1] = other;
            assert_eq!(rec(&public, &mixed), Err(PltError::InvalidShares));
        }
        let (other_public, _) = setup(&config).unwrap();
        assert_eq!(rec(&other_public, &refs), Err(PltError::InvalidShares));
    }
}

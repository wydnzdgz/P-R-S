// This file is part of qFALL-tools, distributed under MPL-2.0.

//! One experiment driver for configurable integer-lattice t-of-k PLT.
//! Run `cargo run --locked -p plt --example plt -- --threshold 3 --participants 5 --output target/plt_3_of_5.json`.

use plt::{PltConfig, participant_image, psamp_pre, ptrap_gen, rec, setup};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    traits::MatrixSetEntry,
};
use serde_json::json;
use std::error::Error;

/// Shared CLI used by this example and the legacy plt_2_of_3 entry point.
pub fn main() -> Result<(), Box<dyn Error>> {
    let (mut threshold, mut participants, mut n, mut q, mut trials) = (2, 3, 1, 257, 8usize);
    let mut output_path = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--threshold" => threshold = args.next().ok_or("missing threshold")?.parse()?,
            "--participants" => {
                participants = args.next().ok_or("missing participant count")?.parse()?
            }
            "--n" => n = args.next().ok_or("missing n")?.parse()?,
            "--modulus" => q = args.next().ok_or("missing modulus")?.parse()?,
            "--trials" => trials = args.next().ok_or("missing trial count")?.parse()?,
            "--output" => output_path = Some(args.next().ok_or("missing output path")?),
            "--help" => {
                println!(
                    "plt [--threshold T] [--participants K] [--n N] [--modulus Q] [--trials COUNT] [--output FILE]"
                );
                return Ok(());
            }
            path if !path.starts_with('-') && output_path.is_none() => {
                output_path = Some(path.to_owned())
            }
            _ => return Err(format!("unrecognized argument: {arg}").into()),
        }
    }
    if trials == 0 {
        return Err("trials must be positive".into());
    }
    let config = PltConfig::with_lattice(n, q, threshold, participants)?;
    let (public, master) = setup(&config)?;
    let partials: Vec<_> = (0..config.participants())
        .map(|id| ptrap_gen(&public, &master, id))
        .collect::<Result<_, _>>()?;
    drop(master);
    let sigma_squared = partials
        .iter()
        .map(|p| p.minimum_sigma_squared())
        .max()
        .unwrap();
    let t = Z::from(config.threshold() as i64);
    let beta_squared = &t * &t * sigma_squared * config.preimage_rows();
    println!(
        "Integer PLT toy experiment: q={q}, n={n}, m={}, t={threshold}, k={participants}",
        config.m()
    );
    println!(
        "A: {}x{}; C_j: {}x{}; U_j: {}x{}",
        config.target_rows(),
        config.preimage_rows(),
        n,
        config.m(),
        config.preimage_rows(),
        config.m()
    );
    println!("Common sigma^2={sigma_squared}; beta^2={beta_squared}");
    let mut partial_records = Vec::new();
    for partial in &partials {
        assert_eq!(
            public.a() * partial.u(),
            participant_image(&config, partial.participant(), partial.c())?
        );
        partial_records.push(json!({
            "auxiliary_sampling_profile": partial.auxiliary_sampling_profile(),
            "participant": partial.participant(), "point": partial.participant()+1,
            "u_integer": partial.u(), "c_mod_q": partial.c(), "au_equals_v_tensor_c": true,
        }));
    }
    let targets: Vec<_> = (0..trials)
        .map(|trial| {
            let mut target = MatZq::new(config.target_rows(), 1, q);
            if trial == 0 {
                for row in 0..config.target_rows() {
                    target.set_entry(row, 0, 12 + 7 * row).unwrap();
                }
            } else if trial == 2 {
                for row in 0..config.target_rows() {
                    target.set_entry(row, 0, q - 1).unwrap();
                }
            } else if trial > 2 {
                target = MatZq::sample_uniform(config.target_rows(), 1, q);
            }
            target
        })
        .collect();
    let mut records = Vec::new();
    let mut set_count = 0usize;
    for set in config.recovery_sets() {
        set_count += 1;
        for (trial, target) in targets.iter().enumerate() {
            let shares = set
                .iter()
                .map(|&id| psamp_pre(&public, &partials[id], &set, target, sigma_squared))
                .collect::<Result<Vec<_>, _>>()?;
            let mut share_records = Vec::new();
            let mut sum = MatZ::new(config.preimage_rows(), 1);
            for share in &shares {
                let partial = &partials[share.participant()];
                assert_eq!(partial.c() * share.d(), *share.z());
                assert_eq!(partial.u() * share.d(), *share.x());
                assert_eq!(
                    public.a() * share.x(),
                    participant_image(&config, share.participant(), share.z())?
                );
                assert!(share.x().norm_eucl_sqrd()? <= sigma_squared * config.preimage_rows());
                sum += share.x();
                share_records.push(json!({
                    "participant": share.participant(), "z_mod_q": share.z(), "d_integer": share.d(), "x_integer": share.x(),
                    "norm_squared": share.x().norm_eucl_sqrd()?, "cd_equals_z": true, "ax_equals_v_tensor_z": true,
                }));
            }
            let x = rec(&public, &shares.iter().collect::<Vec<_>>())?;
            assert_eq!(public.a() * &x, *target);
            assert_eq!(x, sum);
            let norm_squared = x.norm_eucl_sqrd()?;
            assert!(norm_squared <= beta_squared);
            println!("T={set:?}, trial={trial}: all identities OK, ||x||^2={norm_squared}");
            records.push(json!({
                "recovery_set": set, "trial": trial, "target_mod_q": target, "shares": share_records,
                "x_integer": x, "norm_squared": norm_squared, "ax_equals_y": true, "within_beta": true,
            }));
        }
    }
    println!(
        "Passed {} reconstructions across all {set_count} recovery sets.",
        records.len()
    );
    if let Some(path) = output_path {
        let report = json!({
            "master_sampling_profile": public.master_sampling_profile(),
            "scope": "toy algebraic reproduction; sampling distribution and security not certified",
            "parameters": {"q": q, "n": n, "m": config.m(), "t": threshold, "k": participants},
            "trials_per_set": trials,
            "gaussian_convention": "density proportional to exp(-pi*x^T*Sigma^-1*x)",
            "trapdoor_distribution": "existing MP12 ternary PlusMinusOneZero",
            "sigma_squared": sigma_squared, "beta_squared": beta_squared, "a_mod_q": public.a(),
            "partial_trapdoors": partial_records, "reconstructions": records,
        });
        std::fs::write(&path, serde_json::to_string_pretty(&report)?)?;
        println!("Saved integer samples and checks to {path}");
    }
    Ok(())
}

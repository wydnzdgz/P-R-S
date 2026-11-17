// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Repeated A1-A4 acceptance experiment, with every abort counted and no sample retries.
//! `cargo run --locked -p plt --release --example plt_acceptance -- --rounds 1000 --keysets 3 --output target/plt_acceptance.json`

use plt::{PltConfig, PltError, psamp_pre, ptrap_gen, rec, setup};
use qfall_math::{integer_mod_q::MatZq, traits::MatrixSetEntry};
use serde::Serialize;
use serde_json::json;
use std::{collections::BTreeMap, error::Error};

#[derive(Default, Serialize)]
struct Counts {
    setup_attempts: usize,
    issuance_attempts: usize,
    partial_sampling_attempts: usize,
    reconstruction_attempts: usize,
    reconstruction_successes: usize,
    reconstruction_failures: usize,
    skipped_reconstructions: usize,
    norm_failures: usize,
    errors_by_phase_and_reason: BTreeMap<String, usize>,
}
impl Counts {
    fn failure(&mut self, phase: &str, error: &PltError) {
        if *error == PltError::NormExceeded {
            self.norm_failures += 1;
        }
        *self
            .errors_by_phase_and_reason
            .entry(format!("{phase}: {error}"))
            .or_default() += 1;
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let (mut rounds, mut keysets) = (1000usize, 3usize);
    let mut output = "target/plt_acceptance.json".to_string();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--rounds" => rounds = args.next().ok_or("missing rounds")?.parse()?,
            "--keysets" => keysets = args.next().ok_or("missing keysets")?.parse()?,
            "--output" => output = args.next().ok_or("missing output path")?,
            _ => return Err(format!("unknown argument {arg}").into()),
        }
    }
    if rounds == 0 || keysets == 0 {
        return Err("rounds and keysets must be positive".into());
    }
    let mut reports = Vec::new();
    let mut all_passed = true;
    for (t, k) in [(2, 3), (3, 5), (5, 10)] {
        let config = PltConfig::new(t, k)?;
        let sets: Vec<_> = config.recovery_sets().collect();
        let mut counts = Counts::default();
        let mut profiles = Vec::new();
        let mut set_attempts = BTreeMap::new();
        for keyset in 0..keysets {
            counts.setup_attempts += 1;
            let (public, master) = match setup(&config) {
                Ok(value) => value,
                Err(error) => {
                    counts.failure("Setup", &error);
                    counts.skipped_reconstructions += rounds;
                    continue;
                }
            };
            let mut partials = Vec::new();
            for id in 0..k {
                counts.issuance_attempts += 1;
                match ptrap_gen(&public, &master, id) {
                    Ok(partial) => partials.push(partial),
                    Err(error) => {
                        counts.failure("PTrapGen", &error);
                        break;
                    }
                }
            }
            drop(master);
            if partials.len() != k {
                counts.skipped_reconstructions += rounds;
                continue;
            }
            let width = partials
                .iter()
                .map(|p| p.minimum_sigma_squared())
                .max()
                .unwrap();
            profiles.push(
                json!({"keyset":keyset,"master":public.master_sampling_profile(),
                "auxiliary":partials[0].auxiliary_sampling_profile(),"sigma_squared":width}),
            );
            for trial in 0..rounds {
                let set = &sets[trial % sets.len()];
                *set_attempts.entry(format!("{set:?}")).or_insert(0usize) += 1;
                counts.reconstruction_attempts += 1;
                let mut target = MatZq::sample_uniform(config.target_rows(), 1, config.modulus());
                if trial % 3 == 0 {
                    target = MatZq::new(config.target_rows(), 1, config.modulus());
                }
                if trial % 3 == 1 {
                    for row in 0..config.target_rows() {
                        target.set_entry(row, 0, config.modulus() - 1)?;
                    }
                }
                let mut shares = Vec::new();
                for &id in set {
                    counts.partial_sampling_attempts += 1;
                    match psamp_pre(&public, &partials[id], set, &target, width) {
                        Ok(share) => shares.push(share),
                        Err(error) => {
                            counts.failure("PSampPre", &error);
                            break;
                        }
                    }
                }
                if shares.len() != t {
                    counts.reconstruction_failures += 1;
                    continue;
                }
                match rec(&public, &shares.iter().collect::<Vec<_>>()) {
                    Ok(x) if public.a() * &x == target => counts.reconstruction_successes += 1,
                    Ok(_) => {
                        counts.failure("Rec", &PltError::AlgebraicMismatch);
                        counts.reconstruction_failures += 1;
                    }
                    Err(error) => {
                        counts.failure("Rec", &error);
                        counts.reconstruction_failures += 1;
                    }
                }
            }
        }
        let passed = counts.errors_by_phase_and_reason.is_empty()
            && counts.skipped_reconstructions == 0
            && counts.reconstruction_successes == rounds * keysets;
        all_passed &= passed;
        println!(
            "{t}-of-{k}: {}/{} reconstructions passed, {} partial attempts, {} norm failures, {} skipped",
            counts.reconstruction_successes,
            counts.reconstruction_attempts,
            counts.partial_sampling_attempts,
            counts.norm_failures,
            counts.skipped_reconstructions
        );
        reports.push(json!({"parameters":{"n":config.n(),"q":config.modulus(),"m":config.m(),"t":t,"k":k},
            "counts":counts,"key_profiles":profiles,"recovery_set_attempts":set_attempts,"passed":passed}));
    }
    std::fs::write(
        &output,
        serde_json::to_string_pretty(&json!({
            "scope":"A1-A4 engineering acceptance for the selected ternary MP12 instance; not sampling-distribution or security certification",
            "rounds_per_keyset":rounds,"independent_keysets_per_configuration":keysets,
            "failure_policy":"abort the failed operation; count failures and skipped trials; no internal retries",
            "reports":reports,"passed":all_passed,
        }))?,
    )?;
    println!("Saved acceptance counts to {output}");
    if !all_passed {
        return Err("acceptance experiment recorded failures; inspect report".into());
    }
    Ok(())
}

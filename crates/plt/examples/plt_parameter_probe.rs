// This file is part of qFALL-tools, distributed under MPL-2.0.

//! One bounded parameter experiment. Use scripts/sweep_plt_parameters.py for
//! process isolation, time and memory budgets. No security bit target is assumed.

#[path = "support/public_solution.rs"]
mod public_solver;

use plt::{PltConfig, psamp_pre, ptrap_gen, rec, setup, verify_preimage};
use public_solver::{public_solution, public_solution_in_order};
use qfall_math::{integer::Z, integer_mod_q::MatZq, traits::MatrixDimensions};
use rand::seq::SliceRandom;
use serde_json::{Value, json};
use std::{error::Error, path::Path, time::Instant};

fn checkpoint(path: &Path, report: &Value) -> Result<(), Box<dyn Error>> {
    let temporary = path.with_extension("tmp");
    std::fs::write(&temporary, serde_json::to_string_pretty(report)?)?;
    std::fs::rename(temporary, path)?;
    Ok(())
}

fn phase(path: &Path, report: &mut Value, name: &str) -> Result<(), Box<dyn Error>> {
    report["phase"] = json!(name);
    checkpoint(path, report)
}

fn milliseconds(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1000.0
}

fn timing(values: &[f64]) -> Value {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    json!({"count":values.len(), "values_ms":values,
        "mean_ms":values.iter().sum::<f64>() / values.len() as f64,
        "min_ms":sorted[0], "max_ms":sorted[sorted.len()-1]})
}

fn peak_rss_kib() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()?
        .lines()
        .find_map(|line| {
            line.strip_prefix("VmHWM:")?
                .split_whitespace()
                .next()?
                .parse()
                .ok()
        })
}

fn probe(
    config: &PltConfig,
    trials: usize,
    attempts: usize,
    path: &Path,
    report: &mut Value,
) -> Result<(), Box<dyn Error>> {
    phase(path, report, "Setup")?;
    let start = Instant::now();
    let (public, master) = setup(config)?;
    report["setup_ms"] = json!(milliseconds(start));
    report["master_sampling_profile"] = json!(public.master_sampling_profile());
    let mut partials = Vec::new();
    let mut issuance_times = Vec::new();
    for id in 0..config.participants() {
        phase(path, report, &format!("PTrapGen participant {id}"))?;
        let start = Instant::now();
        partials.push(ptrap_gen(&public, &master, id)?);
        issuance_times.push(milliseconds(start));
        report["ptrap_gen"] = timing(&issuance_times);
    }
    drop(master);
    let sigma_squared = partials
        .iter()
        .map(|p| p.minimum_sigma_squared())
        .max()
        .unwrap()
        .clone();
    let assessment = public.assess_parameters(&sigma_squared)?;
    report["sigma_squared"] = json!(sigma_squared);
    report["sigma_approx"] = json!(f64::from(&sigma_squared).sqrt());
    report["beta_approx"] = json!(f64::from(&assessment.accepted_full_norm_bound_squared).sqrt());
    report["assessment"] = json!(assessment);
    report["auxiliary_sampling_profile"] = json!(partials[0].auxiliary_sampling_profile());
    // This probe measures one recovery set; exhaustive A tests are separate.
    let set: Vec<_> = (0..config.threshold()).collect();
    report["recovery_set"] = json!(set);
    let mut share_times = Vec::new();
    let mut rec_times = Vec::new();
    let mut group_times = Vec::new();
    let mut norms = Vec::new();
    let mut outputs = Vec::new();
    for trial in 0..trials {
        phase(path, report, &format!("PSampPre/Rec trial {trial}"))?;
        let target = MatZq::sample_uniform(config.target_rows(), 1, config.modulus());
        let group_start = Instant::now();
        let mut shares = Vec::new();
        for &id in &set {
            let start = Instant::now();
            shares.push(psamp_pre(
                &public,
                &partials[id],
                &set,
                &target,
                &sigma_squared,
            )?);
            share_times.push(milliseconds(start));
        }
        let start = Instant::now();
        let x = rec(&public, &shares.iter().collect::<Vec<_>>())?;
        rec_times.push(milliseconds(start));
        group_times.push(milliseconds(group_start));
        verify_preimage(&public, &target, &x, &sigma_squared)?;
        if public.a() * &x != target {
            return Err("independent final congruence check failed".into());
        }
        norms.push(x.norm_eucl_sqrd()?.to_string());
        outputs.push(json!({"target_mod_q":target,"x_integer":x}));
        report["completed_trials"] = json!(trial + 1);
        report["psamp_pre_per_share"] = timing(&share_times);
        report["rec"] = timing(&rec_times);
        report["online_group_sequential"] = timing(&group_times);
        report["integer_output_norms_squared"] = json!(norms);
        report["online_outputs"] = json!(outputs);
        checkpoint(path, report)?;
    }
    report["algorithm_flow_completed"] = json!(true);
    report["a_mod_q"] = json!(public.a());
    phase(path, report, "sampling diagnostics (first participant)")?;
    let start = Instant::now();
    report["sampling_diagnostics"] = json!(partials[0].sampling_diagnostics(&sigma_squared)?);
    report["sampling_diagnostics_ms"] = json!(milliseconds(start));
    drop(partials);
    phase(
        path,
        report,
        "public sparse-solution attempts (all trapdoors dropped)",
    )?;
    let target = MatZq::sample_uniform(config.target_rows(), 1, config.modulus());
    let mut order: Vec<_> = (0..public.a().get_num_columns() as usize).collect();
    let (mut accepted, mut minimum, mut witness) = (0usize, None::<Z>, None);
    let mut solver_errors = Vec::new();
    let start = Instant::now();
    for attempt in 0..attempts {
        let result = if attempt == 0 {
            public_solution(public.a(), &target)
        } else {
            order.shuffle(&mut rand::rng());
            public_solution_in_order(public.a(), &target, &order)
        };
        let x = match result {
            Ok(x) => x,
            Err(error) => {
                solver_errors.push(error.to_string());
                continue;
            }
        };
        if public.a() * &x != target {
            return Err("public solver returned an invalid congruence".into());
        }
        let norm = x.norm_eucl_sqrd()?;
        if minimum.as_ref().is_none_or(|current| &norm < current) {
            minimum = Some(norm.clone());
            witness = Some(x.clone());
        }
        if verify_preimage(&public, &target, &x, &sigma_squared).is_ok() {
            accepted += 1;
        }
    }
    report["public_attack"] = json!({"attempts":attempts, "accepted":accepted,
        "solver_errors":solver_errors, "elapsed_ms":milliseconds(start),
        "minimum_norm_squared":minimum, "minimum_norm_witness":witness,
        "a_mod_q":public.a(), "target_mod_q":target,
        "secrets_dropped_before_attack":true,
        "interpretation":"zero accepted attempts does not establish hardness"});
    // This gate remains closed even when a finite attack screen finds no witness.
    report["protocol_gate_rejects"] =
        json!(public.require_protocol_parameters(&sigma_squared).is_err());
    report["status"] = json!("completed");
    phase(path, report, "done")
}

fn main() -> Result<(), Box<dyn Error>> {
    let (mut n, mut q, mut t, mut k, mut trials, mut attempts) = (1, 257, 2, 3, 3, 64);
    let mut output = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let value = args.next().ok_or("missing argument value")?;
        match arg.as_str() {
            "--n" => n = value.parse()?,
            "--modulus" => q = value.parse()?,
            "--threshold" => t = value.parse()?,
            "--participants" => k = value.parse()?,
            "--trials" => trials = value.parse()?,
            "--public-attempts" => attempts = value.parse()?,
            "--output" => output = Some(value),
            _ => return Err(format!("unknown argument: {arg}").into()),
        }
    }
    if trials == 0 || attempts == 0 {
        return Err("trials and public attempts must be positive".into());
    }
    let config = PltConfig::with_lattice(n, q, t, k)?;
    let output = output.ok_or("--output is required")?;
    let path = Path::new(&output);
    let start = Instant::now();
    let mut report = json!({"status":"running", "security_bit_target":null,
        "parameters":{"n":n,"q":q,"t":t,"k":k,"m":config.m(),
            "target_rows":config.target_rows(),"preimage_rows":config.preimage_rows()},
        "requested_trials":trials,"requested_public_attempts":attempts,
        "completed_trials":0,"scope":"single fresh key, one recovery set, engineering resource probe"});
    let outcome = probe(&config, trials, attempts, path, &mut report);
    report["elapsed_ms"] = json!(milliseconds(start));
    report["process_peak_rss_kib"] = json!(peak_rss_kib());
    if let Err(ref error) = outcome {
        report["status"] = json!("algorithm_error");
        report["error"] = json!(error.to_string());
    }
    checkpoint(path, &report)?;
    println!(
        "{}: {} ({:.1} ms)",
        output,
        report["status"],
        milliseconds(start)
    );
    outcome
}

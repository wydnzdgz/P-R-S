// This file is part of qFALL-tools, distributed under MPL-2.0.

//! B1-B4 audit: numerical evidence, public-solution counterexample and sampling screens.
//! `cargo run --locked -p plt --release --example plt_sampling_audit -- target/plt_sampling_B.json`

#[path = "support/public_solution.rs"]
mod public_solver;
#[path = "support/sampling_reference.rs"]
mod sampling_reference;
use plt::{PltConfig, psamp_pre, ptrap_gen, setup, verify_preimage};
use public_solver::public_solution;
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::Q,
    traits::{MatrixGetEntry, MatrixSetEntry},
};
use serde_json::{Value, json};
use std::{collections::HashSet, error::Error};

fn as_float(matrix: &MatZ, row: usize) -> f64 {
    let value: Z = matrix.get_entry(row, 0).unwrap();
    f64::from(&Q::from(value))
}

fn summarize(samples: &[Vec<f64>], zero_target: bool, duplicates: usize) -> Value {
    let n = samples.len();
    let dim = samples[0].len();
    let means: Vec<_> = (0..dim)
        .map(|i| samples.iter().map(|s| s[i]).sum::<f64>() / n as f64)
        .collect();
    let variances: Vec<_> = (0..dim)
        .map(|i| {
            samples
                .iter()
                .map(|s| (s[i] - means[i]).powi(2))
                .sum::<f64>()
                / n as f64
        })
        .collect();
    let mut max_cross = 0.0f64;
    let mut max_lag = 0.0f64;
    for i in 0..dim {
        let lag = samples
            .windows(2)
            .map(|s| (s[0][i] - means[i]) * (s[1][i] - means[i]))
            .sum::<f64>()
            / n as f64
            / variances[i];
        max_lag = max_lag.max(lag.abs());
        for j in 0..i {
            let covariance = samples
                .iter()
                .map(|s| (s[i] - means[i]) * (s[j] - means[j]))
                .sum::<f64>()
                / n as f64;
            max_cross = max_cross.max((covariance / (variances[i] * variances[j]).sqrt()).abs());
        }
    }
    let mean_max = means.iter().map(|x| x.abs()).fold(0.0, f64::max);
    let passed = duplicates == 0
        && means.iter().all(|x| x.is_finite())
        && variances.iter().all(|&v| (0.7..=1.3).contains(&v))
        && max_cross <= 0.2
        && max_lag <= 0.2
        && (!zero_target || mean_max <= 0.15);
    json!({"samples":n,"whitened_means":means,"whitened_variances":variances,"max_abs_coordinate_correlation":max_cross,
        "max_abs_lag1_correlation":max_lag,"duplicate_complete_vectors":duplicates,"zero_target_mean_check_applied":zero_target,"screening_passed":passed})
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let output = args
        .next()
        .unwrap_or_else(|| "target/plt_sampling_B.json".into());
    if args.next().is_some() {
        return Err("usage: plt_sampling_audit [output.json]".into());
    }
    let mut diagnostics = Vec::new();
    let mut cohorts = Vec::new();
    let mut counterexamples = Vec::new();
    let mut errors = Vec::new();
    for (t, k) in [(2, 3), (3, 5), (5, 10)] {
        for keyset in 0..3 {
            let config = PltConfig::new(t, k)?;
            let (public, master) = setup(&config)?;
            let partials: Vec<_> = (0..k)
                .map(|id| ptrap_gen(&public, &master, id))
                .collect::<Result<_, _>>()?;
            drop(master);
            let width = partials
                .iter()
                .map(|p| p.minimum_sigma_squared())
                .max()
                .unwrap()
                .clone();
            for partial in &partials {
                diagnostics.push(json!({"t":t,"k":k,"keyset":keyset,"participant":partial.participant(),"metrics":partial.sampling_diagnostics(&width)?}));
            }
            let selected = keyset * (k - 1) / 2;
            let mut set: Vec<_> = (0..k).filter(|&id| id != selected).take(t - 1).collect();
            set.push(selected);
            set.sort_unstable();
            let partial = &partials[selected];
            let gram = partial.u().transpose() * partial.u();
            let l = qfall_math::rational::MatQ::from(&gram).cholesky_decomposition_flint();
            let dim = config.m() as usize;
            let mut gram_l = vec![vec![0.0; dim]; dim];
            for (i, row) in gram_l.iter_mut().enumerate() {
                for (j, value) in row.iter_mut().enumerate() {
                    let entry: Q = l.get_entry(i, j)?;
                    *value = f64::from(&entry);
                }
            }
            let scale = (2.0 * std::f64::consts::PI / f64::from(&width)).sqrt();
            for category in ["zero", "fixed_nonzero", "fixed_random"] {
                let mut target = MatZq::new(config.target_rows(), 1, config.modulus());
                if category == "fixed_nonzero" {
                    for row in 0..config.target_rows() {
                        target.set_entry(row, 0, 12 + 7 * row)?;
                    }
                }
                if category == "fixed_random" {
                    target = MatZq::sample_uniform(config.target_rows(), 1, config.modulus());
                }
                let mut samples = Vec::new();
                let mut seen = HashSet::new();
                let mut duplicates = 0;
                for trial in 0..1000 {
                    match psamp_pre(&public,partial,&set,&target,&width) {
                    Ok(share)=>{
                        if !seen.insert(share.x().to_string()){duplicates+=1;}
                        let d:Vec<_>=(0..dim).map(|i|as_float(share.d(),i)).collect();
                        let w:Vec<_>=(0..dim).map(|i|scale*(0..dim).map(|j|gram_l[j][i]*d[j]).sum::<f64>()).collect();
                        if w.iter().any(|value|!value.is_finite()){return Err("non-finite whitened sample".into());}
                        samples.push(w);
                    }
                    Err(error)=>errors.push(json!({"t":t,"k":k,"keyset":keyset,"category":category,"trial":trial,"error":error.to_string()})),
                }
                }
                if samples.is_empty() {
                    return Err("all samples failed; inspect sampler".into());
                }
                let summary = summarize(&samples, category == "zero", duplicates);
                println!(
                    "PLT {t}-of-{k} keyset={keyset} target={category}: screen={}",
                    summary["screening_passed"]
                );
                cohorts.push(json!({"t":t,"k":k,"keyset":keyset,"participant":selected,"recovery_set":set,"category":category,"metrics":summary}));
            }
            if keyset == 0 {
                let mut target = MatZq::new(config.target_rows(), 1, config.modulus());
                for row in 0..config.target_rows() {
                    target.set_entry(row, 0, 12 + 7 * row)?;
                }
                let assessment = public.assess_parameters(&width)?;
                let gate_rejected = public.require_protocol_parameters(&width).is_err();
                drop(partials); // The solver below receives public data only.
                let x = public_solution(public.a(), &target)?;
                let accepted = verify_preimage(&public, &target, &x, &width).is_ok();
                println!(
                    "B3 {t}-of-{k}: public solution accepted={accepted}; parameter gate rejected={gate_rejected}"
                );
                counterexamples.push(json!({"t":t,"k":k,"a_mod_q":public.a(),"target_mod_q":target,"public_integer_solution":x,
                "public_solution_accepted_by_preimage_checks":accepted,"parameter_gate_rejected":gate_rejected,"assessment":assessment}));
            }
        }
    }
    let reference = sampling_reference::run_reference()?;
    let engineering = errors.is_empty()
        && cohorts
            .iter()
            .all(|cohort| cohort["metrics"]["screening_passed"] == true)
        && reference["screening_passed"] == true;
    std::fs::write(
        &output,
        serde_json::to_string_pretty(&json!({
            "scope":"B1/B2 engineering audit and B4 finite screening; B3 actual protocol parameters rejected",
            "plan":"docs/plt_B_experiment_plan.md","sampling_attempts":27000,"errors":errors,"numerical_diagnostics":diagnostics,
            "plt_cohorts":cohorts,"small_discrete_reference":reference,"b3_counterexamples":counterexamples,
            "engineering_screening_passed":engineering,"protocol_parameters_passed":false,"stage_b_fully_passed":false,
        }))?,
    )?;
    println!(
        "Saved audit to {output}; engineering_screening_passed={engineering}; B3 remains rejected"
    );
    if !engineering {
        return Err("sampling screens reported anomalies; preserve and inspect results".into());
    }
    Ok(())
}

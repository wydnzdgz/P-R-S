// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Defensive public-only LLL + nearest-plane screen on a small column subset.
//! Reads a saved runtime report; never generates or receives a PLT trapdoor.

#[path = "support/public_solution.rs"]
mod public_solver;

use public_solver::public_solution;
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::{MatQ, Q},
    traits::{MatrixDimensions, MatrixGetEntry, MatrixSetEntry},
};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{error::Error, time::Instant};

fn sha256(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

fn public_nearest_plane(a: &MatZq, y: &MatZq) -> Result<(MatZ, usize), Box<dyn Error>> {
    let rows = a.get_num_rows() as usize;
    // Keep this diagnostic small. It is not a full lattice attack estimator.
    if rows > 16 {
        return Err("screen limited to at most 16 equations".into());
    }
    let columns = a.get_num_columns() as usize;
    let support = columns.min(2 * rows + 8);
    let q = Z::from(a.get_mod());
    let mut leading = MatZq::new(rows, rows, &q);
    for i in 0..rows {
        for j in 0..rows {
            let entry: Z = a.get_entry(i, j)?;
            leading.set_entry(i, j, entry)?;
        }
    }
    // Columns q*e_i and (-B^-1*A_j,e_j) form a full-rank basis of
    // the modular kernel restricted to this subset. All entries are public.
    let mut basis = MatZ::new(support, support);
    for i in 0..rows {
        basis.set_entry(i, i, &q)?;
    }
    for j in rows..support {
        let mut column = MatZq::new(rows, 1, &q);
        for i in 0..rows {
            let entry: Z = a.get_entry(i, j)?;
            column.set_entry(i, 0, entry)?;
        }
        let inverse_column = public_solution(&leading, &column)?;
        for i in 0..rows {
            let entry: Z = inverse_column.get_entry(i, 0)?;
            basis.set_entry(i, j, -entry)?;
        }
        basis.set_entry(j, j, 1)?;
    }
    let initial = public_solution(&leading, y)?;
    let mut residual = MatZ::new(support, 1);
    for i in 0..rows {
        let entry: Z = initial.get_entry(i, 0)?;
        residual.set_entry(i, 0, entry)?;
    }
    let reduced = basis.lll(Q::from((99, 100)), Q::from((501, 1000)));
    let gso = MatQ::from(&reduced).gso();
    for j in (0..support).rev() {
        let (mut numerator, mut denominator) = (Q::ZERO, Q::ZERO);
        for i in 0..support {
            let star: Q = gso.get_entry(i, j)?;
            let entry: Z = residual.get_entry(i, 0)?;
            numerator += Q::from(entry) * &star;
            denominator += &star * &star;
        }
        let coefficient = (numerator / denominator).round();
        for i in 0..support {
            let entry: Z = residual.get_entry(i, 0)?;
            let b: Z = reduced.get_entry(i, j)?;
            residual.set_entry(i, 0, entry - &coefficient * b)?;
        }
    }
    let mut x = MatZ::new(columns, 1);
    for i in 0..support {
        let entry: Z = residual.get_entry(i, 0)?;
        x.set_entry(i, 0, entry)?;
    }
    if a * &x != *y {
        return Err("lattice screen congruence failed".into());
    }
    Ok((x, support))
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let input = args.next().ok_or("expected INPUT_REPORT OUTPUT_REPORT")?;
    let output = args.next().ok_or("expected OUTPUT_REPORT")?;
    let input_bytes = std::fs::read(&input)?;
    let report: Value = serde_json::from_slice(&input_bytes)?;
    let mut results = Vec::new();
    for item in report["results"].as_array().ok_or("missing results")? {
        if item["status"] != "completed" && item["probe"]["algorithm_flow_completed"] != true {
            continue;
        }
        let p = &item["probe"];
        let (a_value, y_value, target_source) = if p["public_attack"].is_object() {
            (
                &p["public_attack"]["a_mod_q"],
                &p["public_attack"]["target_mod_q"],
                "public_attack",
            )
        } else {
            (
                &p["a_mod_q"],
                &p["online_outputs"][0]["target_mod_q"],
                "online_trial_0",
            )
        };
        let a: MatZq = serde_json::from_str(&a_value.to_string())?;
        let y: MatZq = serde_json::from_str(&y_value.to_string())?;
        let beta: Q =
            serde_json::from_str(&p["assessment"]["accepted_full_norm_bound_squared"].to_string())?;
        let start = Instant::now();
        let result = match public_nearest_plane(&a, &y) {
            Ok((x, support)) => {
                let norm = x.norm_eucl_sqrd()?;
                json!({"case":item["case"],"status":"completed","support":support,"target_source":target_source,
                    "norm_squared":norm,"x_integer":x,"accepted_by_public_equation_and_norm":norm<=beta,
                    "elapsed_ms":start.elapsed().as_secs_f64()*1000.0})
            }
            Err(error) => json!({"case":item["case"],"status":"error","error":error.to_string()}),
        };
        println!(
            "{}: accepted={}",
            item["case"], result["accepted_by_public_equation_and_norm"]
        );
        results.push(result);
    }
    let output_report = json!({"input_report":input,
        "input_report_sha256":sha256(&input_bytes),
        "source_sha256":sha256(include_str!("plt_public_lattice_screen.rs").as_bytes()),
        "scope":"public-only small-subset LLL and exact rational nearest plane; no bit estimate",
        "lll_delta":"99/100","lll_eta":"501/1000","security_bit_target":null,"results":results});
    std::fs::write(output, serde_json::to_string_pretty(&output_report)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn public_lattice_reduction_preserves_coset_and_shortens_large_solution() {
        use std::str::FromStr;
        let q = 2_305_843_009_213_693_951i64;
        let a = MatZq::from((&MatZ::from_str("[[1,0,2,3],[0,1,5,7]]").unwrap(), q));
        let y = MatZq::from((&MatZ::from_str("[[123456789],[987654321]]").unwrap(), q));
        let (x, _) = public_nearest_plane(&a, &y).unwrap();
        assert_eq!(&a * &x, y);
        assert!(
            x.norm_eucl_sqrd().unwrap()
                <= public_solution(&a, &y).unwrap().norm_eucl_sqrd().unwrap()
        );
    }
}

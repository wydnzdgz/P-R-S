// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Independent finite-box reference for a small conditional discrete Gaussian.

use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::{MatQ, Q},
    traits::{MatrixGetEntry, MatrixSetEntry},
};
use qfall_tools::{
    primitive::psf::PSFPerturbation, sample::g_trapdoor::gadget_parameters::GadgetParameters,
};
use serde_json::{Value, json};
use std::{error::Error, str::FromStr};

pub fn run_reference() -> Result<Value, Box<dyn Error>> {
    let mut gp = GadgetParameters::init_default(1, 3);
    gp.m_bar = Z::ONE;
    let sampler = PSFPerturbation {
        gp,
        r: Q::from(4),
        s: Q::from(64),
    };
    let (a, td) = sampler.try_trap_gen()?;
    let sigma = MatQ::from_str("[[3200,800,0],[800,4800,400],[0,400,6400]]")?;
    let factor = sampler.try_compute_sqrt_sigma_2(&td.0, &(Q::from((1, 16)) * &sigma))?;
    let trapdoor = (td.0, factor, td.2);
    let inverse = sigma.inverse().ok_or("singular reference covariance")?;
    let mut inv = [[0.0; 3]; 3];
    for (i, row) in inv.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            let entry: Q = inverse.get_entry(i, j)?;
            *value = f64::from(&entry);
        }
    }
    let coefficients: Vec<i64> = (0..3)
        .map(|j| {
            let value: Z = a.get_entry(0, j).unwrap();
            i64::try_from(value).unwrap()
        })
        .collect();
    let pivot = coefficients
        .iter()
        .position(|&value| value != 0)
        .ok_or("rank deficient reference")?;
    let free: Vec<_> = (0..3).filter(|&j| j != pivot).collect();
    let (u, v) = (free[0], free[1]);
    let limit = 220i64;
    let bins = (4 * limit + 1) as usize;
    let mut reports = Vec::new();
    for target_value in [0i64, 1] {
        let mut weights = vec![vec![0.0; bins]; 4];
        let mut partition = 0.0;
        let mut points = 0usize;
        for x in -limit..=limit {
            for y in -limit..=limit {
                let wanted = ((target_value - coefficients[u] * x - coefficients[v] * y)
                    * coefficients[pivot])
                    .rem_euclid(3);
                let start = -limit + (wanted + limit).rem_euclid(3); // 1^-1=1, 2^-1=2 modulo 3.
                let constant = inv[u][u] * (x * x) as f64
                    + 2.0 * inv[u][v] * (x * y) as f64
                    + inv[v][v] * (y * y) as f64;
                let linear = 2.0 * (inv[u][pivot] * x as f64 + inv[v][pivot] * y as f64);
                for z in (start..=limit).step_by(3) {
                    let mut d = [0i64; 3];
                    d[u] = x;
                    d[v] = y;
                    d[pivot] = z;
                    let mass = (-std::f64::consts::PI
                        * (constant + linear * z as f64 + inv[pivot][pivot] * (z * z) as f64))
                        .exp();
                    partition += mass;
                    points += 1;
                    for (index, value) in [d[0], d[1], d[2], d[0] + d[1]].into_iter().enumerate() {
                        weights[index][(value + 2 * limit) as usize] += mass;
                    }
                }
            }
        }
        let variance_bound = 6800.0f64;
        let tail = 6.0 * variance_bound / (2.0 * std::f64::consts::PI * limit as f64)
            * (-std::f64::consts::PI * (limit * limit) as f64 / variance_bound).exp()
            * (1.0 + variance_bound.sqrt()).powi(2)
            / partition;
        let mut target = MatZq::new(1, 1, 3);
        target.set_entry(0, 0, target_value)?;
        let mut empirical = vec![vec![0usize; bins]; 4];
        let mut outside = 0usize;
        for _ in 0..10000 {
            let sample: MatZ = sampler.try_samp_p(&a, &trapdoor, &target)?;
            if &a * &sample != target {
                return Err("reference sample violates modular equation".into());
            }
            let mut d = [0i64; 3];
            for (i, value) in d.iter_mut().enumerate() {
                let entry: Z = sample.get_entry(i, 0)?;
                *value = i64::try_from(entry)?;
            }
            for (index, value) in [d[0], d[1], d[2], d[0] + d[1]].into_iter().enumerate() {
                let support_limit = if index < 3 { limit } else { 2 * limit };
                if value < -support_limit || value > support_limit {
                    outside += 1;
                }
                empirical[index][(value + 2 * limit).clamp(0, 4 * limit) as usize] += 1;
            }
        }
        let mut ks = Vec::new();
        for (reference, observed) in weights.iter().zip(&empirical) {
            let (mut cdf, mut empirical_cdf, mut difference) = (0.0f64, 0.0f64, 0.0f64);
            for (mass, count) in reference.iter().zip(observed) {
                cdf += mass / partition;
                empirical_cdf += *count as f64 / 10000.0;
                difference = difference.max((cdf - empirical_cdf).abs());
            }
            ks.push(difference);
        }
        let passed = tail < 1e-8 && outside == 0 && ks.iter().all(|&value| value <= 0.04);
        println!(
            "small discrete reference target={target_value}: KS={ks:?}, tail_bound={tail:e}, passed={passed}"
        );
        reports.push(json!({"target":target_value,"samples":10000,"enumerated_points":points,
            "reference_box_limit":limit,"relative_omitted_mass_upper_bound":tail,"empirical_projection_outside_count":outside,
            "ks_for_d1_d2_d3_d1_plus_d2":ks,"ks_limit":0.04,"reference_projection_weights":weights,"empirical_projection_counts":empirical,
            "partition_function_in_box":partition,"screening_passed":passed}));
    }
    let passed = reports
        .iter()
        .all(|report| report["screening_passed"] == true);
    Ok(
        json!({"modulus":3,"c_mod_q":a,"sigma_squared_matrix":sigma,"reports":reports,"screening_passed":passed}),
    )
}

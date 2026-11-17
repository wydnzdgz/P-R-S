// This file is part of qFALL-tools, distributed under MPL-2.0.

//! PLT phase benchmarks; the baseline sampler has its own package benchmark.

use criterion::criterion_main;

pub mod plt;

criterion_main! {plt::benches}

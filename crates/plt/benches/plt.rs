// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Shared phase benchmarks for configurable integer PLT.
//! Filter with PLT2of3, PLT3of5 or PLT5of10. Each uses its first legal recovery set.

use criterion::{Criterion, criterion_group};
use plt::{PltConfig, psamp_pre, ptrap_gen, rec, setup};
use qfall_math::{integer::MatZ, integer_mod_q::MatZq};
use std::hint::black_box;

/// One benchmark implementation; the configuration is the only difference between scales.
fn bench_plt(c: &mut Criterion) {
    for (t, k) in [(2, 3), (3, 5), (5, 10)] {
        let config = PltConfig::new(t, k).unwrap();
        let label = format!("PLT{t}of{k}");
        let mut group = c.benchmark_group(&label);
        group.bench_function("setup", |b| {
            b.iter(|| black_box(setup(black_box(&config)).unwrap()))
        });
        let (public, master) = setup(&config).unwrap();
        group.bench_function("ptrap_gen_one", |b| {
            b.iter(|| black_box(ptrap_gen(&public, &master, black_box(0)).unwrap()))
        });
        group.bench_function("ptrap_gen_all", |b| {
            b.iter(|| {
                for participant in 0..config.participants() {
                    black_box(ptrap_gen(&public, &master, participant).unwrap());
                }
            })
        });
        let partials: Vec<_> = (0..k)
            .map(|id| ptrap_gen(&public, &master, id).unwrap())
            .collect();
        drop(master);
        let width = partials
            .iter()
            .map(|p| p.minimum_sigma_squared())
            .max()
            .unwrap();
        let target = MatZq::sample_uniform(config.target_rows(), 1, config.modulus());
        let set = config.recovery_sets().next().unwrap();
        let directory = format!("target/criterion/{label}");
        std::fs::create_dir_all(&directory).unwrap();
        std::fs::write(format!("{directory}/context.json"), serde_json::to_string_pretty(&serde_json::json!({
            "parameters": {"q":config.modulus(),"n":config.n(),"m":config.m(),"t":t,"k":k},
            "sigma_squared":width, "recovery_set":set, "target":target,
            "initial_width_rule":"power of two s>=64 with s^2>5*worst_case_full_trapdoor_frobenius_squared+1; width=4*s",
            "online_keys":"one random key set and target reused per configuration",
            "offline_keys":"fresh random trapdoors per iteration",
        })).unwrap()).unwrap();
        for &participant in &set {
            group.bench_function(format!("psamp_pre/party{participant}"), |b| {
                b.iter(|| {
                    black_box(
                        psamp_pre(
                            &public,
                            &partials[participant],
                            &set,
                            black_box(&target),
                            width,
                        )
                        .unwrap(),
                    )
                })
            });
        }
        group.bench_function("online_sequential", |b| {
            b.iter(|| {
                let shares: Vec<_> = set
                    .iter()
                    .map(|&id| psamp_pre(&public, &partials[id], &set, &target, width).unwrap())
                    .collect();
                black_box(rec(&public, &shares.iter().collect::<Vec<_>>()).unwrap())
            })
        });
        let shares: Vec<_> = set
            .iter()
            .map(|&id| psamp_pre(&public, &partials[id], &set, &target, width).unwrap())
            .collect();
        let refs: Vec<_> = shares.iter().collect();
        group.bench_function("rec_checked", |b| {
            b.iter(|| black_box(rec(&public, black_box(&refs)).unwrap()))
        });
        group.bench_function("rec_integer_sum_only", |b| {
            b.iter(|| {
                let mut sum = MatZ::new(config.preimage_rows(), 1);
                for share in black_box(&shares) {
                    sum += share.x();
                }
                black_box(sum)
            })
        });
        group.finish();
    }
}

criterion_group!(benches, bench_plt);

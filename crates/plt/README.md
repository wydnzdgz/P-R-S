# PLT

Configurable integer-lattice partial trapdoor research prototype, extracted from
the maintained qFALL-tools baseline. One implementation of Setup, PTrapGen,
PSampPre and Rec serves all supported threshold configurations.

```rust
use plt::{PltConfig, setup};

let config = PltConfig::new(3, 5)?;
let (public, master) = setup(&config)?;
# Ok::<(), plt::PltError>(())
```

From the workspace root:

```bash
cargo test --locked -p plt
cargo test --locked -p plt --examples
cargo run --locked -p plt --example plt -- --threshold 3 --participants 5
cargo bench --locked -p plt --bench benchmarks -- PLT3of5
```

The library code lives in `src/`, integration tests in `tests/`, experimental
drivers in `examples/`, and phase benchmarks in `benches/`.
The MP12 backend is a dependency on the [baseline member](../qfall-tools/README.md).
No protocol or vehicle-network functionality belongs in this crate.

Functional correctness does not approve cryptographic parameters. The evaluated
instances admit public short-preimage counterexamples; the protocol parameter
gate remains closed. See the [acceptance checklist](../../docs/acceptance_checklist.md),
[parameter report](../../docs/plt_parameter_exploration.md), and
[workspace migration record](../../docs/workspace_migration.md).

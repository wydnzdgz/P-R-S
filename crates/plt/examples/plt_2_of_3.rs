// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Compatibility entry point; all arguments and logic are handled by the generic example.

#[path = "plt.rs"]
mod experiment;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    experiment::main()
}

// This file is part of qFALL-tools, distributed under MPL-2.0.

//! Fresh OS entropy per sampling operation. No cached process/fork RNG state.

use qfall_math::error::MathError;
use rand::{TryCryptoRng, rngs::SysRng};

pub(crate) fn fresh_seed() -> Result<[u8; 32], MathError> {
    seed_from(&mut SysRng)
}

fn seed_from(rng: &mut impl TryCryptoRng) -> Result<[u8; 32], MathError> {
    let mut seed = [0u8; 32];
    rng.try_fill_bytes(&mut seed).map_err(|_| {
        MathError::InvalidIntegerInput(
            "operating-system entropy unavailable; sampling aborted".into(),
        )
    })?;
    Ok(seed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::TryRng;
    struct Failing;
    impl TryRng for Failing {
        type Error = std::io::Error;
        fn try_next_u32(&mut self) -> Result<u32, Self::Error> {
            Err(std::io::Error::other("unavailable"))
        }
        fn try_next_u64(&mut self) -> Result<u64, Self::Error> {
            Err(std::io::Error::other("unavailable"))
        }
        fn try_fill_bytes(&mut self, _: &mut [u8]) -> Result<(), Self::Error> {
            Err(std::io::Error::other("unavailable"))
        }
    }
    impl TryCryptoRng for Failing {}
    #[test]
    fn entropy_failure_is_propagated_without_fallback() {
        assert!(seed_from(&mut Failing).is_err());
    }
    #[test]
    fn concurrent_operations_get_independent_os_seeds() {
        let threads: Vec<_> = (0..4)
            .map(|_| {
                std::thread::spawn(|| (0..16).map(|_| fresh_seed().unwrap()).collect::<Vec<_>>())
            })
            .collect();
        let seeds: std::collections::HashSet<_> = threads
            .into_iter()
            .flat_map(|t| t.join().unwrap())
            .collect();
        assert_eq!(seeds.len(), 64);
    }
}

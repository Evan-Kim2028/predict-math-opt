/// On-chain gas measurement entry functions.
/// Each function calls the target computation in a loop to break out of
/// Sui's minimum gas bucket (1M MIST) and get real differentiated costs.
module predict_math_opt::gas_measure;

use predict_math_opt::original_cdf;
use predict_math_opt::piecewise_cdf;
use predict_math_opt::optimized_math;
use predict_math_opt::math_utils;
use predict_math_opt::compute_nd2;

const FLOAT_SCALING: u64 = 1_000_000_000;

/// 50 calls to original A&S CDF at varied x values.
public entry fun measure_original_cdf_50() {
    let mut x: u64 = 100_000_000;
    let mut i = 0;
    while (i < 50) {
        let _ = original_cdf::normal_cdf(x, i % 3 == 0);
        x = x + 70_000_000;
        if (x > 3_500_000_000) { x = x - 3_400_000_000; };
        i = i + 1;
    };
}

/// 50 calls to piecewise cubic CDF at same x values.
public entry fun measure_piecewise_cdf_50() {
    let mut x: u64 = 100_000_000;
    let mut i = 0;
    while (i < 50) {
        let _ = piecewise_cdf::normal_cdf(x, i % 3 == 0);
        x = x + 70_000_000;
        if (x > 3_500_000_000) { x = x - 3_400_000_000; };
        i = i + 1;
    };
}

/// 50 calls to original ln.
public entry fun measure_original_ln_50() {
    let mut x: u64 = 500_000_000;
    let mut i = 0;
    while (i < 50) {
        let (_val, _neg) = original_cdf::ln(x);
        x = x + 100_000_000;
        if (x > 5_000_000_000) { x = 500_000_000; };
        i = i + 1;
    };
}

/// 50 calls to Horner ln.
public entry fun measure_optimized_ln_50() {
    let mut x: u64 = 500_000_000;
    let mut i = 0;
    while (i < 50) {
        let (_val, _neg) = optimized_math::ln(x);
        x = x + 100_000_000;
        if (x > 5_000_000_000) { x = 500_000_000; };
        i = i + 1;
    };
}

/// 50 calls to original sqrt.
public entry fun measure_original_sqrt_50() {
    let mut x: u64 = 100_000_000;
    let mut i = 0;
    while (i < 50) {
        let _ = math_utils::sqrt(x, FLOAT_SCALING);
        x = x + 200_000_000;
        if (x > 10_000_000_000) { x = 100_000_000; };
        i = i + 1;
    };
}

/// 50 calls to optimized sqrt.
public entry fun measure_optimized_sqrt_50() {
    let mut x: u64 = 100_000_000;
    let mut i = 0;
    while (i < 50) {
        let _ = optimized_math::sqrt(x, FLOAT_SCALING);
        x = x + 200_000_000;
        if (x > 10_000_000_000) { x = 100_000_000; };
        i = i + 1;
    };
}

/// 20 calls to full original compute_nd2 pipeline.
public entry fun measure_original_nd2_20() {
    let svi = compute_nd2::new_svi(
        50_000_000, 200_000_000, 300_000_000, true, 10_000_000, false, 100_000_000,
    );
    let forward = 50_000_000_000_000;
    let mut i = 0u64;
    while (i < 20) {
        let strike = forward - 5_000_000_000_000 + i * 500_000_000_000;
        let _ = compute_nd2::compute_nd2_original(forward, strike, &svi, i % 2 == 0);
        i = i + 1;
    };
}

/// 20 calls to full optimized compute_nd2 pipeline.
public entry fun measure_optimized_nd2_20() {
    let svi = compute_nd2::new_svi(
        50_000_000, 200_000_000, 300_000_000, true, 10_000_000, false, 100_000_000,
    );
    let forward = 50_000_000_000_000;
    let mut i = 0u64;
    while (i < 20) {
        let strike = forward - 5_000_000_000_000 + i * 500_000_000_000;
        let _ = compute_nd2::compute_nd2_optimized(forward, strike, &svi, i % 2 == 0);
        i = i + 1;
    };
}

/// 100 calls to optimized compute_nd2 (the vault rebalance scenario).
public entry fun measure_optimized_nd2_100() {
    let svi = compute_nd2::new_svi(
        50_000_000, 200_000_000, 300_000_000, true, 10_000_000, false, 100_000_000,
    );
    let forward = 50_000_000_000_000;
    let mut i = 0u64;
    while (i < 100) {
        let strike = forward - 5_000_000_000_000 + i * 100_000_000_000;
        let _ = compute_nd2::compute_nd2_optimized(forward, strike, &svi, i % 2 == 0);
        i = i + 1;
    };
}

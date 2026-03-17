/// Full compute_nd2 pipeline — original vs fully optimized.
///
/// This module implements the complete SVI + Black-Scholes N(d2) computation
/// from deepbook_predict::oracle::compute_nd2, in both original and optimized form.
///
/// The original uses: loop-based ln, loop-based sqrt, A&S CDF with Taylor exp
/// The optimized uses: Horner ln, unrolled sqrt, piecewise cubic CDF
///
/// Forked from: https://github.com/MystenLabs/deepbookv3/blob/main/packages/predict/sources/oracle.move#L378
module predict_math_opt::compute_nd2;

use predict_math_opt::math_utils;
use predict_math_opt::original_cdf;
use predict_math_opt::piecewise_cdf;
use predict_math_opt::optimized_math;

const FLOAT_SCALING: u64 = 1_000_000_000;

// === SVI Parameters (passed as struct for readability) ===

public struct SVIParams has copy, drop {
    a: u64,
    b: u64,
    rho: u64,
    rho_negative: bool,
    m: u64,
    m_negative: bool,
    sigma: u64,
}

public fun new_svi(
    a: u64, b: u64, rho: u64, rho_negative: bool,
    m: u64, m_negative: bool, sigma: u64,
): SVIParams {
    SVIParams { a, b, rho, rho_negative, m, m_negative, sigma }
}

// ============================================================
// Original compute_nd2 (using A&S CDF, loop-based ln/sqrt)
// ============================================================

/// Original compute_nd2 — exact port from deepbook_predict::oracle.
/// Uses original_cdf::normal_cdf (A&S), original_cdf::ln, math_utils::sqrt.
public fun compute_nd2_original(
    forward: u64,
    strike: u64,
    svi: &SVIParams,
    is_up: bool,
): u64 {
    // SVI: compute total variance from log-moneyness
    let (k, k_neg) = original_cdf::ln(math_utils::div(strike, forward));

    let (k_minus_m, km_neg) = math_utils::sub_signed(k, k_neg, svi.m, svi.m_negative);

    let sq = math_utils::sqrt(
        math_utils::mul(k_minus_m, k_minus_m) + math_utils::mul(svi.sigma, svi.sigma),
        FLOAT_SCALING,
    );

    let (rho_km, rho_km_neg) = math_utils::mul_signed(
        svi.rho, svi.rho_negative, k_minus_m, km_neg,
    );

    let (inner, _inner_neg) = math_utils::add_signed(rho_km, rho_km_neg, sq, false);
    let total_var = svi.a + math_utils::mul(svi.b, inner);

    // d2 = (-k - total_var/2) / sqrt(total_var)
    let sqrt_var = math_utils::sqrt(total_var, FLOAT_SCALING);
    let (d2, d2_neg) = math_utils::sub_signed(k, !k_neg, total_var / 2, false);
    let d2 = math_utils::div(d2, sqrt_var);
    let cdf_neg = if (is_up) { d2_neg } else { !d2_neg };

    original_cdf::normal_cdf(d2, cdf_neg)
}

// ============================================================
// Optimized compute_nd2 (piecewise CDF, Horner ln, fast sqrt)
// ============================================================

/// Optimized compute_nd2 — uses all three optimizations:
/// 1. Piecewise cubic CDF (no Taylor exp)
/// 2. Horner-form ln (no loop)
/// 3. Unrolled Newton sqrt (no loop)
public fun compute_nd2_optimized(
    forward: u64,
    strike: u64,
    svi: &SVIParams,
    is_up: bool,
): u64 {
    // SVI: same math, different implementations
    let (k, k_neg) = optimized_math::ln(math_utils::div(strike, forward));

    let (k_minus_m, km_neg) = math_utils::sub_signed(k, k_neg, svi.m, svi.m_negative);

    let sq = optimized_math::sqrt(
        math_utils::mul(k_minus_m, k_minus_m) + math_utils::mul(svi.sigma, svi.sigma),
        FLOAT_SCALING,
    );

    let (rho_km, rho_km_neg) = math_utils::mul_signed(
        svi.rho, svi.rho_negative, k_minus_m, km_neg,
    );

    let (inner, _inner_neg) = math_utils::add_signed(rho_km, rho_km_neg, sq, false);
    let total_var = svi.a + math_utils::mul(svi.b, inner);

    // d2 = (-k - total_var/2) / sqrt(total_var)
    let sqrt_var = optimized_math::sqrt(total_var, FLOAT_SCALING);
    let (d2, d2_neg) = math_utils::sub_signed(k, !k_neg, total_var / 2, false);
    let d2 = math_utils::div(d2, sqrt_var);
    let cdf_neg = if (is_up) { d2_neg } else { !d2_neg };

    piecewise_cdf::normal_cdf(d2, cdf_neg)
}

// ============================================================
// Tests: full pipeline regression
// ============================================================

#[test]
/// Original and optimized should produce similar results for realistic SVI params.
/// Test with BTC-like parameters.
fun test_nd2_original_vs_optimized_btc() {
    // Realistic SVI params (BTC-like)
    let svi = new_svi(
        50_000_000,    // a = 0.05 (base variance)
        200_000_000,   // b = 0.2 (smile slope)
        300_000_000,   // rho = -0.3 (skew)
        true,          // rho negative
        10_000_000,    // m = 0.01 (shift)
        false,         // m positive
        100_000_000,   // sigma = 0.1 (smoothness)
    );

    let forward = 50_000_000_000_000; // $50,000 in some scaling

    // Test at various strikes (0.8x to 1.2x of forward = ITM to OTM)
    let strikes: vector<u64> = vector[
        40_000_000_000_000, // 0.8x (deep ITM for UP)
        45_000_000_000_000, // 0.9x
        50_000_000_000_000, // 1.0x (ATM)
        55_000_000_000_000, // 1.1x
        60_000_000_000_000, // 1.2x (deep OTM for UP)
    ];

    let mut i = 0;
    while (i < strikes.length()) {
        let strike = strikes[i];
        let orig = compute_nd2_original(forward, strike, &svi, true);
        let opt = compute_nd2_optimized(forward, strike, &svi, true);

        // Both should be in valid CDF range [0, FLOAT_SCALING]
        assert!(orig <= FLOAT_SCALING, 100 + i);
        assert!(opt <= FLOAT_SCALING, 200 + i);

        // Results should be close (< 2bp difference)
        let diff = if (orig > opt) { orig - opt } else { opt - orig };
        assert!(diff < 200_000, 300 + i); // < 2bp

        i = i + 1;
    };
}

#[test]
/// Test with ETH-like parameters.
fun test_nd2_original_vs_optimized_eth() {
    let svi = new_svi(
        80_000_000,    // a = 0.08
        150_000_000,   // b = 0.15
        250_000_000,   // rho = -0.25
        true,          // rho negative
        5_000_000,     // m = 0.005
        true,          // m negative
        150_000_000,   // sigma = 0.15
    );

    let forward = 3_000_000_000_000; // $3,000

    let strikes: vector<u64> = vector[
        2_400_000_000_000, // 0.8x
        2_700_000_000_000, // 0.9x
        3_000_000_000_000, // ATM
        3_300_000_000_000, // 1.1x
        3_600_000_000_000, // 1.2x
    ];

    let mut i = 0;
    while (i < strikes.length()) {
        let strike = strikes[i];
        let orig_up = compute_nd2_original(forward, strike, &svi, true);
        let opt_up = compute_nd2_optimized(forward, strike, &svi, true);
        let orig_down = compute_nd2_original(forward, strike, &svi, false);
        let opt_down = compute_nd2_optimized(forward, strike, &svi, false);

        // UP and DOWN should be close between implementations
        let diff_up = if (orig_up > opt_up) { orig_up - opt_up } else { opt_up - orig_up };
        let diff_down = if (orig_down > opt_down) { orig_down - opt_down } else { opt_down - orig_down };
        assert!(diff_up < 200_000, 100 + i);
        assert!(diff_down < 200_000, 200 + i);

        // UP + DOWN should be close to 1.0 (complementary)
        let sum_orig = orig_up + orig_down;
        let sum_diff = if (sum_orig > FLOAT_SCALING) { sum_orig - FLOAT_SCALING } else { FLOAT_SCALING - sum_orig };
        assert!(sum_diff < 1_000_000, 300 + i); // < 10bp tolerance for complement

        i = i + 1;
    };
}

#[test]
/// ATM option: strike = forward → d2 depends only on variance.
fun test_nd2_atm() {
    let svi = new_svi(
        100_000_000, 100_000_000, 200_000_000, true,
        0, false, 100_000_000,
    );
    let forward = 1_000_000_000_000;
    let strike = forward; // ATM

    let orig = compute_nd2_original(forward, strike, &svi, true);
    let opt = compute_nd2_optimized(forward, strike, &svi, true);

    // ATM: both should be near 0.5 (since d2 is small)
    assert!(orig > 200_000_000 && orig < 800_000_000, 0);
    assert!(opt > 200_000_000 && opt < 800_000_000, 1);

    let diff = if (orig > opt) { orig - opt } else { opt - orig };
    assert!(diff < 200_000, 2);
}

// ============================================================
// Gas benchmarks: full pipeline
// ============================================================

#[test]
fun bench_compute_nd2_original_single() {
    let svi = new_svi(50_000_000, 200_000_000, 300_000_000, true, 10_000_000, false, 100_000_000);
    let result = compute_nd2_original(50_000_000_000_000, 55_000_000_000_000, &svi, true);
    assert!(result > 0 && result < FLOAT_SCALING);
}

#[test]
fun bench_compute_nd2_optimized_single() {
    let svi = new_svi(50_000_000, 200_000_000, 300_000_000, true, 10_000_000, false, 100_000_000);
    let result = compute_nd2_optimized(50_000_000_000_000, 55_000_000_000_000, &svi, true);
    assert!(result > 0 && result < FLOAT_SCALING);
}

#[test]
/// Simulate 20-call vault rebalance with original.
/// (100 calls times out — confirming the performance problem Aslan described)
fun bench_compute_nd2_original_20() {
    let svi = new_svi(50_000_000, 200_000_000, 300_000_000, true, 10_000_000, false, 100_000_000);
    let forward = 50_000_000_000_000;
    let mut sum = 0u64;
    let mut i = 0u64;
    while (i < 20) {
        let strike = forward - 5_000_000_000_000 + i * 500_000_000_000;
        sum = sum + compute_nd2_original(forward, strike, &svi, i % 2 == 0);
        i = i + 1;
    };
    assert!(sum > 0);
}

#[test]
/// Simulate 100-call vault rebalance with optimized.
fun bench_compute_nd2_optimized_100() {
    let svi = new_svi(50_000_000, 200_000_000, 300_000_000, true, 10_000_000, false, 100_000_000);
    let forward = 50_000_000_000_000;
    let mut sum = 0u64;
    let mut i = 0u64;
    while (i < 100) {
        let strike = forward - 5_000_000_000_000 + i * 100_000_000_000;
        sum = sum + compute_nd2_optimized(forward, strike, &svi, i % 2 == 0);
        i = i + 1;
    };
    assert!(sum > 0);
}

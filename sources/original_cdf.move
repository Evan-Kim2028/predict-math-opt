/// Original CDF implementation — exact port of deepbook_predict::math.
/// Abramowitz & Stegun approximation 26.2.17.
///
/// Forked from: https://github.com/MystenLabs/deepbookv3/blob/main/packages/predict/sources/helper/math.move
///
/// This is the EXACT algorithm currently used in DeepBook Predict production.
/// Ported here as a baseline for gas comparison against optimized alternatives.
/// 674 bytecode instructions per call (from `sui move disassemble`).
module predict_math_opt::original_cdf;

use predict_math_opt::math_utils;

const FLOAT_SCALING: u64 = 1_000_000_000;

/// Standard normal CDF using Abramowitz & Stegun (26.2.17).
/// Ported directly from deepbook_predict::math::normal_cdf.
public fun normal_cdf(x: u64, x_negative: bool): u64 {
    if (x > 8 * FLOAT_SCALING) {
        return if (x_negative) { 0 } else { FLOAT_SCALING }
    };

    let t = cdf_t(x);
    let poly = cdf_poly(t);
    let pdf = cdf_pdf(x);
    let complement = math_utils::mul(pdf, poly);

    let cdf = if (FLOAT_SCALING > complement) {
        FLOAT_SCALING - complement
    } else {
        0
    };

    if (x_negative) { FLOAT_SCALING - cdf } else { cdf }
}

/// Natural logarithm of x (in FLOAT_SCALING).
/// Ported from deepbook_predict::math::ln.
public fun ln(x: u64): (u64, bool) {
    assert!(x > 0, 0);
    if (x == FLOAT_SCALING) return (0, false);

    if (x < FLOAT_SCALING) {
        let inv = math_utils::div(FLOAT_SCALING, x);
        let (result, _) = ln(inv);
        return (result, true)
    };

    let (y, n) = normalize(x);
    let z = log_ratio(y);
    let ln_y = ln_series(z);
    let result = n * 693_147_181 + ln_y;

    (result, false)
}

/// Exponential function e^(+-x).
/// Ported from deepbook_predict::math::exp.
public fun exp(x: u64, x_negative: bool): u64 {
    if (x == 0) return FLOAT_SCALING;

    let (r, mut n) = reduce_exp(x);
    let exp_r = exp_series(r);

    if (x_negative) {
        let mut result = math_utils::div(FLOAT_SCALING, exp_r);
        if (n >= 32) { result = result >> 32; if (result == 0) return 0; n = n - 32; };
        if (n >= 16) { result = result >> 16; if (result == 0) return 0; n = n - 16; };
        if (n >= 8) { result = result >> 8; if (result == 0) return 0; n = n - 8; };
        if (n >= 4) { result = result >> 4; if (result == 0) return 0; n = n - 4; };
        if (n >= 2) { result = result >> 2; if (result == 0) return 0; n = n - 2; };
        if (n >= 1) { result = result >> 1; };
        result
    } else {
        let mut result = exp_r;
        if (n >= 32) { result = result << 32; n = n - 32; };
        if (n >= 16) { result = result << 16; n = n - 16; };
        if (n >= 8) { result = result << 8; n = n - 8; };
        if (n >= 4) { result = result << 4; n = n - 4; };
        if (n >= 2) { result = result << 2; n = n - 2; };
        if (n >= 1) { result = result << 1; };
        result
    }
}

// === Internal helpers (ported verbatim) ===

fun cdf_t(x: u64): u64 {
    math_utils::div(FLOAT_SCALING, FLOAT_SCALING + math_utils::mul(231_641_900, x))
}

fun cdf_pdf(x: u64): u64 {
    let x_sq_half = math_utils::mul(x, x) / 2;
    math_utils::mul(exp(x_sq_half, true), 398_942_280)
}

fun cdf_poly(t: u64): u64 {
    let t2 = math_utils::mul(t, t);
    let t3 = math_utils::mul(t2, t);
    let t4 = math_utils::mul(t3, t);
    let t5 = math_utils::mul(t4, t);

    let pos =
        math_utils::mul(319_381_530, t)
        + math_utils::mul(1_781_477_937, t3)
        + math_utils::mul(1_330_274_429, t5);

    let neg = math_utils::mul(356_563_782, t2)
        + math_utils::mul(1_821_255_978, t4);

    pos - neg
}

fun reduce_exp(x: u64): (u64, u64) {
    let n = x / 693_147_181;
    let r = x - n * 693_147_181;
    (r, n)
}

fun exp_series(r: u64): u64 {
    let mut sum = FLOAT_SCALING;
    let mut term = FLOAT_SCALING;
    let mut k: u64 = 1;
    while (k <= 12) {
        term = math_utils::div(math_utils::mul(term, r), k * FLOAT_SCALING);
        if (term == 0) break;
        sum = sum + term;
        k = k + 1;
    };
    sum
}

fun ln_series(z: u64): u64 {
    let z2 = math_utils::mul(z, z);
    let mut term = z;
    let mut sum = 0;
    let mut k: u64 = 1;
    while (k <= 13) {
        sum = sum + math_utils::div(term, k * FLOAT_SCALING);
        term = math_utils::mul(term, z2);
        k = k + 2;
    };
    math_utils::mul(2 * FLOAT_SCALING, sum)
}

fun log_ratio(y: u64): u64 {
    math_utils::div(y - FLOAT_SCALING, y + FLOAT_SCALING)
}

fun normalize(x: u64): (u64, u64) {
    let mut y = x;
    let mut n: u64 = 0;
    let scale = FLOAT_SCALING;

    if (y >> 32 >= scale) { y = y >> 32; n = n + 32; };
    if (y >> 16 >= scale) { y = y >> 16; n = n + 16; };
    if (y >> 8 >= scale) { y = y >> 8; n = n + 8; };
    if (y >> 4 >= scale) { y = y >> 4; n = n + 4; };
    if (y >> 2 >= scale) { y = y >> 2; n = n + 2; };
    if (y >> 1 >= scale) { y = y >> 1; n = n + 1; };

    (y, n)
}

// === Tests ===

#[test]
fun test_cdf_at_zero() {
    let result = normal_cdf(0, false);
    let diff = if (result > 500_000_000) { result - 500_000_000 } else { 500_000_000 - result };
    assert!(diff < 1_000, 0); // < 0.000001
}

#[test]
fun test_cdf_at_one() {
    let result = normal_cdf(1_000_000_000, false);
    // Phi(1) ≈ 0.8413
    let expected = 841_344_746;
    let diff = if (result > expected) { result - expected } else { expected - result };
    assert!(diff < 1_000_000, 1); // < 0.001 tolerance
}

#[test]
fun test_cdf_at_two() {
    let result = normal_cdf(2_000_000_000, false);
    // Phi(2) ≈ 0.9772
    let expected = 977_249_868;
    let diff = if (result > expected) { result - expected } else { expected - result };
    assert!(diff < 1_000_000, 2);
}

#[test]
fun test_cdf_negative() {
    let pos = normal_cdf(1_000_000_000, false);
    let neg = normal_cdf(1_000_000_000, true);
    // Phi(-x) + Phi(x) = 1
    let sum = pos + neg;
    let diff = if (sum > FLOAT_SCALING) { sum - FLOAT_SCALING } else { FLOAT_SCALING - sum };
    assert!(diff < 1_000, 3);
}

#[test]
fun test_exp_zero() {
    assert!(exp(0, false) == FLOAT_SCALING);
}

#[test]
fun test_exp_one() {
    let result = exp(FLOAT_SCALING, false);
    // e ≈ 2.71828
    let expected = 2_718_281_828;
    let diff = if (result > expected) { result - expected } else { expected - result };
    assert!(diff < 1_000_000, 0);
}

#[test]
fun test_ln_e() {
    // ln(e) = 1.0; e ≈ 2_718_281_828
    let (result, neg) = ln(2_718_281_828);
    assert!(!neg);
    let diff = if (result > FLOAT_SCALING) { result - FLOAT_SCALING } else { FLOAT_SCALING - result };
    assert!(diff < 5_000_000, 0); // ~0.5% tolerance (ln precision)
}

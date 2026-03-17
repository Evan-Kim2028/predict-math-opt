/// Standalone fixed-point math helpers (ported from deepbook::math).
/// FLOAT_SCALING = 1e9 throughout.
///
/// Forked from: https://github.com/MystenLabs/deepbookv3/blob/main/packages/deepbook/sources/helper/math.move
/// Self-contained so this package has no dependency on deepbook.
#[allow(unused_const)]
module predict_math_opt::math_utils;

const FLOAT_SCALING: u64 = 1_000_000_000;

/// Fixed-point multiply: (a * b) / FLOAT_SCALING
public fun mul(a: u64, b: u64): u64 {
    ((a as u128) * (b as u128) / (FLOAT_SCALING as u128) as u64)
}

/// Fixed-point divide: (a * FLOAT_SCALING) / b
public fun div(a: u64, b: u64): u64 {
    ((a as u128) * (FLOAT_SCALING as u128) / (b as u128) as u64)
}

/// Integer square root with scaling factor.
/// Returns sqrt(x) adjusted by sqrt(scaling).
public fun sqrt(x: u64, scaling: u64): u64 {
    let scaled = (x as u128) * (scaling as u128);
    (sqrt_u128(scaled) as u64)
}

/// Integer square root via Newton-Raphson.
fun sqrt_u128(x: u128): u128 {
    if (x == 0) return 0;
    let mut guess = (x as u256);
    let mut result = (x as u256);
    // Initial guess: half the bit length
    result = (result + 1) / 2;
    let mut i = 0;
    while (i < 128) {
        if (result >= guess) break;
        guess = result;
        result = (result + (x as u256) / result) / 2;
        i = i + 1;
    };
    (guess as u128)
}

public fun float_scaling(): u64 { FLOAT_SCALING }

/// Signed addition: (a, a_neg) + (b, b_neg) = (result, result_neg)
public fun add_signed(a: u64, a_neg: bool, b: u64, b_neg: bool): (u64, bool) {
    if (a_neg == b_neg) {
        let sum = a + b;
        if (sum == 0) (0, false) else (sum, a_neg)
    } else {
        if (a >= b) {
            let diff = a - b;
            if (diff == 0) (0, false) else (diff, a_neg)
        } else {
            let diff = b - a;
            if (diff == 0) (0, false) else (diff, b_neg)
        }
    }
}

/// Signed subtraction: (a, a_neg) - (b, b_neg)
public fun sub_signed(a: u64, a_neg: bool, b: u64, b_neg: bool): (u64, bool) {
    add_signed(a, a_neg, b, !b_neg)
}

/// Signed multiplication: (a, a_neg) * (b, b_neg)
public fun mul_signed(a: u64, a_neg: bool, b: u64, b_neg: bool): (u64, bool) {
    let product = mul(a, b);
    if (product == 0) return (0, false);
    (product, a_neg != b_neg)
}

// === Tests ===

#[test]
fun test_mul() {
    // 2.0 * 3.0 = 6.0
    assert!(mul(2_000_000_000, 3_000_000_000) == 6_000_000_000);
    // 0.5 * 0.5 = 0.25
    assert!(mul(500_000_000, 500_000_000) == 250_000_000);
}

#[test]
fun test_div() {
    // 6.0 / 3.0 = 2.0
    assert!(div(6_000_000_000, 3_000_000_000) == 2_000_000_000);
    // 1.0 / 2.0 = 0.5
    assert!(div(1_000_000_000, 2_000_000_000) == 500_000_000);
}

#[test]
fun test_sqrt() {
    // sqrt(4.0) = 2.0
    let result = sqrt(4_000_000_000, FLOAT_SCALING);
    let diff = if (result > 2_000_000_000) { result - 2_000_000_000 } else { 2_000_000_000 - result };
    assert!(diff < 2); // Allow 1 unit rounding
}

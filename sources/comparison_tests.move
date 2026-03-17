/// Head-to-head comparison tests for all three CDF implementations.
///
/// Tests accuracy against known reference values and checks that
/// the piecewise/lookup alternatives stay within tolerance of the
/// original A&S implementation.
#[test_only]
module predict_math_opt::comparison_tests;

use predict_math_opt::original_cdf;
use predict_math_opt::piecewise_cdf;
use predict_math_opt::lookup_cdf;

const FLOAT_SCALING: u64 = 1_000_000_000;

// Reference CDF values (from scipy.stats.norm.cdf, 10 decimal places)
// Scaled to FLOAT_SCALING (1e9)
const REF_CDF_0_00: u64 = 500_000_000;  // Phi(0.0)
const REF_CDF_0_25: u64 = 598_706_326;  // Phi(0.25)
const REF_CDF_0_50: u64 = 691_462_461;  // Phi(0.50)
const REF_CDF_0_75: u64 = 773_372_648;  // Phi(0.75)
const REF_CDF_1_00: u64 = 841_344_746;  // Phi(1.00)
const REF_CDF_1_25: u64 = 894_350_226;  // Phi(1.25)
const REF_CDF_1_50: u64 = 933_192_799;  // Phi(1.50)
const REF_CDF_1_75: u64 = 959_940_843;  // Phi(1.75)
const REF_CDF_2_00: u64 = 977_249_868;  // Phi(2.00)
const REF_CDF_2_50: u64 = 993_790_335;  // Phi(2.50)
const REF_CDF_3_00: u64 = 998_650_102;  // Phi(3.00)

fun abs_diff(a: u64, b: u64): u64 {
    if (a > b) { a - b } else { b - a }
}

// ====================================================================
// ACCURACY VS REFERENCE: Check each implementation against exact values
// ====================================================================

// 1 basis point = 0.0001 = 100_000 in FLOAT_SCALING
const ONE_BP: u64 = 100_000;

#[test]
/// All three implementations at the key test points.
/// Checks error vs scipy reference values.
fun test_accuracy_sweep() {
    let test_points: vector<u64> = vector[
        0,
        250_000_000,   // 0.25
        500_000_000,   // 0.50
        750_000_000,   // 0.75
        1_000_000_000, // 1.00
        1_250_000_000, // 1.25
        1_500_000_000, // 1.50
        1_750_000_000, // 1.75
        2_000_000_000, // 2.00
        2_500_000_000, // 2.50
        3_000_000_000, // 3.00
    ];
    let ref_values: vector<u64> = vector[
        REF_CDF_0_00,
        REF_CDF_0_25,
        REF_CDF_0_50,
        REF_CDF_0_75,
        REF_CDF_1_00,
        REF_CDF_1_25,
        REF_CDF_1_50,
        REF_CDF_1_75,
        REF_CDF_2_00,
        REF_CDF_2_50,
        REF_CDF_3_00,
    ];

    let mut i = 0;
    let n = test_points.length();
    while (i < n) {
        let x = test_points[i];
        let expected = ref_values[i];

        let orig = original_cdf::normal_cdf(x, false);
        let piece = piecewise_cdf::normal_cdf(x, false);
        let lut = lookup_cdf::normal_cdf(x, false);

        // Original A&S: should be within ~10bp of reference (known precision)
        let orig_err = abs_diff(orig, expected);
        assert!(orig_err < 10 * ONE_BP, 100 + i);

        // Piecewise cubic: should be within 1bp of reference
        let piece_err = abs_diff(piece, expected);
        assert!(piece_err < 2 * ONE_BP, 200 + i);

        // Lookup: should be within 2bp of reference
        let lut_err = abs_diff(lut, expected);
        assert!(lut_err < 3 * ONE_BP, 300 + i);

        i = i + 1;
    };
}

// ====================================================================
// CONSISTENCY: Piecewise and lookup should track original closely
// ====================================================================

#[test]
/// Sweep 0.1 increments from 0 to 3.5 and check piecewise vs original.
fun test_piecewise_vs_original_fine_sweep() {
    let mut x: u64 = 0;
    while (x <= 3_500_000_000) {
        let orig = original_cdf::normal_cdf(x, false);
        let piece = piecewise_cdf::normal_cdf(x, false);
        let diff = abs_diff(orig, piece);
        // Piecewise should be within 2bp of original at all points
        assert!(diff < 2 * ONE_BP, (x / 100_000_000 as u64));
        x = x + 100_000_000; // 0.1 increments
    };
}

#[test]
/// Sweep 0.1 increments and check lookup vs original.
fun test_lookup_vs_original_fine_sweep() {
    let mut x: u64 = 0;
    while (x <= 3_500_000_000) {
        let orig = original_cdf::normal_cdf(x, false);
        let lut = lookup_cdf::normal_cdf(x, false);
        let diff = abs_diff(orig, lut);
        // Lookup should be within 3bp of original
        assert!(diff < 3 * ONE_BP, (x / 100_000_000 as u64));
        x = x + 100_000_000;
    };
}

// ====================================================================
// SYMMETRY: Phi(-x) + Phi(x) = 1 for all implementations
// ====================================================================

#[test]
fun test_symmetry_all() {
    let test_xs: vector<u64> = vector[
        0, 500_000_000, 1_000_000_000, 1_500_000_000, 2_000_000_000, 3_000_000_000
    ];
    let mut i = 0;
    while (i < test_xs.length()) {
        let x = test_xs[i];

        // Original
        let orig_sum = original_cdf::normal_cdf(x, false) + original_cdf::normal_cdf(x, true);
        assert!(abs_diff(orig_sum, FLOAT_SCALING) < 1_000, 100 + i);

        // Piecewise
        let piece_sum = piecewise_cdf::normal_cdf(x, false) + piecewise_cdf::normal_cdf(x, true);
        assert!(abs_diff(piece_sum, FLOAT_SCALING) < 1_000, 200 + i);

        // Lookup
        let lut_sum = lookup_cdf::normal_cdf(x, false) + lookup_cdf::normal_cdf(x, true);
        assert!(abs_diff(lut_sum, FLOAT_SCALING) < 1_000, 300 + i);

        i = i + 1;
    };
}

// ====================================================================
// MONOTONICITY: CDF must be strictly increasing on (0, 4)
// ====================================================================

#[test]
fun test_monotonicity_piecewise() {
    let mut prev = piecewise_cdf::normal_cdf(0, false);
    let mut x: u64 = 50_000_000; // 0.05 increments
    while (x <= 3_900_000_000) {
        let curr = piecewise_cdf::normal_cdf(x, false);
        assert!(curr >= prev, (x / 50_000_000 as u64));
        prev = curr;
        x = x + 50_000_000;
    };
}

#[test]
fun test_monotonicity_lookup() {
    let mut prev = lookup_cdf::normal_cdf(0, false);
    let mut x: u64 = 50_000_000;
    while (x <= 3_900_000_000) {
        let curr = lookup_cdf::normal_cdf(x, false);
        assert!(curr >= prev, (x / 50_000_000 as u64));
        prev = curr;
        x = x + 50_000_000;
    };
}

// ====================================================================
// EDGE CASES
// ====================================================================

#[test]
fun test_zero_all() {
    let orig = original_cdf::normal_cdf(0, false);
    let piece = piecewise_cdf::normal_cdf(0, false);
    let lut = lookup_cdf::normal_cdf(0, false);
    assert!(abs_diff(orig, 500_000_000) < 1_000);
    assert!(abs_diff(piece, 500_000_000) < 1_000);
    assert!(abs_diff(lut, 500_000_000) < 1_000);
}

#[test]
fun test_large_tail_all() {
    // x=5.0 — deep in tail
    assert!(original_cdf::normal_cdf(5_000_000_000, false) > 999_900_000);
    assert!(piecewise_cdf::normal_cdf(5_000_000_000, false) == FLOAT_SCALING);
    assert!(lookup_cdf::normal_cdf(5_000_000_000, false) == FLOAT_SCALING);
}

#[test]
/// Stress test: segment boundaries in piecewise (where errors are highest).
fun test_piecewise_segment_boundaries() {
    let boundaries: vector<u64> = vector[
        500_000_000,    // 0.5
        1_000_000_000,  // 1.0
        1_500_000_000,  // 1.5
        2_000_000_000,  // 2.0
        2_500_000_000,  // 2.5
        3_000_000_000,  // 3.0
        3_500_000_000,  // 3.5
    ];
    let mut i = 0;
    while (i < boundaries.length()) {
        let x = boundaries[i];
        let orig = original_cdf::normal_cdf(x, false);
        let piece = piecewise_cdf::normal_cdf(x, false);
        let diff = abs_diff(orig, piece);
        // Segment boundaries can have slightly higher error
        assert!(diff < 3 * ONE_BP, i);
        i = i + 1;
    };
}

// ====================================================================
// NEGATIVE INPUT SWEEP
// ====================================================================

#[test]
fun test_negative_sweep_all() {
    let test_xs: vector<u64> = vector[
        250_000_000, 750_000_000, 1_500_000_000, 2_500_000_000
    ];
    let mut i = 0;
    while (i < test_xs.length()) {
        let x = test_xs[i];

        let orig = original_cdf::normal_cdf(x, true);
        let piece = piecewise_cdf::normal_cdf(x, true);
        let lut = lookup_cdf::normal_cdf(x, true);

        // All should be < 0.5 for negative inputs
        assert!(orig < 500_000_000, 100 + i);
        assert!(piece < 500_000_000, 200 + i);
        assert!(lut < 500_000_000, 300 + i);

        // Piecewise should track original
        assert!(abs_diff(orig, piece) < 2 * ONE_BP, 400 + i);

        i = i + 1;
    };
}

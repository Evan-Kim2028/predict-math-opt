/// Lookup table CDF with linear interpolation.
///
/// Alternative approach: 65 precomputed CDF values at spacing=0.0625 covering [0, 4].
/// Linear interpolation between entries.
///
/// Max interpolation error: 1.18 basis points (worse than piecewise cubic's 0.17 bp).
/// The binary-search if-else tree for 65 entries adds significant branch overhead,
/// making this approach slower than piecewise cubic in practice despite fewer arithmetic ops.
/// Included for comparison purposes.
#[allow(unused_const)]
module predict_math_opt::lookup_cdf;

const FLOAT_SCALING: u64 = 1_000_000_000;
const TABLE_STEP: u64 = 62_500_000;  // 0.0625 in FLOAT_SCALING
const TABLE_MAX_X: u64 = 4_000_000_000;
const TABLE_ENTRIES: u64 = 64;

/// Standard normal CDF using lookup table + linear interpolation.
/// 65 entries covering [0, 4] at spacing 0.0625.
/// Max error: ~1.18 basis points.
public fun normal_cdf(x: u64, x_negative: bool): u64 {
    if (x >= TABLE_MAX_X) {
        return if (x_negative) { 0 } else { FLOAT_SCALING }
    };

    let idx = x / TABLE_STEP;
    let frac = x - idx * TABLE_STEP; // remainder within interval

    let lo = table_value(idx);
    let hi = table_value(idx + 1);

    // Linear interpolation: lo + (hi - lo) * frac / TABLE_STEP
    let result = if (hi >= lo) {
        lo + ((hi - lo) as u128 * (frac as u128) / (TABLE_STEP as u128) as u64)
    } else {
        // Shouldn't happen (CDF is monotonic) but handle gracefully
        lo
    };

    if (x_negative) { FLOAT_SCALING - result } else { result }
}

/// Lookup table: 65 CDF values at x = 0, 0.0625, 0.125, ..., 4.0
/// Using if-else chain since Move const vectors have allocation overhead.
/// Binary search pattern: ~6 comparisons for 65 entries.
fun table_value(idx: u64): u64 {
    // Split into groups of 16 for binary search
    if (idx < 16) {
        if (idx < 8) {
            if (idx < 4) {
                if (idx == 0) { 500_000_001 }
                else if (idx == 1) { 524_917_740 }
                else if (idx == 2) { 549_738_266 }
                else { 574_365_676 } // 3
            } else {
                if (idx == 4) { 598_706_274 }
                else if (idx == 5) { 622_669_654 }
                else if (idx == 6) { 646_169_713 }
                else { 669_125_585 } // 7
            }
        } else {
            if (idx < 12) {
                if (idx == 8) { 691_462_468 }
                else if (idx == 9) { 713_112_336 }
                else if (idx == 10) { 734_014_532 }
                else { 754_116_223 } // 11
            } else {
                if (idx == 12) { 773_372_721 }
                else if (idx == 13) { 791_747_669 }
                else if (idx == 14) { 809_213_090 }
                else { 825_749_307 } // 15
            }
        }
    } else if (idx < 32) {
        if (idx < 24) {
            if (idx < 20) {
                if (idx == 16) { 841_344_740 }
                else if (idx == 17) { 855_995_592 }
                else if (idx == 18) { 869_705_436 }
                else { 882_484_712 } // 19
            } else {
                if (idx == 20) { 894_350_161 }
                else if (idx == 21) { 905_324_193 }
                else if (idx == 22) { 915_434_221 }
                else { 924_711_970 } // 23
            }
        } else {
            if (idx < 28) {
                if (idx == 24) { 933_192_771 }
                else if (idx == 25) { 940_914_868 }
                else if (idx == 26) { 947_918_730 }
                else { 954_246_403 } // 27
            } else {
                if (idx == 28) { 959_940_886 }
                else if (idx == 29) { 965_045_569 }
                else if (idx == 30) { 969_603_703 }
                else { 973_657_943 } // 31
            }
        }
    } else if (idx < 48) {
        if (idx < 40) {
            if (idx < 36) {
                if (idx == 32) { 977_249_938 }
                else if (idx == 33) { 980_419_988 }
                else if (idx == 34) { 983_206_754 }
                else { 985_647_029 } // 35
            } else {
                if (idx == 36) { 987_775_567 }
                else if (idx == 37) { 989_624_954 }
                else if (idx == 38) { 991_225_538 }
                else { 992_605_392 } // 39
            }
        } else {
            if (idx < 44) {
                if (idx == 40) { 993_790_320 }
                else if (idx == 41) { 994_803_894 }
                else if (idx == 42) { 995_667_514 }
                else { 996_400_497 } // 43
            } else {
                if (idx == 44) { 997_020_181 }
                else if (idx == 45) { 997_542_037 }
                else if (idx == 46) { 997_979_797 }
                else { 998_345_581 } // 47
            }
        }
    } else {
        if (idx < 56) {
            if (idx < 52) {
                if (idx == 48) { 998_650_033 }
                else if (idx == 49) { 998_902_449 }
                else if (idx == 50) { 999_110_907 }
                else { 999_282_393 } // 51
            } else {
                if (idx == 52) { 999_422_914 }
                else if (idx == 53) { 999_537_612 }
                else if (idx == 54) { 999_630_868 }
                else { 999_706_396 } // 55
            }
        } else {
            if (idx < 60) {
                if (idx == 56) { 999_767_327 }
                else if (idx == 57) { 999_816_290 }
                else if (idx == 58) { 999_855_484 }
                else { 999_886_735 } // 59
            } else {
                if (idx == 60) { 999_911_555 }
                else if (idx == 61) { 999_931_192 }
                else if (idx == 62) { 999_946_667 }
                else if (idx == 63) { 999_958_815 }
                else { 999_968_314 } // 64
            }
        }
    }
}

// === Tests ===

#[test]
fun test_cdf_at_zero() {
    let result = normal_cdf(0, false);
    let diff = if (result > 500_000_000) { result - 500_000_000 } else { 500_000_000 - result };
    assert!(diff < 100, 0);
}

#[test]
fun test_cdf_at_one() {
    // x=1.0 is at index 16, should be exact table hit
    let result = normal_cdf(1_000_000_000, false);
    let expected = 841_344_740;
    let diff = if (result > expected) { result - expected } else { expected - result };
    assert!(diff < 1_000, 1); // Exact table lookup, should be very close
}

#[test]
fun test_cdf_interpolated() {
    // x=0.03125 (halfway between idx 0 and 1)
    let result = normal_cdf(31_250_000, false);
    // Should be roughly midpoint of Phi(0) and Phi(0.0625)
    assert!(result > 500_000_000);
    assert!(result < 530_000_000);
}

#[test]
fun test_cdf_at_two() {
    let result = normal_cdf(2_000_000_000, false);
    let expected = 977_249_938;
    let diff = if (result > expected) { result - expected } else { expected - result };
    assert!(diff < 1_000, 2);
}

#[test]
fun test_cdf_symmetry() {
    let pos = normal_cdf(1_000_000_000, false);
    let neg = normal_cdf(1_000_000_000, true);
    let sum = pos + neg;
    let diff = if (sum > FLOAT_SCALING) { sum - FLOAT_SCALING } else { FLOAT_SCALING - sum };
    assert!(diff < 1_000, 3);
}

#[test]
fun test_cdf_monotonic() {
    let v0 = normal_cdf(0, false);
    let v1 = normal_cdf(250_000_000, false);
    let v2 = normal_cdf(500_000_000, false);
    let v3 = normal_cdf(1_000_000_000, false);
    let v4 = normal_cdf(2_000_000_000, false);
    assert!(v1 > v0);
    assert!(v2 > v1);
    assert!(v3 > v2);
    assert!(v4 > v3);
}

#[test]
fun test_cdf_tail() {
    assert!(normal_cdf(5_000_000_000, false) == FLOAT_SCALING);
    assert!(normal_cdf(5_000_000_000, true) == 0);
}

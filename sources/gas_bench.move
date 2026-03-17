/// Gas benchmarking: call each CDF implementation in a loop
/// to measure relative gas cost.
///
/// Run with: sui move test --gas-report gas_bench
#[test_only]
module predict_math_opt::gas_bench;

use predict_math_opt::original_cdf;
use predict_math_opt::piecewise_cdf;
use predict_math_opt::lookup_cdf;

/// Baseline: 10 calls to original A&S CDF at various x values.
#[test]
fun bench_original_10_calls() {
    let xs: vector<u64> = vector[
        100_000_000, 300_000_000, 500_000_000, 700_000_000, 900_000_000,
        1_100_000_000, 1_500_000_000, 2_000_000_000, 2_500_000_000, 3_000_000_000
    ];
    let mut i = 0;
    let mut sum = 0u64;
    while (i < 10) {
        sum = sum + original_cdf::normal_cdf(xs[i], false);
        i = i + 1;
    };
    assert!(sum > 0);
}

/// 10 calls to piecewise cubic CDF.
#[test]
fun bench_piecewise_10_calls() {
    let xs: vector<u64> = vector[
        100_000_000, 300_000_000, 500_000_000, 700_000_000, 900_000_000,
        1_100_000_000, 1_500_000_000, 2_000_000_000, 2_500_000_000, 3_000_000_000
    ];
    let mut i = 0;
    let mut sum = 0u64;
    while (i < 10) {
        sum = sum + piecewise_cdf::normal_cdf(xs[i], false);
        i = i + 1;
    };
    assert!(sum > 0);
}

/// 10 calls to lookup table CDF.
#[test]
fun bench_lookup_10_calls() {
    let xs: vector<u64> = vector[
        100_000_000, 300_000_000, 500_000_000, 700_000_000, 900_000_000,
        1_100_000_000, 1_500_000_000, 2_000_000_000, 2_500_000_000, 3_000_000_000
    ];
    let mut i = 0;
    let mut sum = 0u64;
    while (i < 10) {
        sum = sum + lookup_cdf::normal_cdf(xs[i], false);
        i = i + 1;
    };
    assert!(sum > 0);
}

/// Single call benchmarks for cleaner comparison.
#[test]
fun bench_original_single() {
    let result = original_cdf::normal_cdf(1_500_000_000, false);
    assert!(result > 900_000_000);
}

#[test]
fun bench_piecewise_single() {
    let result = piecewise_cdf::normal_cdf(1_500_000_000, false);
    assert!(result > 900_000_000);
}

#[test]
fun bench_lookup_single() {
    let result = lookup_cdf::normal_cdf(1_500_000_000, false);
    assert!(result > 900_000_000);
}

/// Simulate vault rebalance: 100 CDF calls with varied inputs.
#[test]
fun bench_original_100_calls() {
    let mut x: u64 = 50_000_000;
    let mut sum = 0u64;
    let mut i = 0;
    while (i < 100) {
        sum = sum + original_cdf::normal_cdf(x, i % 3 == 0);
        x = x + 35_000_000;
        if (x > 3_500_000_000) { x = x - 3_400_000_000; };
        i = i + 1;
    };
    assert!(sum > 0);
}

#[test]
fun bench_piecewise_100_calls() {
    let mut x: u64 = 50_000_000;
    let mut sum = 0u64;
    let mut i = 0;
    while (i < 100) {
        sum = sum + piecewise_cdf::normal_cdf(x, i % 3 == 0);
        x = x + 35_000_000;
        if (x > 3_500_000_000) { x = x - 3_400_000_000; };
        i = i + 1;
    };
    assert!(sum > 0);
}

#[test]
fun bench_lookup_100_calls() {
    let mut x: u64 = 50_000_000;
    let mut sum = 0u64;
    let mut i = 0;
    while (i < 100) {
        sum = sum + lookup_cdf::normal_cdf(x, i % 3 == 0);
        x = x + 35_000_000;
        if (x > 3_500_000_000) { x = x - 3_400_000_000; };
        i = i + 1;
    };
    assert!(sum > 0);
}

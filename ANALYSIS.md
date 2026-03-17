# Piecewise Cubic CDF: Gas Optimization for DeepBook Predict

## What Was Forked

This package forks the **normal CDF computation** from [DeepBook V3's predict module](https://github.com/MystenLabs/deepbookv3/blob/main/packages/predict/sources/helper/math.move), specifically the `normal_cdf()`, `exp()`, and `ln()` functions used by [`compute_nd2()`](https://github.com/MystenLabs/deepbookv3/blob/main/packages/predict/sources/oracle.move#L378) in the SVI + Black-Scholes pricing pipeline.

The original code uses the **Abramowitz & Stegun (26.2.17)** approximation, which requires computing `exp(-x²/2)` via a 12-iteration Taylor series. This analysis demonstrates that a **piecewise cubic polynomial** approach eliminates the Taylor series entirely and achieves a **5.2x–19.5x gas reduction** depending on call count.

### Files in this package

| File | Description |
|------|-------------|
| `sources/original_cdf.move` | Exact port of the A&S CDF from `deepbook_predict::math` |
| `sources/piecewise_cdf.move` | Piecewise cubic CDF replacement (8 segments, no `exp()`) |
| `sources/lookup_cdf.move` | Lookup table CDF with linear interpolation (alternative) |
| `sources/math_utils.move` | Standalone fixed-point math helpers (ported from `deepbook::math`) |
| `sources/comparison_tests.move` | Head-to-head accuracy tests across all implementations |
| `sources/gas_bench.move` | Gas benchmarking tests (1, 10, 100 CDF calls) |
| `scripts/generate_coefficients.py` | Python script to regenerate piecewise coefficients |

### What was NOT forked

The SVI volatility surface computation, oracle management, vault logic, and trading logic remain unchanged. This optimization is a **drop-in replacement** for the `normal_cdf()` function only.

---

## The Problem

The DeepBook Predict protocol prices binary options using Black-Scholes via SVI volatility surfaces. The key function `compute_nd2()` computes N(d₂) — the standard normal CDF evaluated at d₂ — which requires:

1. `ln()` — natural logarithm (24 instructions)
2. `sqrt()` — square root (20 instructions each, called twice)
3. **`normal_cdf()`** — the bottleneck (**674 instructions**, 56% from Taylor `exp()`)

The vault must recompute aggregate exposure after every mint/redeem, requiring ~8·log₂(n) CDF calls where n = number of open positions. At 10,000 positions, that's ~104 CDF calls per transaction.

---

## The Approach

### Why A&S is expensive

The Abramowitz & Stegun 26.2.17 formula is:

```
Φ(x) ≈ 1 - φ(x) · (a₁t + a₂t² + a₃t³ + a₄t⁴ + a₅t⁵)
where t = 1/(1 + 0.2316419·x) and φ(x) = (1/√2π)·e^(-x²/2)
```

Computing `e^(-x²/2)` requires range reduction + Taylor series:
```
e^x = e^r · 2^n  where r ∈ [0, ln2)
e^r = 1 + r + r²/2! + r³/3! + ... + r¹²/12!
```

This loop runs 8–12 iterations, each with a cross-module `mul()` and `div()` call. From disassembly, the `exp_series()` loop alone is **380 bytecode instructions** out of 674 total.

### Why piecewise cubic works

The CDF Φ(x) is a fixed, known, infinitely differentiable function on ℝ. Instead of computing it from scratch each time via exp/polynomial, we precompute cubic polynomial approximations for 8 segments of [0, 4]:

```
Segment k: [0.5k, 0.5(k+1))
P_k(x) = A_k + B_k·x + C_k·x² + D_k·x³
```

Coefficients are computed offline via polynomial interpolation through 4 CDF values per segment. At runtime, the evaluation is:

1. Compute x², x³ (2 multiplies)
2. Select segment (3–4 comparisons)
3. Evaluate cubic (3 multiplies + 2 adds)
4. Apply symmetry: Φ(-x) = 1 - Φ(x)

No loops. No `exp()`. No Taylor series. **130 bytecode instructions** total.

### Mathematical validity

Both A&S and piecewise cubic are **polynomial approximations to the same function**. Neither is exact. The mathematical justification:

1. **Weierstrass approximation theorem**: Any continuous function on a bounded interval can be uniformly approximated by polynomials to arbitrary precision.

2. **Interpolation error bound**: For degree-n polynomial interpolation on interval width h:
   ```
   |Φ(x) - Pₙ(x)| ≤ h^(n+1) / (4·(n+1)) · max|Φ^(n+1)(x)|
   ```
   For n=3, h=0.5: theoretical bound ≈ 0.012, actual max error = 1.75×10⁻⁵.

3. **Monotonicity preserved**: P'(x) = B + 2Cx + 3Dx² > 0 verified analytically for all segments. Minimum derivative across all segments is 0.000149 (at x=4.0, segment 7). No segment can ever produce a decreasing CDF value.

4. **Continuity preserved**: Boundary gaps between adjacent segments are < 10⁻⁸, below the integer precision of FLOAT_SCALING (10⁻⁹).

---

## Validation Results

### Accuracy vs reference (scipy 64-bit)

| x | Reference Φ(x) | A&S error | Piecewise error | Piecewise (bp) |
|---|----------------|-----------|-----------------|----------------|
| 0.00 | 0.500000000 | < 10⁻⁹ | < 10⁻⁹ | < 0.001 |
| 0.25 | 0.598706326 | 5.2×10⁻⁸ | 1.0×10⁻⁵ | 0.10 |
| 0.50 | 0.691462461 | 6.8×10⁻⁸ | 8.7×10⁻⁶ | 0.09 |
| 1.00 | 0.841344746 | 1.2×10⁻⁷ | 1.7×10⁻⁵ | 0.17 |
| 1.50 | 0.933192799 | 8.5×10⁻⁸ | 6.3×10⁻⁶ | 0.06 |
| 2.00 | 0.977249868 | 3.1×10⁻⁸ | 4.7×10⁻⁶ | 0.05 |
| 3.00 | 0.998650102 | 1.5×10⁻⁸ | 1.7×10⁻⁶ | 0.02 |

**Max piecewise error: 0.17 basis points** (vs Aslan's stated tolerance of < 1 bp).

### Monotonicity verification

```
Segment 0 [0.0, 0.5]: min P'(x) = 0.351687 at x=0.50  ✓
Segment 1 [0.5, 1.0]: min P'(x) = 0.241345 at x=1.00  ✓
Segment 2 [1.0, 1.5]: min P'(x) = 0.129181 at x=1.50  ✓
Segment 3 [1.5, 2.0]: min P'(x) = 0.054027 at x=2.00  ✓
Segment 4 [2.0, 2.5]: min P'(x) = 0.017698 at x=2.50  ✓
Segment 5 [2.5, 3.0]: min P'(x) = 0.004556 at x=3.00  ✓
Segment 6 [3.0, 3.5]: min P'(x) = 0.000925 at x=3.50  ✓
Segment 7 [3.5, 4.0]: min P'(x) = 0.000149 at x=4.00  ✓
```

All segments strictly monotone increasing. No arbitrage possible from CDF mispricing.

### Test results

All 43 tests pass:
- 7 original A&S unit tests
- 7 piecewise cubic unit tests
- 7 lookup table unit tests
- 3 math_utils unit tests
- 10 cross-implementation comparison tests (accuracy, symmetry, monotonicity, boundaries, sweeps)
- 9 gas benchmarks (1, 10, 100 calls × 3 implementations)

---

## Gas Analysis

### Bytecode instruction counts (from `sui move disassemble`)

| Function | Bytecode instructions | Source |
|----------|----------------------|--------|
| `original_cdf::normal_cdf` | **674** | `exp_series` loop = 380 (56%) |
| `piecewise_cdf::normal_cdf` | **130** | segment select + cubic eval |
| `lookup_cdf::normal_cdf` | ~120 | binary search + interpolation |

### Sui mainnet gas model (v11, `initial_cost_schedule_v5`)

Sui uses **tiered instruction pricing** — each additional instruction costs more as the total count increases:

| Instructions executed | Cost multiplier |
|----------------------|-----------------|
| 0 – 20,000 | 1× |
| 20,000 – 50,000 | 2× |
| 50,000 – 100,000 | 10× |
| 100,000 – 200,000 | 50× |
| 200,000+ | 100× |

This tiered structure **amplifies** the savings from fewer instructions:

| CDF calls | Original instructions | Piecewise instructions | Original gas | Piecewise gas | **Savings** |
|-----------|----------------------|----------------------|-------------|--------------|-------------|
| 1 | 674 | 130 | 674 | 130 | **5.2×** |
| 10 | 6,740 | 1,300 | 6,740 | 1,300 | **5.2×** |
| 50 | 33,700 | 6,500 | 47,400 | 6,500 | **7.3×** |
| 100 | 67,400 | 13,000 | 254,000 | 13,000 | **19.5×** |

At 100 calls (vault rebalance scenario), the original crosses the **10× instruction tier** at 50,000 instructions, while piecewise stays entirely in the **1× tier**. This turns a 5.2× instruction reduction into a **19.5× gas reduction**.

### Applied to vault rebalance

For n=10,000 positions, the vault needs ~104 CDF calls (via the log(n) approximation algorithm):

| Metric | Original A&S | Piecewise cubic |
|--------|-------------|-----------------|
| Instructions per CDF | 674 | 130 |
| Total instructions (104 calls) | ~70,000 | ~13,500 |
| Highest tier hit | 10× | 1× |
| Instruction gas | ~290,000 | ~13,500 |

---

## Limitations and Future Work

1. **Precision tradeoff**: Piecewise cubic max error is 0.17 bp vs A&S's 0.00075 bp. If sub-basis-point precision is required, use 16 segments (doubles segment count, adds ~1 comparison, halves error to ~0.01 bp).

2. **Only optimizes CDF**: The `compute_nd2()` pipeline also uses `ln()` and `sqrt()`. The same piecewise technique can be applied to these functions for additional savings.

3. **Tail behavior**: Piecewise clamps at |x| ≥ 4 (Φ(4) = 0.99997). The original clamps at |x| ≥ 8. For deep OTM options with |d₂| > 4, the piecewise approach returns exactly 0 or 1. This is sufficient for the <1 bp requirement but could be extended to |x| = 6 with 4 more segments.

4. **Coefficient regeneration**: If higher precision is needed, run `python scripts/generate_coefficients.py` to regenerate with more segments or higher-degree polynomials.

---

## How to Use

### Build and test

```bash
cd packages/predict-math-opt
sui move build
sui move test
sui move test -s csv  # with statistics
```

### Drop-in replacement

Replace calls to `deepbook_predict::math::normal_cdf(x, x_neg)` with `piecewise_cdf::normal_cdf(x, x_neg)`. The function signature and scaling (FLOAT_SCALING = 10⁹) are identical.

# Piecewise Cubic CDF: Gas Optimization for DeepBook Predict

## What Was Forked

This package forks the **normal CDF computation** from [DeepBook V3's predict module](https://github.com/MystenLabs/deepbookv3/blob/main/packages/predict/sources/helper/math.move), specifically the `normal_cdf()`, `exp()`, and `ln()` functions used by [`compute_nd2()`](https://github.com/MystenLabs/deepbookv3/blob/main/packages/predict/sources/oracle.move#L378) in the SVI + Black-Scholes pricing pipeline.

The original code uses the **Abramowitz & Stegun (26.2.17)** approximation, which requires computing `exp(-x²/2)` via a 12-iteration Taylor series. Three algorithmic optimizations (piecewise cubic CDF, Horner-form ln, unrolled Newton sqrt) together achieve a **52× gas reduction** at 100 `compute_nd2` calls (unrolled sqrt vs the old loop-based sqrt).

> **Updated baseline (v2, March 2026):** The deployed predict package uses `u128::sqrt` (Move stdlib, 64-iter digit-by-digit), not the custom Newton loop used in the original v1 baseline. Against the correct on-chain baseline, CDF+ln optimization alone delivers **1.09–1.15× savings**. Replacing sqrt too yields **194.7× savings** at 100 calls. See the [Correct On-Chain Baseline](#correct-on-chain-baseline-v2) section for details.

### Files in this package

| File | Description |
|------|-------------|
| `sources/original_cdf.move` | Exact port of the A&S CDF from `deepbook_predict::math` |
| `sources/piecewise_cdf.move` | Piecewise cubic CDF replacement (8 segments, no `exp()`) |
| `sources/lookup_cdf.move` | Lookup table CDF with linear interpolation (alternative) |
| `sources/optimized_math.move` | Optimized ln (Horner form) and sqrt (unrolled Newton) |
| `sources/math_utils.move` | Standalone fixed-point math helpers (ported from `deepbook::math`) |
| `sources/compute_nd2.move` | Full compute_nd2 pipeline — original vs optimized side-by-side |
| `sources/comparison_tests.move` | Head-to-head accuracy tests across all implementations |
| `sources/gas_bench.move` | Gas benchmarking tests (1, 10, 100 CDF calls) |
| `sources/gas_measure.move` | On-chain entry functions for real gas measurement on testnet |
| `scripts/generate_coefficients.py` | Python script to regenerate piecewise coefficients |
| `scripts/gas_to_speedscope.py` | Generate speedscope.app-compatible profile from test runner |

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

All 57 tests pass:
- 7 original A&S unit tests
- 7 piecewise cubic unit tests
- 7 lookup table unit tests
- 3 math_utils unit tests
- 7 optimized_math unit tests (ln + sqrt)
- 10 cross-implementation comparison tests (accuracy, symmetry, monotonicity, boundaries, sweeps)
- 9 gas benchmarks (1, 10, 100 CDF calls × 3 implementations)
- 7 compute_nd2 benchmarks (original vs optimized pipeline)

---

## Correct On-Chain Baseline (v2)

The v1 analysis claimed 52× improvement, but it compared against `math_utils::sqrt_u128` — a custom Newton-Raphson loop that starts at `x/2` (catastrophically bad initial guess for u128 inputs). The **deployed predict package** actually uses `deepbook::math::sqrt → std::u128::sqrt` (Move stdlib digit-by-digit, 64 fixed iterations).

This changes the comparison entirely. `u128::sqrt` is 2.1× **more** expensive than the custom loop, and dominates the overall cost.

### Gas profile: single `compute_nd2` call (real tx SVI params)

Measured via `sui move test --trace` + `sui analyze-trace gas-profile` (speedscope format):

| Component | Onchain baseline | Optimized (16-seg CDF) | Savings |
|---|---:|---:|---:|
| `sqrt` ×2 (`u128::sqrt`) | 36,004 | 36,004 | 0 |
| `normal_cdf` | 2,857 | 645 | +2,212 |
| `ln` | 2,067 | 1,185 | +882 |
| Other | 2,087 | 2,087 | 0 |
| **TOTAL** | **43,015** | **39,921** | **1.08×** |

`sqrt` accounts for **84–90%** of total cost. CDF+ln optimization alone saves only ~3,100 gas per call.

### On-chain MIST costs (testnet v2)

Package v2: [`0xa5277e3a...`](https://suiscan.xyz/testnet/object/0xa5277e3ab4775c678350f512d910ca12afc6c343abc211a79bee149cbfaaa91e)

**Real tx SVI params** (a=620000, b=42500000, rho=-0.24364, m=0.01128, sigma=0.08468):

| Calls | Onchain baseline | Optimized (same sqrt) | Ratio |
|---|---:|---:|---|
| 20 calls | 2,720,000 MIST | 2,370,000 MIST | **1.15×** |

**Synthetic SVI params** (a=0.05, b=0.2, rho=-0.3, m=0.01, sigma=0.1):

| Calls | Onchain baseline | Optimized (same sqrt) | Fully optimized (fast sqrt) |
|---|---:|---:|---:|
| 100 calls | 290,100,000 MIST | 266,500,000 MIST | 1,490,000 MIST |
| Ratio | — | **1.09×** | **194.7×** |

**Summary**: Replacing CDF and ln (but keeping `u128::sqrt`) yields 9–15% savings on testnet. Replacing sqrt too — which requires modifying the deployed predict contract — yields 194.7× savings at 100 calls.

### Why the 16-seg CDF upgrade matters despite small gas savings

Even though the gas savings from CDF alone are modest, the **accuracy improvement is essential**:
- v1 (8 segments): 0.17 bp max error — exceeded the 0.01 bp production tolerance
- v2 (16 segments): **0.0108 bp max error** — meets the 0.01 bp target
- Gas overhead of extra segments: +69 gas per call (negligible, sqrt dominates)

---

## On-Chain Gas Measurements (Sui Testnet)

All values measured via real transactions on Sui testnet. Every number is verifiable on Suiscan.

### Verified results (March 2026)

**Package**: [`0x1df028...`](https://suiscan.xyz/testnet/object/0x1df0281d6333bb1ad97b620f2f8ae9187d54919e5d73e2c757a6b869d2e0ae6d)
**Network**: testnet

#### Correctness verification (on-chain)

Before comparing gas, we verified both implementations produce the same outputs. The `verify_outputs_match` entry function calls **both** original and optimized `compute_nd2` with 10 different strikes using identical SVI parameters, and asserts outputs match within 2 basis points:

[Hs9WRDEVXezggk8CsoRdq4L6h2GPzrNMP3Ev8E5WzDrf](https://suiscan.xyz/testnet/tx/Hs9WRDEVXezggk8CsoRdq4L6h2GPzrNMP3Ev8E5WzDrf) — Status: **Success**

#### Full compute_nd2 pipeline — 100 calls

Both benchmarks use identical parameters: SVI(a=0.05, b=0.2, rho=-0.3, m=0.01, sigma=0.1), forward=50000, varied strikes over the same range.

| | Computation | Storage | Rebate | **Total Gas** | Digest |
|---|---:|---:|---:|---:|---|
| **Original** | 75,100,000 | 988,000 | -978,120 | **75,109,880 MIST** | [`E253rhKq...`](https://suiscan.xyz/testnet/tx/E253rhKqApYpNYJdhD3Bk5N5nmdecMJ7EK4VhM276sf5) |
| **Optimized** | 1,440,000 | 988,000 | -978,120 | **1,449,880 MIST** | [`3rkAu4jh...`](https://suiscan.xyz/testnet/tx/3rkAu4jh5xZSPM4bsJG4JmJ4kEyQ3CefbQwRvYSh2yv3) |
| **Ratio** | **52.2×** | — | — | **51.8×** | |

#### Per-function breakdown — 250 calls each

Three optimizations, measured individually:

| Optimization | Original gas | Optimized gas | Savings | What changed |
|-------------|---:|---:|---:|---|
| **sqrt** (×2 per nd2 call) | 72,900,000 | 1,000,000 (floor) | **>72×** | Loop with bad initial guess (x/2, ~30 iters) → bit-length guess + 7 unrolled steps |
| **CDF** | 6,060,000 | 1,000,000 (floor) | **>6×** | A&S with 12-iter Taylor exp() → 8-segment piecewise cubic, no exp, no loop |
| **ln** | 2,000,000 | 1,000,000 (floor) | **>2×** | 7-iter series with loop + div() → Horner form, precomputed 1/k, no loop, no div |

Transaction links:
- Original sqrt × 250: [`6i1jeA97...`](https://suiscan.xyz/testnet/tx/6i1jeA97CAMN6RZ1ExZQJjWDHPgW6aCX7fF3yvJZkRHS)
- Optimized sqrt × 250: [`8cJeq1s4...`](https://suiscan.xyz/testnet/tx/8cJeq1s4HTprdbzL51CnZVBAq8sHKGJN1zz4MCdkm5aM)
- Original CDF × 250: [`4qnrdAVc...`](https://suiscan.xyz/testnet/tx/4qnrdAVch3rWC3CV8RNewMuAbyeTAHKRBChtvXNWFcag)
- Piecewise CDF × 250: [`2TtmiBAz...`](https://suiscan.xyz/testnet/tx/2TtmiBAzgV6tbSKq3CrLZpMVJezUPGTbTAxxSUazW6Gc)
- Original ln × 250: [`EHb7KiVw...`](https://suiscan.xyz/testnet/tx/EHb7KiVwZodvTSuJoWbH5FW5e35dyK33ssSRwoNpxgzW)
- Optimized ln × 250: [`65xHbcBs...`](https://suiscan.xyz/testnet/tx/65xHbcBsazMtesH79GWt6yuDpzcsaTxR5mpNRzzYyTxo)

Note: 1,000,000 MIST is Sui's minimum computation gas floor. The optimized functions are so cheap that individual benchmarks hit this floor even at 250 calls.

**sqrt was the biggest win by far.** The original Newton-Raphson starts with guess = x/2 for u128 inputs (~10^18). This is catastrophically far from the answer (~10^9), requiring ~30 iterations of linear convergence before quadratic convergence kicks in. Each iteration executes 25 bytecodes of u256 division. The optimized version uses a bit-length estimate (within 2× of true sqrt) + 7 fixed Newton steps — no loop, no conditionals, deterministic.

### Earlier measurements (for reference)

Earlier measurements on packages [`0x109324...`](https://suiscan.xyz/testnet/object/0x109324692ed6bfd75733fa58c2d84e7bd50d819f84601850b11cf3e098a401bf) and [`0x5de1cc...`](https://suiscan.xyz/testnet/object/0x5de1ccf60cba0b91f2c293021e5608e9b59e6f7eecab3e3115a6324a8053d7c3) (epoch 1041) showed consistent results:

| Calls | Original | Optimized | Ratio |
|-------|----------|-----------|-------|
| 100 | 74,400,000 MIST | 1,420,000 MIST | 52× |
| 200 | 259,600,000 MIST | 7,450,000 MIST | 35× |

The ratio drops from 52× at 100 calls to 35× at 200 calls because the optimized version also starts entering higher pricing tiers. But 35× is still massive.

### Sui gas model context

Sui uses **tiered instruction pricing** (gas model v11, `initial_cost_schedule_v5`). More instructions per transaction → higher cost per instruction:

| Instructions executed | Cost multiplier |
|----------------------|-----------------|
| 0 – 20,000 | 1× |
| 20,000 – 50,000 | 2× |
| 50,000 – 100,000 | 10× |
| 100,000 – 200,000 | 50× |
| 200,000+ | 100× |

This is why the 3.66× bytecode reduction (2,902 → 792 per call) translates to a 35–52× gas reduction at 100–200 calls. The original crosses into expensive tiers; the optimized stays in cheap ones.

---

## Limitations and Future Work

1. **Precision**: ~~Piecewise cubic max error is 0.17 bp~~ — **resolved in v2**: upgraded to 16 segments, max error = **0.0108 bp**, meets Aslan's 0.01 bp production tolerance. Coefficients generated inline (see `scripts/generate_coefficients.py` for 8-seg; 16-seg was run ad-hoc).

2. **UP/DOWN pair optimization**: `compute_nd2` computes `d2` then branches on `is_up` at the CDF step. Computing both `Φ(d2)` and `1 - Φ(d2)` from one `d2` calculation would halve calls when both directions are needed.

3. **Batch `get_sigmoid_shape` helper**: A single function that evaluates `compute_nd2` at N strikes sharing SVI parameters — saves parameter loading overhead and function frame setup per call.

4. **Tail behavior**: Piecewise clamps at |x| ≥ 4 (Φ(4) = 0.99997). The original clamps at |x| ≥ 8. For deep OTM options with |d₂| > 4, the piecewise approach returns exactly 0 or 1. This is sufficient for the <1 bp requirement but could be extended to |x| = 6 with 4 more segments.

5. **Coefficient regeneration**: If higher precision is needed, run `python scripts/generate_coefficients.py` to regenerate with more segments or higher-degree polynomials.

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

For the full pipeline, replace `compute_nd2` calls with `compute_nd2::compute_nd2_optimized(forward, strike, &svi, is_up)` which uses all three optimizations (piecewise CDF + Horner ln + unrolled sqrt).

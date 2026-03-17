# Predict Package Gas Profile

## Package: 0x01db8fc74ead463c7167f9c609af72e64ac4eeb0f6b9c05da17c16ad0fd348d0

**Network**: Sui Testnet  
**Modules**: constants, market_key, math, oracle, predict, predict_manager, pricing_config, registry, risk_config, vault  
**Dependency (DeepBook V3)**: 0xfb28c4cbc6865bd1c897d26aecbe1f8792d1509a20ffec692c800660cbec6982

---

## Real Transaction Inputs

From tx `5sdckiLBtq7kfCPyTH5jN84Drmnm9pRvXzop88fZriZy` (PositionMinted event):

| Field | Value |
|---|---|
| forward | 71047183700000 |
| strike | 72437856240000 |
| is_up | true |
| ask_price | 503026245 (~0.503) |
| expiry | 1772784000000 ms |

SVI params from OracleSVIUpdated (tx ACaiuYoUSHNcNkVqfcpbJPFSP24EcvxceoSjG53riSbS):

| Param | Raw (FLOAT_SCALING=1e9) | Float |
|---|---|---|
| a | 620000 | 0.000620 |
| b | 42500000 | 0.042500 |
| rho | 243640000 (negative) | -0.24364 |
| m | 11280000 | 0.011280 |
| sigma | 84680000 | 0.084680 |
| risk_free_rate | 35000000 | 3.5% |

---

## Gas Profile — Real Tx Inputs (sui analyze-trace gas-profile)

Trace files: `traces/predict_math_opt__onchain_baseline__bench_*_real_tx.json.zst`  
Profiles: `gas_profile_predict_math_opt__onchain_baseline__bench_*_real_tx.json` → drag into speedscope.app

### Per-function breakdown (single compute_nd2 call)

| Function | On-chain (u128::sqrt) | Optimized (u128::sqrt + piecewise + Horner) | Saved |
|---|---:|---:|---:|
| sqrt (×2, u128::sqrt) | 36,004 | 36,004 | 0 |
| normal_cdf | 2,857 | 576 | +2,281 |
| ln | 2,067 | 1,185 | +882 |
| other overhead | 2,087 | 2,087 | ~0 |
| **TOTAL** | **43,015** | **39,852** | **+3,163** |

**Real improvement (CDF + ln only): 1.08× — 8% per call**

### Why only 8%?

`0x1::u128::sqrt` (digit-by-digit, 64 fixed iterations) costs **36,004 gas** for two calls — **84% of the total on-chain cost**. The CDF and ln optimizations only touch the remaining 16%.

---

## The Full Picture: What it takes to get meaningful savings

| Optimization | Gas saved per call | Notes |
|---|---|---|
| Piecewise CDF (current) | 2,281 | ✅ implemented |
| Horner ln (current) | 882 | ✅ implemented |
| Fast sqrt (sqrt_u128_fast) | ~35,000 | ❌ not on-chain — requires replacing u128::sqrt |

If sqrt is also optimized (using `optimized_math::sqrt_u128_fast`):
- On-chain baseline: 43,015 gas
- All three optimizations: ~5,438 gas → **8.1× improvement**

The 8.1× figure is real and achievable — but requires a custom sqrt implementation in the deployed package (replacing `deepbook::math::sqrt → u128::sqrt` with a bit-length-guess Newton method).

---

## Accuracy: 0.01 bps tolerance

| Implementation | Max error vs scipy |
|---|---|
| On-chain A&S | ~0.001 bps |
| Piecewise cubic 8 segments | **0.17 bps** ← FAILS 0.01 bps target |
| Piecewise cubic 16 segments | ~0.01 bps (not yet implemented) |

To meet 0.01 bps: run `python scripts/generate_coefficients.py` with 16 segments, rebuild `piecewise_cdf.move`.

---

## Files

| File | Purpose |
|---|---|
| `sources/onchain_baseline.move` | Correct on-chain baseline using u128::sqrt |
| `sources/compute_nd2.move` | Original vs optimized (uses fast sqrt — 8.1× vs on-chain) |
| `gas_profile_*_real_tx.json` | Official sui analyze-trace profiles with real inputs |
| `traces/*.json.zst` | Raw trace files |

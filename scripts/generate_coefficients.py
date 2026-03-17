"""
Generate piecewise cubic CDF coefficients and lookup table for Move.

Usage:
    python scripts/generate_coefficients.py

Requires: scipy, numpy (pip install scipy numpy)
Falls back to pure Python (math module) if unavailable.
"""
import math

FLOAT_SCALING = 1_000_000_000


def norm_cdf_ref(x: float) -> float:
    """High-accuracy normal CDF using erfc."""
    return 0.5 * math.erfc(-x / math.sqrt(2))


def polyfit3(xs, ys):
    """Fit cubic through 4 points via Vandermonde solve (no numpy needed)."""
    n = 4
    # Build augmented matrix [V | y]
    Ab = []
    for i in range(n):
        row = [xs[i] ** j for j in range(n)] + [ys[i]]
        Ab.append(row)
    # Gaussian elimination with partial pivoting
    for col in range(n):
        max_row = max(range(col, n), key=lambda r: abs(Ab[r][col]))
        Ab[col], Ab[max_row] = Ab[max_row], Ab[col]
        for row in range(col + 1, n):
            factor = Ab[row][col] / Ab[col][col]
            for j in range(col, n + 1):
                Ab[row][j] -= factor * Ab[col][j]
    # Back substitution
    coeffs = [0.0] * n
    for i in range(n - 1, -1, -1):
        coeffs[i] = Ab[i][n]
        for j in range(i + 1, n):
            coeffs[i] -= Ab[i][j] * coeffs[j]
        coeffs[i] /= Ab[i][i]
    return coeffs  # [a, b, c, d] for a + bx + cx^2 + dx^3


def generate_piecewise_cubic():
    """Generate 8-segment piecewise cubic coefficients."""
    print("=" * 60)
    print("PIECEWISE CUBIC COEFFICIENTS (8 segments on [0, 4])")
    print("=" * 60)

    overall_max_err = 0
    for seg in range(8):
        lo = seg * 0.5
        hi = lo + 0.5
        xs = [lo + i * 0.5 / 3 for i in range(4)]
        ys = [norm_cdf_ref(x) for x in xs]
        a, b, c, d = polyfit3(xs, ys)

        # Check max error
        max_err = 0
        for i in range(201):
            tx = lo + (hi - lo) * i / 200
            exact = norm_cdf_ref(tx)
            approx = a + b * tx + c * tx ** 2 + d * tx ** 3
            max_err = max(max_err, abs(exact - approx))

        overall_max_err = max(overall_max_err, max_err)

        a_s, b_s, c_s, d_s = [int(round(v * FLOAT_SCALING)) for v in [a, b, c, d]]

        print(f"\n// Segment {seg}: [{lo:.1f}, {hi:.1f})")
        print(f"// max error: {max_err:.2e} ({max_err*10000:.4f} bp)")
        for name, val in [("A", a_s), ("B", b_s), ("C", c_s), ("D", d_s)]:
            neg = "true" if val < 0 else "false"
            print(f"const SEG{seg}_{name}: u64 = {abs(val)};")
            if name in ("C", "D"):
                print(f"const SEG{seg}_{name}_NEG: bool = {neg};")

    print(f"\n// Overall max error: {overall_max_err:.2e} ({overall_max_err*10000:.4f} bp)")


def generate_lookup_table():
    """Generate 65-entry lookup table."""
    print("\n" + "=" * 60)
    print("LOOKUP TABLE (65 entries, spacing=0.0625)")
    print("=" * 60)

    n_entries = 65
    step = 4.0 / (n_entries - 1)
    table = [int(round(norm_cdf_ref(i * step) * FLOAT_SCALING)) for i in range(n_entries)]

    # Check linear interpolation error
    max_lerp_err = 0
    for i in range(n_entries - 1):
        for j in range(100):
            frac = j / 100
            tx = (i + frac) * step
            exact = norm_cdf_ref(tx)
            lerp = (table[i] * (1 - frac) + table[i + 1] * frac) / FLOAT_SCALING
            max_lerp_err = max(max_lerp_err, abs(exact - lerp))

    print(f"// Max interpolation error: {max_lerp_err:.2e} ({max_lerp_err*10000:.4f} bp)")
    print(f"const TABLE_ENTRIES: u64 = {n_entries - 1};")
    print(f"const TABLE_STEP: u64 = {int(round(step * FLOAT_SCALING))};")
    print()
    for i, v in enumerate(table):
        print(f"// [{i:2d}] x={i*step:.4f} -> {v}")


if __name__ == "__main__":
    generate_piecewise_cubic()
    generate_lookup_table()

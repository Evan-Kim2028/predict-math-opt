#!/usr/bin/env python3
"""
Run Sui Move gas benchmarks and generate a speedscope-compatible JSON profile.

Usage:
    cd packages/predict-math-opt
    python scripts/gas_to_speedscope.py

Opens in browser:
    python scripts/gas_to_speedscope.py --open

Output: gas_profile.speedscope.json  (drag into https://www.speedscope.app/)
"""

import subprocess
import json
import csv
import io
import sys
import os
import re
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent.parent


def run_benchmarks() -> list[dict]:
    """Run sui move test -s csv and parse the output."""
    result = subprocess.run(
        ["sui", "move", "test", "-s", "csv", "--silence-warnings"],
        capture_output=True,
        text=True,
        cwd=PACKAGE_DIR,
    )
    output = result.stdout + result.stderr

    # Find the CSV section (starts with "name,nanos,gas")
    lines = output.strip().split("\n")
    csv_start = None
    for i, line in enumerate(lines):
        if line.startswith("name,nanos,gas"):
            csv_start = i
            break

    if csv_start is None:
        print("ERROR: Could not find CSV output. Raw output:")
        print(output)
        sys.exit(1)

    # Parse CSV rows until we hit a non-CSV line
    csv_lines = []
    for line in lines[csv_start:]:
        if "," in line and not line.startswith("Test result"):
            csv_lines.append(line)
        else:
            break

    reader = csv.DictReader(io.StringIO("\n".join(csv_lines)))
    return [row for row in reader]


def classify_test(name: str) -> tuple[str, str, str]:
    """Classify a test into (category, implementation, label).

    Returns (category, impl, short_label) where:
      - category groups tests for side-by-side comparison
      - impl is 'original' or 'optimized'
      - short_label is the human-readable name
    """
    # gas_bench tests: CDF implementations
    if "gas_bench" in name:
        m = re.search(r"bench_(original|piecewise|lookup)_(\w+)", name)
        if m:
            impl_raw, detail = m.group(1), m.group(2)
            impl = "original" if impl_raw == "original" else "optimized"
            variant = "piecewise" if impl_raw == "piecewise" else (
                "lookup" if impl_raw == "lookup" else "A&S"
            )
            return (f"CDF {detail}", impl, f"{variant} CDF ({detail})")

    # compute_nd2 benchmarks
    if "bench_compute_nd2" in name:
        m = re.search(r"bench_compute_nd2_(original|optimized)_(\w+)", name)
        if m:
            impl, detail = m.group(1), m.group(2)
            return (f"nd2 {detail}", impl, f"compute_nd2 ({detail})")

    # comparison_tests — group as "validation"
    if "comparison_tests" in name:
        short = name.split("::")[-1]
        return ("validation", "validation", short)

    # Per-module unit tests
    if "original_cdf::" in name:
        short = name.split("::")[-1]
        return ("unit_original", "original", short)
    if "piecewise_cdf::" in name:
        short = name.split("::")[-1]
        return ("unit_piecewise", "optimized", short)
    if "lookup_cdf::" in name:
        short = name.split("::")[-1]
        return ("unit_lookup", "optimized", short)
    if "optimized_math::" in name:
        short = name.split("::")[-1]
        return ("unit_opt_math", "optimized", short)
    if "math_utils::" in name:
        short = name.split("::")[-1]
        return ("unit_math_utils", "original", short)

    return ("other", "other", name.split("::")[-1])


def build_speedscope(rows: list[dict]) -> dict:
    """Build speedscope JSON with two profiles: Original vs Optimized.

    Uses the 'evented' profile format where the horizontal axis = gas cost,
    giving a flame-chart view of relative gas consumption.
    """
    frames = []
    frame_index = {}

    def get_frame(name: str, file: str = "") -> int:
        key = (name, file)
        if key not in frame_index:
            frame_index[key] = len(frames)
            f = {"name": name}
            if file:
                f["file"] = file
            frames.append(f)
        return frame_index[key]

    # Separate benchmarks into profiles
    bench_rows = []
    for row in rows:
        name = row["name"]
        nanos = int(row["nanos"])
        gas = int(row["gas"])
        category, impl, label = classify_test(name)
        bench_rows.append({
            "name": name,
            "nanos": nanos,
            "gas": gas,
            "category": category,
            "impl": impl,
            "label": label,
        })

    # --- Profile 1: Gas comparison (evented, unit=gas) ---
    # Show only benchmark tests side by side
    original_events = []
    optimized_events = []
    orig_cursor = 0
    opt_cursor = 0

    # Group by category for paired comparison
    from collections import defaultdict
    categories = defaultdict(list)
    for r in bench_rows:
        if r["category"].startswith("unit_") or r["category"] == "validation":
            continue
        categories[r["category"]].append(r)

    for cat in sorted(categories.keys()):
        items = categories[cat]
        for r in items:
            gas = r["gas"]
            # Category frame (parent)
            cat_frame = get_frame(r["category"])
            # Function frame (child)
            fn_frame = get_frame(r["label"], r["name"].replace("::", "/") + ".move")

            if r["impl"] == "original":
                original_events.append({"type": "O", "frame": cat_frame, "at": orig_cursor})
                original_events.append({"type": "O", "frame": fn_frame, "at": orig_cursor})
                orig_cursor += gas
                original_events.append({"type": "C", "frame": fn_frame, "at": orig_cursor})
                original_events.append({"type": "C", "frame": cat_frame, "at": orig_cursor})
            else:
                optimized_events.append({"type": "O", "frame": cat_frame, "at": opt_cursor})
                optimized_events.append({"type": "O", "frame": fn_frame, "at": opt_cursor})
                opt_cursor += gas
                optimized_events.append({"type": "C", "frame": fn_frame, "at": opt_cursor})
                optimized_events.append({"type": "C", "frame": cat_frame, "at": opt_cursor})

    # --- Profile 2: Wall-time comparison (nanoseconds) ---
    orig_time_events = []
    opt_time_events = []
    orig_t = 0
    opt_t = 0

    for cat in sorted(categories.keys()):
        items = categories[cat]
        for r in items:
            nanos = r["nanos"]
            cat_frame = get_frame(r["category"])
            fn_frame = get_frame(r["label"], r["name"].replace("::", "/") + ".move")

            if r["impl"] == "original":
                orig_time_events.append({"type": "O", "frame": cat_frame, "at": orig_t})
                orig_time_events.append({"type": "O", "frame": fn_frame, "at": orig_t})
                orig_t += nanos
                orig_time_events.append({"type": "C", "frame": fn_frame, "at": orig_t})
                orig_time_events.append({"type": "C", "frame": cat_frame, "at": orig_t})
            else:
                opt_time_events.append({"type": "O", "frame": cat_frame, "at": opt_t})
                opt_time_events.append({"type": "O", "frame": fn_frame, "at": opt_t})
                opt_t += nanos
                opt_time_events.append({"type": "C", "frame": fn_frame, "at": opt_t})
                opt_time_events.append({"type": "C", "frame": cat_frame, "at": opt_t})

    # --- Profile 3: All tests by gas (single view) ---
    all_events = []
    cursor = 0
    for r in sorted(bench_rows, key=lambda x: (x["impl"], x["category"])):
        gas = r["gas"]
        impl_frame = get_frame(f"[{r['impl']}]")
        cat_frame = get_frame(r["category"])
        fn_frame = get_frame(r["label"], r["name"].replace("::", "/") + ".move")

        all_events.append({"type": "O", "frame": impl_frame, "at": cursor})
        all_events.append({"type": "O", "frame": cat_frame, "at": cursor})
        all_events.append({"type": "O", "frame": fn_frame, "at": cursor})
        cursor += gas
        all_events.append({"type": "C", "frame": fn_frame, "at": cursor})
        all_events.append({"type": "C", "frame": cat_frame, "at": cursor})
        all_events.append({"type": "C", "frame": impl_frame, "at": cursor})

    # --- Profile 4: On-chain gas from testnet transactions (hardcoded from ANALYSIS.md) ---
    # These are the REAL gas measurements that show 35-52x improvement.
    onchain_events_orig = []
    onchain_events_opt = []
    oc_orig = 0
    oc_opt = 0

    onchain_data = [
        # (label, original_mist, optimized_mist)
        ("compute_nd2 ×100", 74_400_000, 1_420_000),
        ("compute_nd2 ×200", 259_600_000, 7_450_000),
    ]
    for label, orig_gas, opt_gas in onchain_data:
        fr = get_frame(label)
        onchain_events_orig.append({"type": "O", "frame": fr, "at": oc_orig})
        oc_orig += orig_gas
        onchain_events_orig.append({"type": "C", "frame": fr, "at": oc_orig})

        onchain_events_opt.append({"type": "O", "frame": fr, "at": oc_opt})
        oc_opt += opt_gas
        onchain_events_opt.append({"type": "C", "frame": fr, "at": oc_opt})

    # Per-function on-chain breakdown (250 calls each)
    fn_breakdown = [
        ("sqrt ×250", 72_900_000, 1_000_000),
        ("CDF ×250", 6_180_000, 1_000_000),
        ("ln ×250", 2_020_000, 1_000_000),
    ]
    for label, orig_gas, opt_gas in fn_breakdown:
        fr = get_frame(label)
        onchain_events_orig.append({"type": "O", "frame": fr, "at": oc_orig})
        oc_orig += orig_gas
        onchain_events_orig.append({"type": "C", "frame": fr, "at": oc_orig})

        onchain_events_opt.append({"type": "O", "frame": fr, "at": oc_opt})
        oc_opt += opt_gas
        onchain_events_opt.append({"type": "C", "frame": fr, "at": oc_opt})

    profiles = [
        {
            "type": "evented",
            "name": "On-Chain Original (Testnet MIST)",
            "unit": "none",
            "startValue": 0,
            "endValue": oc_orig,
            "events": onchain_events_orig,
        },
        {
            "type": "evented",
            "name": "On-Chain Optimized (Testnet MIST)",
            "unit": "none",
            "startValue": 0,
            "endValue": oc_opt,
            "events": onchain_events_opt,
        },
        {
            "type": "evented",
            "name": "Original — Wall Time (test runner)",
            "unit": "nanoseconds",
            "startValue": 0,
            "endValue": orig_t,
            "events": orig_time_events,
        },
        {
            "type": "evented",
            "name": "Optimized — Wall Time (test runner)",
            "unit": "nanoseconds",
            "startValue": 0,
            "endValue": opt_t,
            "events": opt_time_events,
        },
        {
            "type": "evented",
            "name": "All Tests — Overview",
            "unit": "none",
            "startValue": 0,
            "endValue": cursor,
            "events": all_events,
        },
    ]

    return {
        "$schema": "https://www.speedscope.app/file-format-schema.json",
        "shared": {"frames": frames},
        "profiles": profiles,
        "name": "predict-math-opt gas profile",
        "activeProfileIndex": 0,
        "exporter": "predict-math-opt gas_to_speedscope.py",
    }


def print_summary(rows: list[dict]):
    """Print a quick comparison table to stdout."""
    benchmarks = {}
    for row in rows:
        name = row["name"]
        gas = int(row["gas"])
        nanos = int(row["nanos"])
        _, impl, label = classify_test(name)
        cat = classify_test(name)[0]
        if cat.startswith("unit_") or cat == "validation":
            continue
        benchmarks.setdefault(cat, {})[impl] = (gas, nanos, label)

    # Note: test-mode gas is capped at ~998k (Sui minimum bucket), so
    # wall-clock nanos are the meaningful metric from the test runner.
    # Real on-chain gas differences (35-52x) are measured via testnet txns.
    print("\n" + "=" * 78)
    print("BENCHMARK SUMMARY (wall-clock ns — test-mode gas is capped at ~998k floor)")
    print("=" * 78)

    print(f"\n{'Category':<20} {'Orig (ns)':>14} {'Opt (ns)':>14} {'Speedup':>10} {'Gas (test)':>12}")
    print("-" * 74)
    for cat in sorted(benchmarks.keys()):
        items = benchmarks[cat]
        if "original" in items and "optimized" in items:
            og, on, _ = items["original"]
            optg, optn, _ = items["optimized"]
            time_ratio = on / optn if optn > 0 else float("inf")
            print(f"{cat:<20} {on:>14,} {optn:>14,} {time_ratio:>9.1f}× {og:>12,}")
        else:
            for impl, (g, n, label) in items.items():
                print(f"{cat:<20} {label}: {n:>14,} ns")

    print("\nNote: On-chain gas (testnet) shows 35-52x improvement at 100-200 nd2 calls.")
    print("      Test-mode gas uses flat budget. Wall-clock ns reflects real compute cost.")
    print()


def main():
    os.chdir(PACKAGE_DIR)
    print("Running sui move test -s csv ...")
    rows = run_benchmarks()
    print(f"Parsed {len(rows)} test results.")

    print_summary(rows)

    profile = build_speedscope(rows)
    out_path = PACKAGE_DIR / "gas_profile.speedscope.json"
    with open(out_path, "w") as f:
        json.dump(profile, f, indent=2)

    print(f"Speedscope profile written to: {out_path}")
    print(f"Open https://www.speedscope.app/ and drag in the file, or:")
    print(f"  open https://www.speedscope.app/#profileURL=file://{out_path}")

    if "--open" in sys.argv:
        import webbrowser
        webbrowser.open(f"https://www.speedscope.app/")
        print("\nOpened speedscope.app in browser. Drag the JSON file into it.")


if __name__ == "__main__":
    main()

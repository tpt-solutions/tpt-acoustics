#!/usr/bin/env python3
"""Benchmark regression gate for criterion results.

Run after `cargo bench -- --baseline <name>`. Criterion then writes
`<criterion-dir>/<bench>/change/estimates.json` for every benchmark, with the
relative change of the mean (0.05 means 5 % slower) and its 95 % confidence
interval.

A benchmark fails the gate only if its mean change is above `--threshold`
*and* the lower bound of its confidence interval is above zero. The second
condition keeps shared-runner noise from failing a pull request.

The gate also fails when no comparison is found. A gate that passes because
nothing was measured is not a gate.

Usage:
    python3 scripts/bench_gate.py [--criterion-dir target/criterion] [--threshold 0.15]

Run `python3 scripts/bench_gate.py --self-test` to check the logic.
"""

import argparse
import json
import sys
import unittest
from pathlib import Path


def load_changes(criterion_dir):
    """Map benchmark name to (mean change, lower CI bound) for every comparison."""
    changes = {}
    for path in sorted(Path(criterion_dir).rglob("change/estimates.json")):
        name = str(path.parent.parent.relative_to(criterion_dir)).replace("\\", "/")
        data = json.loads(path.read_text())
        mean = data["mean"]
        changes[name] = (mean["point_estimate"], mean["confidence_interval"]["lower_bound"])
    return changes


def find_regressions(changes, threshold):
    """Names whose slowdown exceeds `threshold` with a significant lower bound."""
    return sorted(
        name
        for name, (point, lower) in changes.items()
        if point > threshold and lower > 0.0
    )


def run(criterion_dir, threshold):
    changes = load_changes(criterion_dir)
    if not changes:
        print(f"bench gate: no comparisons found under {criterion_dir}; run with --baseline first")
        return 1
    print(f"{'benchmark':<40} {'mean change':>12} {'95% lower':>11}")
    for name, (point, lower) in changes.items():
        print(f"{name:<40} {point * 100:>+11.2f}% {lower * 100:>+10.2f}%")
    regressions = find_regressions(changes, threshold)
    if regressions:
        print(
            f"bench gate: FAILED, {len(regressions)} regression(s) above "
            f"{threshold * 100:.0f}% with significance: {', '.join(regressions)}"
        )
        return 1
    print(f"bench gate: passed ({len(changes)} benchmarks, threshold {threshold * 100:.0f}%)")
    return 0


class GateLogicTests(unittest.TestCase):
    def test_flags_significant_slowdown_only(self):
        changes = {
            "slow_significant": (0.30, 0.10),
            "slow_noisy": (0.30, -0.02),
            "small": (0.05, 0.01),
            "faster": (-0.20, -0.25),
        }
        self.assertEqual(find_regressions(changes, 0.15), ["slow_significant"])

    def test_threshold_is_exclusive(self):
        self.assertEqual(find_regressions({"edge": (0.15, 0.05)}, 0.15), [])

    def test_loads_nested_benchmark_names(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, point, lower in [("group/a", 0.2, 0.1), ("b", -0.1, -0.2)]:
                folder = root / name / "change"
                folder.mkdir(parents=True)
                payload = {"mean": {"point_estimate": point,
                                    "confidence_interval": {"lower_bound": lower}}}
                (folder / "estimates.json").write_text(json.dumps(payload))
            loaded = load_changes(root)
            self.assertEqual(sorted(loaded), ["b", "group/a"])
            self.assertEqual(loaded["group/a"], (0.2, 0.1))

    def test_empty_directory_fails(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(run(Path(tmp), 0.15), 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--criterion-dir", default="target/criterion")
    parser.add_argument("--threshold", type=float, default=0.15,
                        help="allowed mean slowdown as a fraction (default 0.15)")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(GateLogicTests)
        return 0 if unittest.TextTestRunner(verbosity=1).run(suite).wasSuccessful() else 1
    return run(Path(args.criterion_dir), args.threshold)


if __name__ == "__main__":
    sys.exit(main())

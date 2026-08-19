#!/usr/bin/env python3
"""Visualize one or more comma-separated perf-stat reports."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path, help="perf CSV reports")
    parser.add_argument("--output", type=Path, required=True, help="output image")
    return parser.parse_args()


def read_perf(path: Path) -> dict[str, float]:
    counters: dict[str, float] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.reader(line for line in stream if not line.startswith("#")):
            if len(row) < 3 or not row[0].strip() or "not" in row[0]:
                continue
            try:
                value = float(row[0].strip())
            except ValueError:
                continue
            event = row[2].strip().split(":", 1)[0]
            counters[event] = value
    return counters


def kernel_repetitions(path: Path) -> int:
    result_path = path.with_suffix(".result.txt")
    if not result_path.is_file():
        raise SystemExit(f"missing numerical result beside perf report: {result_path}")
    repetitions = 0
    with result_path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("BENCH_RESULT,"):
                fields = line.rstrip().split(",")
                repetitions = int(fields[7])
    if repetitions < 1:
        raise SystemExit(f"no valid BENCH_RESULT in {result_path}")
    return repetitions


def main() -> None:
    args = arguments()
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "qrlinalg-matplotlib"))
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as error:
        raise SystemExit("matplotlib and numpy are required for perf plots") from error

    reports = []
    for path in args.inputs:
        repetitions = kernel_repetitions(path)
        reports.append((path, {event: value / repetitions for event, value in read_perf(path).items()}))
    if not any(counters for _, counters in reports):
        raise SystemExit("no numeric perf counters found")

    labels = [path.parent.name.replace("_", " ") for path, _ in reports]
    count_events = ["cycles", "instructions", "branches", "branch-misses", "cache-references", "cache-misses"]
    x_values = np.arange(len(labels))
    width = 0.8 / len(count_events)
    figure = plt.figure(figsize=(max(11, 1.5 * len(labels)), 9), constrained_layout=True)
    layout = figure.add_gridspec(2, 3, height_ratios=(1.5, 1.0), hspace=0.36, wspace=0.25)
    count_axis = figure.add_subplot(layout[0, :])
    rate_axes = [figure.add_subplot(layout[1, index]) for index in range(3)]

    for event_index, event in enumerate(count_events):
        values = [counters.get(event, float("nan")) for _, counters in reports]
        offset = (event_index - (len(count_events) - 1) / 2.0) * width
        count_axis.bar(x_values + offset, values, width=width, label=event)
    count_axis.set_yscale("log")
    count_axis.set_ylabel("Events per numerical operation")
    count_axis.set_title("Internally gated perf hardware counters")
    count_axis.set_xticks(x_values, labels)
    count_axis.grid(True, axis="y", which="both", alpha=0.25)
    count_axis.legend(ncol=3)

    ipc = []
    branch_miss = []
    cache_miss = []
    for _, counters in reports:
        ipc.append(counters.get("instructions", 0.0) / max(counters.get("cycles", 0.0), 1.0))
        branch_miss.append(100.0 * counters.get("branch-misses", 0.0) / max(counters.get("branches", 0.0), 1.0))
        cache_miss.append(100.0 * counters.get("cache-misses", 0.0) / max(counters.get("cache-references", 0.0), 1.0))
    rate_series = [
        (ipc, "Instructions per cycle", "IPC", "#4c78a8"),
        (branch_miss, "Branch misses", "Percent", "#f58518"),
        (cache_miss, "Cache misses", "Percent", "#54a24b"),
    ]
    for axis, (values, title, ylabel, color) in zip(rate_axes, rate_series):
        axis.bar(x_values, values, width=0.68, color=color)
        axis.set_title(title)
        axis.set_ylabel(ylabel)
        axis.set_xticks(x_values, labels, rotation=40, ha="right")
        axis.grid(True, axis="y", alpha=0.25)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    print(args.output)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Plot verbose benchmark CSV records without rerunning numerical code."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from collections import defaultdict
from pathlib import Path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path, help="QR or LDLT result CSV files")
    parser.add_argument("--output", type=Path, required=True, help="output directory")
    return parser.parse_args()


def main() -> None:
    args = arguments()
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "qrlinalg-matplotlib"))
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise SystemExit("matplotlib is required: python3 -m pip install matplotlib") from error

    records: list[dict[str, str]] = []
    for input_path in args.inputs:
        with input_path.open(newline="", encoding="utf-8") as stream:
            records.extend(csv.DictReader(stream))
    if not records:
        raise SystemExit("no benchmark records found")
    args.output.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for record in records:
        grouped[(record["implementation"], record["kind"])].append(record)

    for (implementation, kind), group in sorted(grouped.items()):
        figure, axis = plt.subplots(figsize=(9, 5.5))
        operations = sorted({record["operation"] for record in group})
        for operation in operations:
            selected = sorted(
                (record for record in group if record["operation"] == operation),
                key=lambda record: int(record["n"]),
            )
            sizes = [int(record["n"]) for record in selected]
            milliseconds = [1000.0 * float(record["seconds_per_operation"]) for record in selected]
            axis.plot(sizes, milliseconds, marker="o", label=operation.replace("_", " "))
        axis.set_title(f"{implementation.upper()} {kind} operation timings")
        axis.set_xlabel("Matrix order")
        axis.set_ylabel("Milliseconds per operation")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend()
        figure.tight_layout()
        output_file = args.output / f"{implementation}_{kind}_timings.png"
        figure.savefig(output_file, dpi=180)
        plt.close(figure)
        print(output_file)

    by_key = {
        (record["implementation"], record["kind"], int(record["n"]), record["operation"]): record
        for record in records
    }
    comparison_pairs = {
        "fresh factorization": ("factorization", "factorization"),
        "factorization + solve": ("full_solve", "full_solve"),
        "solve with factors": ("solve", "solve"),
        "append / LDLT append": ("append", "append_factorization"),
        "replacement / LDLT rebuild": ("replacement", "factorization"),
        "deletion / LDLT rebuild": ("deletion", "factorization"),
    }
    for kind in sorted({record["kind"] for record in records}):
        figure, axis = plt.subplots(figsize=(9, 5.5))
        plotted = False
        sizes = sorted({int(record["n"]) for record in records if record["kind"] == kind})
        for label, (qr_operation, ldlt_operation) in comparison_pairs.items():
            points = []
            for matrix_n in sizes:
                qr_record = by_key.get(("qr", kind, matrix_n, qr_operation))
                ldlt_record = by_key.get(("ldlt", kind, matrix_n, ldlt_operation))
                if qr_record and ldlt_record:
                    qr_time = float(qr_record["seconds_per_operation"])
                    ldlt_time = float(ldlt_record["seconds_per_operation"])
                    points.append((matrix_n, ldlt_time / qr_time))
            if points:
                plotted = True
                axis.plot([point[0] for point in points], [point[1] for point in points], marker="o", label=label)
        if plotted:
            axis.axhline(1.0, color="black", linewidth=1.0, linestyle="--")
            axis.text(sizes[0], 1.08, "QR faster above this line", fontsize=9)
            axis.set_title(f"QR speed relative to LDLT — {kind}")
            axis.set_xlabel("Matrix order")
            axis.set_ylabel("LDLT time / QR time")
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.legend(fontsize=8, ncol=2)
            figure.tight_layout()
            output_file = args.output / f"qr_vs_ldlt_{kind}_speedup.png"
            figure.savefig(output_file, dpi=200)
            print(output_file)
        plt.close(figure)

    print("\nRecorded convergence metrics:")
    for record in sorted(records, key=lambda item: (item["implementation"], item["kind"], int(item["n"]), item["operation"])):
        milliseconds = 1000.0 * float(record["seconds_per_operation"])
        print(
            f"{record['implementation']:4s} {record['kind']:7s} n={int(record['n']):4d} "
            f"{record['operation']:22s} {milliseconds:11.4f} ms  "
            f"iter={int(record['iterations']):2d} lambda={float(record['eigenvalue']): .12e} "
            f"residual={float(record['residual']):.3e} status={record['status'].strip()}"
        )


if __name__ == "__main__":
    main()

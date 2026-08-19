#!/usr/bin/env python3
"""Generate an auditable Markdown report from independent QR and LDLT CSVs."""

from __future__ import annotations

import argparse
import csv
import platform
from datetime import datetime, timezone
from pathlib import Path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qr", type=Path, required=True)
    parser.add_argument("--ldlt", type=Path, required=True)
    parser.add_argument("--plots", type=Path, required=True)
    parser.add_argument("--perf-root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
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
            counters[row[2].strip().split(":", 1)[0]] = value
    return counters


def perf_repetitions(path: Path) -> int:
    result_path = path.with_suffix(".result.txt")
    repetitions = 0
    if result_path.is_file():
        for line in result_path.read_text(encoding="utf-8").splitlines():
            if line.startswith("BENCH_RESULT,"):
                repetitions = int(line.split(",")[7])
    return repetitions


def perf_table(perf_root: Path, implementation: str, kind: str) -> list[str]:
    lines = [
        "| Operation | Instructions/op | Cycles/op | IPC | Branch miss | Cache miss |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    reports: list[tuple[str, dict[str, float]]] = []
    for path in sorted(perf_root.glob(f"{implementation}/*/{kind}_1000.csv")):
        repetitions = perf_repetitions(path)
        if repetitions < 1:
            continue
        reports.append((path.parent.name, {key: value / repetitions for key, value in read_perf(path).items()}))
    for operation, counters in reports:
        instructions = counters.get("instructions", 0.0)
        cycles = counters.get("cycles", 0.0)
        branches = counters.get("branches", 0.0)
        branch_misses = counters.get("branch-misses", 0.0)
        cache_references = counters.get("cache-references", 0.0)
        cache_misses = counters.get("cache-misses", 0.0)
        ipc = instructions / max(cycles, 1.0)
        branch_rate = 100.0 * branch_misses / max(branches, 1.0)
        cache_rate = 100.0 * cache_misses / max(cache_references, 1.0)
        lines.append(
            f"| {operation.replace('_', ' ')} | {instructions:.3e} | {cycles:.3e} | "
            f"{ipc:.3f} | {branch_rate:.3f}% | {cache_rate:.3f}% |"
        )
    return lines


def read_records(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def milliseconds(record: dict[str, str]) -> float:
    return 1000.0 * float(record["seconds_per_operation"])


def timing_table(records: list[dict[str, str]], implementation: str, kind: str) -> list[str]:
    selected = [record for record in records if record["implementation"] == implementation and record["kind"] == kind]
    operations = sorted({record["operation"] for record in selected})
    sizes = sorted({int(record["n"]) for record in selected})
    indexed = {(int(record["n"]), record["operation"]): record for record in selected}
    lines = [
        "| n | " + " | ".join(operation.replace("_", " ") for operation in operations) + " |",
        "|---:|" + "|".join("---:" for _ in operations) + "|",
    ]
    for matrix_n in sizes:
        values = []
        for operation in operations:
            record = indexed.get((matrix_n, operation))
            values.append(f"{milliseconds(record):.4f}" if record else "—")
        lines.append(f"| {matrix_n} | " + " | ".join(values) + " |")
    return lines


def comparison_table(records: list[dict[str, str]], kind: str) -> list[str]:
    indexed = {
        (record["implementation"], int(record["n"]), record["operation"]): record
        for record in records
        if record["kind"] == kind
    }
    pairs = [
        ("Fresh factorization", "factorization", "factorization"),
        ("Solve with existing factors", "solve", "solve"),
        ("Factorization plus solve", "full_solve", "full_solve"),
        ("Append", "append", "append_factorization"),
        ("Replacement / LDLT rebuild", "replacement", "factorization"),
        ("Deletion / LDLT rebuild", "deletion", "factorization"),
    ]
    sizes = sorted({key[1] for key in indexed})
    lines = [
        "| n | Comparison | QR (ms) | LDLT baseline (ms) | LDLT/QR | Faster |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for matrix_n in sizes:
        for label, qr_operation, ldlt_operation in pairs:
            qr_record = indexed.get(("qr", matrix_n, qr_operation))
            ldlt_record = indexed.get(("ldlt", matrix_n, ldlt_operation))
            if not qr_record or not ldlt_record:
                continue
            qr_time = milliseconds(qr_record)
            ldlt_time = milliseconds(ldlt_record)
            ratio = ldlt_time / qr_time
            faster = "QR" if ratio > 1.0 else "LDLT"
            lines.append(f"| {matrix_n} | {label} | {qr_time:.4f} | {ldlt_time:.4f} | {ratio:.3f}× | {faster} |")
    return lines


def main() -> None:
    args = arguments()
    records = read_records(args.qr) + read_records(args.ldlt)
    if not records:
        raise SystemExit("no benchmark records found")
    statuses = {record["status"].strip() for record in records}
    maximum_residual = max(float(record["residual"]) for record in records)
    iterations = sorted({int(record["iterations"]) for record in records})
    precision = sorted({record["precision"] for record in records})
    plot_path = Path(args.plots)

    lines = [
        "# QR and LDLT benchmark report",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        "",
        f"Platform: `{platform.platform()}`",
        "",
        f"Working precision kind(s): `{', '.join(precision)}`",
        "",
        f"Statuses observed: `{', '.join(sorted(statuses))}`; inverse-iteration counts: `{iterations}`; maximum normalized residual: `{maximum_residual:.3e}`.",
        "",
        "All times below are milliseconds per operation. Input, allocation, warm-up, destructive-input copies, and post-operation validation are outside the timed region.",
        "",
    ]
    for implementation in ("qr", "ldlt"):
        lines.extend([f"## {implementation.upper()} timings", ""])
        for kind in ("real", "complex"):
            lines.extend([f"### {kind.capitalize()}", ""])
            lines.extend(timing_table(records, implementation, kind))
            lines.extend(["", f"![{implementation} {kind} timings]({plot_path.name}/{implementation}_{kind}_timings.png)", ""])

    lines.extend(["## Explicit QR-versus-LDLT comparisons", ""])
    for kind in ("real", "complex"):
        lines.extend([f"### {kind.capitalize()}", ""])
        lines.extend(comparison_table(records, kind))
        lines.extend(["", f"![QR versus LDLT {kind}]({plot_path.name}/qr_vs_ldlt_{kind}_speedup.png)", ""])

    if args.perf_root:
        lines.extend([
            "## Controlled hardware-counter measurements",
            "",
            "These measurements use order 1000. Counters are enabled and disabled at the same internal boundaries as the wall-clock timer, then normalized by the number of kernel operations.",
            "",
        ])
        for implementation in ("qr", "ldlt"):
            for kind in ("real", "complex"):
                lines.extend([f"### {implementation.upper()} {kind}", ""])
                lines.extend(perf_table(args.perf_root, implementation, kind))
                lines.extend(["", f"![{implementation} {kind} perf]({plot_path.name}/perf_{implementation}_{kind}.png)", ""])

    lines.extend([
        "## Limitations",
        "",
        "- These are single-machine, single-process results. CPU affinity, background load, thermal state, and frequency policy were not controlled by the harness.",
        "- Five wall-clock repetitions provide a useful exploratory mean but not a confidence interval or robust distribution. Repeat longer runs before making production decisions.",
        "- Fresh QR includes forming `H-shift*S` and explicitly generating dense `Q`. LDLT factorization receives an already formed shifted matrix and does not construct an equivalent orthogonal factor. Fresh-factorization timings therefore describe the public operations, not equal low-level work.",
        "- QR append and deletion copy a prepared state outside timing and measure the final index. Other indices and long uninterrupted update sequences can have different cache and numerical behavior.",
        "- Replacement uses deterministic perturbations derived from one matrix column (`1e-6*H(:,i)` and `1e-7*S(:,i)`) and alternates their signs. It measures the update kernel, not a distribution of physical changes.",
        "- The pristine LDLT implementation has no replacement or deletion primitive. Those comparison rows use full refactorization as an explicitly labelled rebuild baseline.",
        "- All dataset files store binary64 values. Builds with `wp=10` or `wp=16` exercise wider arithmetic but do not gain wider-precision input data.",
        "- Shifts, all-ones starting vectors, tolerance `1e-12`, iteration limit 30, and Euclidean normalization reproduce the original Claude driver. Different shifts or starting vectors change iteration counts and solve costs.",
        "- The bundled generic BLAS/LAPACK is measured as built by this repository; results do not predict performance with a tuned vendor BLAS or a redesigned implicit-Q backend.",
        "- Perf counters are hardware- and kernel-dependent. The perf script gates counters at internal timer boundaries, but event availability, multiplexing, and kernel permissions still affect measurements.",
        "- Numerical validation is intentionally outside timing. Successful status and a small residual establish a usable eigenpair but do not prove bitwise equivalence between QR and LDLT factors.",
        "",
        "## Raw data",
        "",
        f"- QR: `{args.qr}`",
        f"- LDLT: `{args.ldlt}`",
    ])
    if args.perf_root:
        lines.append(f"- Perf: `{args.perf_root}`")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()

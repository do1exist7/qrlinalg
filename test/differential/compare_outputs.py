#!/usr/bin/env python3
"""Run and compare qrlinalg and original inverse-iteration implementations."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path
import re
import subprocess
import sys


REQUIRED_COLUMNS = {
    "case_id",
    "kind",
    "h_file",
    "s_file",
    "shift",
    "tol",
    "max_iter",
    "norm_mode",
}


@dataclass(frozen=True)
class Result:
    kind: str
    order: int
    info: int
    iterations: int
    eigenvalue: Decimal
    relative_accuracy: Decimal
    residual: Decimal
    epsilon: Decimal
    vector: tuple[tuple[Decimal, Decimal], ...]


def decimal_value(text: str, field: str) -> Decimal:
    try:
        value = Decimal(text.replace("D", "E").replace("d", "e"))
    except InvalidOperation as error:
        raise ValueError(f"invalid {field}: {text!r}") from error
    if not value.is_finite():
        raise ValueError(f"non-finite {field}: {text!r}")
    return value


def parse_result(output: str, implementation: str) -> Result:
    header: list[str] | None = None
    components: dict[int, tuple[Decimal, Decimal]] = {}

    for line in output.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "DIFF_RESULT":
            if header is not None:
                raise ValueError(f"{implementation} emitted multiple results")
            header = fields
        elif fields[0] == "DIFF_VECTOR":
            if len(fields) != 4:
                raise ValueError(f"malformed {implementation} vector record: {line}")
            index = int(fields[1])
            if index in components:
                raise ValueError(
                    f"{implementation} emitted duplicate vector index {index}"
                )
            components[index] = (
                decimal_value(fields[2], "vector real component"),
                decimal_value(fields[3], "vector imaginary component"),
            )

    if header is None:
        raise ValueError(f"{implementation} emitted no DIFF_RESULT record")
    if len(header) != 9:
        raise ValueError(f"malformed {implementation} result record")

    kind = header[1]
    order = int(header[2])
    vector = tuple(components.get(index, ()) for index in range(1, order + 1))
    if any(len(component) != 2 for component in vector):
        raise ValueError(f"{implementation} emitted an incomplete eigenvector")
    if len(components) != order:
        raise ValueError(f"{implementation} emitted duplicate/out-of-range indices")

    return Result(
        kind=kind,
        order=order,
        info=int(header[3]),
        iterations=int(header[4]),
        eigenvalue=decimal_value(header[5], "eigenvalue"),
        relative_accuracy=decimal_value(header[6], "relative accuracy"),
        residual=decimal_value(header[7], "residual"),
        epsilon=decimal_value(header[8], "epsilon"),
        vector=vector,
    )


def run_driver(command: list[str], label: str, output_base: Path) -> Result:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    stdout_path = Path(f"{output_base}.{label}.out")
    stderr_path = Path(f"{output_base}.{label}.err")
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"{label} driver exited with status {completed.returncode}; "
            f"see {stderr_path}"
        )
    return parse_result(completed.stdout, label)


def status_name(implementation: str, info: int) -> str:
    mappings = {
        "new": {0: "success", 5: "singular", 6: "no_convergence"},
        "orig": {0: "success", 1: "singular", 2: "no_convergence"},
    }
    return mappings[implementation].get(info, f"error_{info}")


def normalized(vector: tuple[tuple[Decimal, Decimal], ...]):
    norm_squared = sum((real * real + imag * imag for real, imag in vector), Decimal(0))
    if norm_squared <= 0:
        raise ValueError("cannot compare a zero eigenvector")
    norm = norm_squared.sqrt()
    return tuple((real / norm, imag / norm) for real, imag in vector)


def aligned_vector_error(reference: Result, candidate: Result) -> Decimal:
    reference_vector = normalized(reference.vector)
    candidate_vector = normalized(candidate.vector)

    if reference.kind == "real":
        overlap = sum(
            (reference_part[0] * candidate_part[0]
             for reference_part, candidate_part
             in zip(reference_vector, candidate_vector)),
            Decimal(0),
        )
        sign = Decimal(-1) if overlap < 0 else Decimal(1)
        error_squared = sum(
            ((candidate_part[0] * sign - reference_part[0]) ** 2
             for reference_part, candidate_part
             in zip(reference_vector, candidate_vector)),
            Decimal(0),
        )
        return error_squared.sqrt()

    overlap_real = Decimal(0)
    overlap_imag = Decimal(0)
    for (reference_real, reference_imag), (candidate_real, candidate_imag) in zip(
        reference_vector, candidate_vector
    ):
        overlap_real += reference_real * candidate_real + reference_imag * candidate_imag
        overlap_imag += reference_real * candidate_imag - reference_imag * candidate_real

    overlap_magnitude = (overlap_real**2 + overlap_imag**2).sqrt()
    if overlap_magnitude == 0:
        return Decimal(2).sqrt()
    phase_real = overlap_real / overlap_magnitude
    phase_imag = overlap_imag / overlap_magnitude

    error_squared = Decimal(0)
    for (reference_real, reference_imag), (candidate_real, candidate_imag) in zip(
        reference_vector, candidate_vector
    ):
        aligned_real = candidate_real * phase_real + candidate_imag * phase_imag
        aligned_imag = candidate_imag * phase_real - candidate_real * phase_imag
        error_squared += (aligned_real - reference_real) ** 2
        error_squared += (aligned_imag - reference_imag) ** 2
    return error_squared.sqrt()


def optional_tolerance(row: dict[str, str], name: str, default: Decimal) -> Decimal:
    text = (row.get(name) or "").strip()
    if not text:
        return default
    value = decimal_value(text, name)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def scientific(value: Decimal) -> str:
    if value == 0:
        return "0.000E+00"
    return f"{value:.3E}"


def compare_case(row: dict[str, str], new_result: Result, orig_result: Result) -> str:
    expected_kind = row["kind"].strip().lower()
    if new_result.kind != expected_kind or orig_result.kind != expected_kind:
        raise ValueError("driver result kind does not match the manifest")
    if new_result.order != orig_result.order:
        raise ValueError(
            f"matrix-order mismatch: new={new_result.order}, orig={orig_result.order}"
        )

    new_status = status_name("new", new_result.info)
    orig_status = status_name("orig", orig_result.info)
    if new_status != orig_status:
        raise ValueError(
            f"solver-status mismatch: new={new_status}, orig={orig_status}"
        )
    if new_status == "singular":
        return "both implementations detected a singular shifted matrix"
    if new_status not in {"success", "no_convergence"}:
        raise ValueError(f"both drivers returned unsupported status {new_status}")

    requested_tolerance = abs(decimal_value(row["tol"], "tol"))
    arithmetic_floor = (
        Decimal(10000) * Decimal(new_result.order)
        * max(new_result.epsilon, orig_result.epsilon)
    )
    default_tolerance = max(Decimal(100) * requested_tolerance, arithmetic_floor)
    eigen_tolerance = optional_tolerance(row, "eigen_rtol", default_tolerance)
    vector_tolerance = optional_tolerance(row, "vector_rtol", default_tolerance)
    residual_tolerance = optional_tolerance(row, "residual_tol", default_tolerance)

    eigen_scale = max(
        Decimal(1), abs(new_result.eigenvalue), abs(orig_result.eigenvalue)
    )
    eigen_error = abs(new_result.eigenvalue - orig_result.eigenvalue) / eigen_scale
    vector_error = aligned_vector_error(orig_result, new_result)

    failures: list[str] = []
    if eigen_error > eigen_tolerance:
        failures.append(f"eigenvalue relative error {eigen_error} > {eigen_tolerance}")
    if vector_error > vector_tolerance:
        failures.append(f"aligned vector error {vector_error} > {vector_tolerance}")
    if new_result.residual > residual_tolerance:
        failures.append(
            f"new residual {new_result.residual} > {residual_tolerance}"
        )
    if orig_result.residual > residual_tolerance:
        failures.append(
            f"original residual {orig_result.residual} > {residual_tolerance}"
        )
    if failures:
        raise ValueError("; ".join(failures))

    return (
        f"status={new_status}, eigen_rel={scientific(eigen_error)}, "
        f"vector_rel={scientific(vector_error)}, "
        f"residuals={scientific(new_result.residual)}/"
        f"{scientific(orig_result.residual)}, "
        f"rel_acc={scientific(new_result.relative_accuracy)}/"
        f"{scientific(orig_result.relative_accuracy)}, "
        f"iterations={new_result.iterations}/{orig_result.iterations}"
    )


def load_cases(manifest: Path) -> list[dict[str, str]]:
    with manifest.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(
            line for line in stream if line.strip() and not line.lstrip().startswith("#")
        )
        if reader.fieldnames is None:
            raise ValueError("manifest has no header")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"manifest lacks columns: {', '.join(sorted(missing))}")
        cases = []
        for row_number, row in enumerate(reader, start=2):
            missing_values = [
                column
                for column in REQUIRED_COLUMNS
                if not (row.get(column) or "").strip()
            ]
            if missing_values:
                raise ValueError(
                    f"manifest row {row_number} has blank fields: "
                    f"{', '.join(sorted(missing_values))}"
                )
            cases.append(dict(row))
    if not cases:
        raise ValueError("manifest contains no cases")
    return cases


def case_arguments(row: dict[str, str], manifest: Path) -> list[str]:
    kind = row["kind"].strip().lower()
    if kind not in {"real", "complex"}:
        raise ValueError(f"kind must be real or complex, got {kind!r}")

    paths: list[str] = []
    for field in ("h_file", "s_file"):
        path = Path(row[field].strip())
        if not path.is_absolute():
            path = manifest.parent / path
        path = path.resolve()
        if not path.is_file():
            raise ValueError(f"{field} does not exist: {path}")
        paths.append(str(path))

    decimal_value(row["shift"], "shift")
    decimal_value(row["tol"], "tol")
    max_iter = int(row["max_iter"])
    if max_iter <= 0:
        raise ValueError("max_iter must be positive")
    norm_mode = int(row["norm_mode"])
    return [kind, *paths, row["shift"], row["tol"], str(max_iter), str(norm_mode)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--new-exe", type=Path, required=True)
    parser.add_argument("--orig-exe", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()

    manifest = arguments.manifest.resolve()
    if not manifest.is_file():
        parser.error(f"manifest does not exist: {manifest}")

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    failures = 0
    try:
        cases = load_cases(manifest)
    except (OSError, ValueError) as error:
        parser.error(str(error))

    with localcontext() as context:
        context.prec = 80
        for row in cases:
            case_id = row["case_id"].strip()
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", case_id):
                print(f"FAIL: invalid case_id {case_id!r}", file=sys.stderr)
                failures += 1
                continue
            try:
                driver_arguments = case_arguments(row, manifest)
                output_base = arguments.output_dir / case_id
                new_result = run_driver(
                    [str(arguments.new_exe.resolve()), *driver_arguments],
                    "new",
                    output_base,
                )
                orig_result = run_driver(
                    [str(arguments.orig_exe.resolve()), *driver_arguments],
                    "orig",
                    output_base,
                )
                summary = compare_case(row, new_result, orig_result)
                print(f"PASS: {case_id}: {summary}")
            except (OSError, RuntimeError, ValueError) as error:
                print(f"FAIL: {case_id}: {error}", file=sys.stderr)
                failures += 1

    if failures:
        print(f"FAIL: differential comparison failed for {failures} case(s)")
        return 1
    print(f"PASS: differential comparison completed for {len(cases)} case(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Bounded, sequential PLT runtime probes, without a security-bit target.

Every point runs in a fresh Linux process with wall/CPU/address-space limits.
Timeout/error/skip results remain in the report. Builds are excluded from timing.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
Q31 = 2**31 - 1
Q61 = 2**61 - 1
POINTS = [
    (1, 257, 2, 3), (1, 257, 3, 5),
    (1, 65537, 2, 3), (1, 16777213, 2, 3),
    (1, Q31, 2, 3), (1, Q61, 2, 3),
    (2, Q31, 2, 3), (2, Q61, 2, 3),
    (4, Q31, 2, 3), (4, Q61, 2, 3), (8, Q31, 2, 3),
    (1, Q31, 3, 5), (2, Q31, 3, 5), (4, Q31, 3, 5),
    (1, Q31, 5, 10), (2, Q31, 5, 10),
    (8, Q61, 3, 5), (16, Q31, 2, 3), (64, Q61, 2, 3),
]


def plan(n, q, t, k):
    ell = (q - 1).bit_length()
    m = n * (ell + 1)
    dimension = 2 * m * t
    return {"n": n, "q": q, "t": t, "k": k, "m": m,
            "target_rows": n * t, "preimage_rows": dimension,
            "single_dense_f64_square_mib": 8 * dimension**2 / 2**20,
            "all_u_i64_payload_lower_mib": 8 * k * dimension * m / 2**20,
            "memory_estimate_note": "payload only; ignores big integers, rationals, copies and allocator overhead"}


def version(command):
    return subprocess.check_output(command, cwd=ROOT, text=True).strip()


def process_peak(pid):
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError):
        pass
    return 0


def save(path, report):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    temporary.replace(path)


def run_point(args, point):
    n, q, t, k = point
    case = f"n{n}_q{q}_{t}of{k}"
    result = {"case": case, "plan": plan(*point)}
    if result["plan"]["preimage_rows"] > args.max_dimension:
        return result | {"status": "skipped_dimension_budget", "probe": None}
    output = args.output_dir / f"{case}.json"
    # A partial checkpoint belongs only to this invocation.
    if output.exists():
        output.unlink()
    command = [str(ROOT / "target/release/examples/plt_parameter_probe"),
               "--n", str(n), "--modulus", str(q), "--threshold", str(t),
               "--participants", str(k), "--trials", str(args.trials),
               "--public-attempts", str(args.public_attempts), "--output", str(output)]

    def limits():
        memory = args.memory_mib * 2**20
        resource.setrlimit(resource.RLIMIT_AS, (memory, memory))
        cpu = math.ceil(args.timeout_seconds) + 2
        resource.setrlimit(resource.RLIMIT_CPU, (cpu, cpu))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))

    start = time.monotonic()
    process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True, preexec_fn=limits)
    timed_out, peak = False, 0
    while process.poll() is None:
        peak = max(peak, process_peak(process.pid))
        if time.monotonic() - start > args.timeout_seconds:
            process.kill()
            timed_out = True
            break
        time.sleep(0.05)
    stdout, stderr = process.communicate()
    probe = json.loads(output.read_text()) if output.exists() else None
    status = ("timeout" if timed_out else
              probe["status"] if probe and process.returncode == 0 else
              "algorithm_error" if probe and probe["status"] == "algorithm_error" else
              "process_error")
    return result | {"status": status, "wall_seconds": time.monotonic() - start,
                     "observed_peak_rss_kib": peak, "returncode": process.returncode,
                     "stdout": stdout, "stderr": stderr, "probe": probe}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout-seconds", type=float, default=30)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--max-dimension", type=int, default=1024)
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--public-attempts", type=int, default=64)
    parser.add_argument("--case", action="append", help="case label; repeat to select points")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "target/plt_parameter_probe")
    parser.add_argument("--report", type=Path, default=ROOT / "target/plt_parameter_exploration.json")
    parser.add_argument("--no-build", action="store_true")
    args = parser.parse_args()
    if min(args.timeout_seconds, args.memory_mib, args.max_dimension,
           args.trials, args.public_attempts) <= 0:
        parser.error("all budgets and sample counts must be positive")
    if platform.system() != "Linux":
        parser.error("this resource monitor currently supports Linux only")
    points = [p for p in POINTS if not args.case or f"n{p[0]}_q{p[1]}_{p[2]}of{p[3]}" in args.case]
    if args.case and len(points) != len(set(args.case)):
        parser.error("unknown or repeated case selection")
    args.output_dir = args.output_dir.resolve()
    args.report = args.report.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    if not args.no_build:
        subprocess.run(["cargo", "build", "--locked", "-p", "plt", "--release", "--example", "plt_parameter_probe"],
                       cwd=ROOT, check=True)
    sources = sorted(set(ROOT.glob("crates/*/src/**/*.rs"))
                     | set(ROOT.glob("crates/*/Cargo.toml")) | {
        ROOT / "Cargo.lock", ROOT / "Cargo.toml", Path(__file__).resolve(),
        ROOT / "crates/plt/examples/plt_parameter_probe.rs",
        ROOT / "crates/plt/examples/support/public_solution.rs"})
    cpu_lines = Path("/proc/cpuinfo").read_text().splitlines()
    report = {"scope": "runtime range, not security certification", "security_bit_target": None,
              "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "host": {"platform": platform.platform(), "logical_cpus": os.cpu_count(),
                       "cpu": next((line.split(":", 1)[1].strip() for line in cpu_lines
                                    if line.startswith("model name")), "unknown"),
                       "rustc": version(["rustc", "--version"]), "cargo": version(["cargo", "--version"])},
              "limits": {"timeout_seconds_per_point": args.timeout_seconds,
                         "address_space_mib": args.memory_mib, "max_dimension": args.max_dimension,
                         "trials_per_point": args.trials, "public_attempts_per_point": args.public_attempts},
              "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources}, "results": []}
    save(args.report, report)
    for point in points:
        result = run_point(args, point)
        report["results"].append(result)
        save(args.report, report)
        print(f"{result['case']}: {result['status']} "
              f"{result.get('wall_seconds', 0):.2f}s; "
              f"phase={(result.get('probe') or {}).get('phase', 'not started')}", flush=True)
    report["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    save(args.report, report)
    print(f"Saved {args.report}")


if __name__ == "__main__":
    main()

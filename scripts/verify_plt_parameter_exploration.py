#!/usr/bin/env python3
"""Independently check saved runtime probes and optionally exhaust two-column
public solutions when A has exactly two rows. No secret inputs are used.
"""

import argparse
import hashlib
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path
import time

from verify_plt_report import matrix, multiply, congruent, norm_squared, verify_profile


def verify(report):
    assert report["security_bit_target"] is None
    statuses = {"completed", "algorithm_error", "process_error", "timeout", "skipped_dimension_budget"}
    seen = set()
    for item in report["results"]:
        assert item["case"] not in seen
        seen.add(item["case"])
        assert item["status"] in statuses
        plan = item["plan"]
        n, q, t, k, m = (plan[key] for key in ("n", "q", "t", "k", "m"))
        assert m == n * ((q-1).bit_length()+1)
        assert plan["target_rows"] == n*t and plan["preimage_rows"] == 2*m*t
        if item["status"] == "skipped_dimension_budget":
            assert plan["preimage_rows"] > report["limits"]["max_dimension"]
            assert item["probe"] is None
            continue
        probe = item["probe"]
        if item["status"] != "completed":
            # Incomplete runs remain incomplete; a checkpoint is not success.
            assert not probe or probe["status"] != "completed"
            if not probe or not probe.get("algorithm_flow_completed"):
                continue
        else:
            assert item["returncode"] == 0 and probe["status"] == "completed"
            assert probe["protocol_gate_rejects"]
            assert "sampling_diagnostics" in probe and "public_attack" in probe
        assert probe["completed_trials"] == report["limits"]["trials_per_point"]
        p = probe["parameters"]
        assert all(p[key] == plan[key] for key in p)
        verify_profile(probe["master_sampling_profile"], n*t, 2*m*t, q)
        verify_profile(probe["auxiliary_sampling_profile"], n, m, q)
        width = Fraction(probe["sigma_squared"]["value"])
        assessment = probe["assessment"]
        beta = Fraction(assessment["accepted_full_norm_bound_squared"]["value"])
        bound = n*t*(q//2)**2
        assert beta == t*t*2*m*t*width
        assert int(assessment["public_linear_solution_norm_bound_squared"]["value"]) == bound
        assert assessment["accepts_public_linear_solutions"] == (beta >= bound)
        assert Fraction(assessment["table2_issuance_spectral_upper_bound_squared"]["value"]) == Fraction(q*q, m)
        initial = Fraction(probe["master_sampling_profile"]["initial_width"]["value"])
        assert assessment["table2_issuance_upper_bound_satisfied"] == (initial**2 <= Fraction(q*q, m))
        assert assessment["table2_m_greater_than_participants"] == (m > k)
        for flag in ("protocol_security_approved", "concrete_hardness_established",
                     "trapdoor_uniformity_established", "smoothing_conditions_established"):
            assert not assessment[flag]
        assert all(0 <= int(norm) <= beta for norm in probe["integer_output_norms_squared"])
        assert probe["algorithm_flow_completed"]
        a_online = matrix(probe["a_mod_q"], n*t, 2*m*t, q)
        assert len(probe["online_outputs"]) == probe["completed_trials"]
        for output, recorded in zip(probe["online_outputs"], probe["integer_output_norms_squared"]):
            y = matrix(output["target_mod_q"], n*t, 1, q)
            x = matrix(output["x_integer"], 2*m*t, 1)
            congruent(multiply(a_online, x), y, q)
            assert norm_squared(x) == int(recorded) <= beta
        if "sampling_diagnostics" in probe:
            d = probe["sampling_diagnostics"]
            assert d["pullback_metric_exact"] and 0 <= d["total_squared_width_relative_residual"] <= 1e-8
            assert all(math.isfinite(d[key]) and d[key] > 0 for key in
                       ("gram_condition_inf", "covariance_condition_inf", "perturbation_condition_inf"))
            assert d["perturbation_condition_inf"] * 2**-52 <= 1e-8
        for key, count in (("ptrap_gen", k), ("psamp_pre_per_share", t*probe["completed_trials"]),
                           ("rec", probe["completed_trials"]), ("online_group_sequential", probe["completed_trials"])):
            timings = probe[key]
            values = timings["values_ms"]
            assert timings["count"] == count == len(values)
            assert all(value >= 0 and math.isfinite(value) for value in values)
            assert math.isclose(timings["mean_ms"], sum(values)/count, rel_tol=1e-12)
        if "public_attack" not in probe:
            continue
        attack = probe["public_attack"]
        assert attack["secrets_dropped_before_attack"]
        assert attack["attempts"] == report["limits"]["public_attempts_per_point"]
        assert attack["solver_errors"] == []
        a = matrix(attack["a_mod_q"], n*t, 2*m*t, q)
        target = matrix(attack["target_mod_q"], n*t, 1, q)
        x = matrix(attack["minimum_norm_witness"], 2*m*t, 1)
        congruent(multiply(a, x), target, q)
        norm = norm_squared(x)
        assert norm == int(attack["minimum_norm_squared"]["value"]) <= bound
        assert 0 <= attack["accepted"] <= attack["attempts"]
        assert (attack["accepted"] > 0) == (norm <= beta)
    return len(seen)


def exhaustive_pairs(item):
    probe = item["probe"]
    p = probe["parameters"]
    q, columns = p["q"], p["preimage_rows"]
    public = probe["public_attack"]
    a = matrix(public["a_mod_q"], 2, columns, q)
    y = matrix(public["target_mod_q"], 2, 1, q)
    beta = Fraction(probe["assessment"]["accepted_full_norm_bound_squared"]["value"])
    best, witness = None, None
    invertible = accepted = 0
    start = time.monotonic()
    for i, j in combinations(range(columns), 2):
        determinant = (a[0][i]*a[1][j]-a[0][j]*a[1][i]) % q
        if determinant == 0:
            continue
        invertible += 1
        inverse = pow(determinant, -1, q)
        left = ((y[0][0]*a[1][j]-a[0][j]*y[1][0])*inverse) % q
        right = ((a[0][i]*y[1][0]-y[0][0]*a[1][i])*inverse) % q
        left = left-q if left > q//2 else left
        right = right-q if right > q//2 else right
        norm = left*left + right*right
        if norm <= beta:
            accepted += 1
        if best is None or norm < best:
            best, witness = norm, {"columns": [i, j], "integer_entries": [left, right]}
    # Independently reconstruct and verify the minimum witness.
    x = [[0] for _ in range(columns)]
    for column, entry in zip(witness["columns"], witness["integer_entries"]):
        x[column][0] = entry
    congruent(multiply(a, x), y, q)
    assert norm_squared(x) == best
    return {"case": item["case"], "total_pairs": columns*(columns-1)//2,
            "invertible_pairs": invertible, "accepted_pairs": accepted,
            "minimum_norm_squared": str(best), "beta_squared": str(beta),
            "minimum_witness": witness, "elapsed_seconds": time.monotonic()-start,
            "scope": "all two-column pivot solutions for one saved public A,y; does not cover arbitrary integer solutions"}


def verify_lattice_report(report, lattice):
    cases = {item["case"]: item for item in report["results"] if item["status"] == "completed"
             or (item.get("probe") or {}).get("algorithm_flow_completed")}
    assert lattice["security_bit_target"] is None
    assert {item["case"] for item in lattice["results"]} == set(cases)
    assert len(lattice["results"]) == len(cases)
    for result in lattice["results"]:
        item = cases[result["case"]]
        assert result["status"] == "completed", "public lattice screen error retained"
        p = item["plan"]
        rows, columns, q = p["target_rows"], p["preimage_rows"], p["q"]
        probe = item["probe"]
        if "public_attack" in probe:
            assert result["target_source"] == "public_attack"
            a_value, y_value = probe["public_attack"]["a_mod_q"], probe["public_attack"]["target_mod_q"]
        else:
            assert result["target_source"] == "online_trial_0"
            a_value, y_value = probe["a_mod_q"], probe["online_outputs"][0]["target_mod_q"]
        a = matrix(a_value, rows, columns, q)
        y = matrix(y_value, rows, 1, q)
        x = matrix(result["x_integer"], columns, 1)
        congruent(multiply(a, x), y, q)
        norm = norm_squared(x)
        assert norm == int(result["norm_squared"]["value"])
        beta = Fraction(item["probe"]["assessment"]["accepted_full_norm_bound_squared"]["value"])
        assert result["accepted_by_public_equation_and_norm"] == (norm <= beta)
        assert 1 <= result["support"] <= columns
        assert all(row == [0] for row in x[result["support"]:])
    return len(cases)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--exhaustive-pairs-output", type=Path)
    parser.add_argument("--lattice-screen", type=Path)
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    count = verify(report)
    print(f"Verified {count} retained results; no protocol parameter approval.")
    if args.lattice_screen:
        lattice = json.loads(args.lattice_screen.read_text())
        assert lattice["input_report_sha256"] == hashlib.sha256(args.report.read_bytes()).hexdigest()
        print(f"Verified {verify_lattice_report(report, lattice)} public lattice witnesses independently.")
    if args.exhaustive_pairs_output:
        results = [exhaustive_pairs(item) for item in report["results"]
                   if item["status"] == "completed" and item["plan"]["target_rows"] == 2]
        output = {"input_report": str(args.report), "security_bit_target": None, "results": results}
        args.exhaustive_pairs_output.write_text(json.dumps(output, indent=2)+"\n")
        for result in results:
            print(f"{result['case']}: {result['accepted_pairs']}/{result['invertible_pairs']} public pairs accepted")


if __name__ == "__main__":
    main()

"""Verify saved B audit evidence; a valid report still rejects B3 protocol parameters.

Usage: python3 scripts/verify_sampling_audit.py target/plt_sampling_B.json
"""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from verify_plt_report import matrix, multiply, congruent, norm_squared


def verify(report):
    assert len(report["numerical_diagnostics"]) == 54
    assert len(report["plt_cohorts"]) == 27 and report["sampling_attempts"] == 27000
    assert report["errors"] == []
    for item in report["numerical_diagnostics"]:
        metrics = item["metrics"]
        assert metrics["pullback_metric_exact"]
        assert 0 <= metrics["total_squared_width_relative_residual"] <= 1e-8
        for key in ("gram_condition_inf", "covariance_condition_inf", "perturbation_condition_inf"):
            assert math.isfinite(metrics[key]) and metrics[key] > 0
        assert metrics["perturbation_condition_inf"] * 2**-52 <= 1e-8
    for item in report["plt_cohorts"]:
        m = item["metrics"]
        assert m["samples"] == 1000 and m["duplicate_complete_vectors"] == 0
        assert all(math.isfinite(x) for x in m["whitened_means"])
        assert all(0.7 <= x <= 1.3 for x in m["whitened_variances"])
        assert 0 <= m["max_abs_coordinate_correlation"] <= 0.2
        assert 0 <= m["max_abs_lag1_correlation"] <= 0.2
        if m["zero_target_mean_check_applied"]:
            assert max(map(abs, m["whitened_means"])) <= 0.15
    for item in report["b3_counterexamples"]:
        t, k = item["t"], item["k"]
        a = matrix(item["a_mod_q"], t, 20*t, 257)
        y = matrix(item["target_mod_q"], t, 1, 257)
        x = matrix(item["public_integer_solution"], 20*t, 1)
        congruent(multiply(a, x), y, 257)
        norm = norm_squared(x)
        assessment = item["assessment"]
        assert norm <= t*128**2
        assert norm <= Fraction(assessment["accepted_full_norm_bound_squared"]["value"])
        assert assessment["accepts_public_linear_solutions"] and not assessment["protocol_security_approved"]
        assert Fraction(assessment["table2_issuance_spectral_upper_bound_squared"]["value"]) == Fraction(257**2, 10)
        assert not assessment["table2_issuance_upper_bound_satisfied"]
        assert assessment["table2_m_greater_than_participants"] == (k < 10)
        assert item["public_solution_accepted_by_preimage_checks"] and item["parameter_gate_rejected"]
    reference = report["small_discrete_reference"]
    assert reference["modulus"] == 3
    for item in reference["reports"]:
        partition = item["partition_function_in_box"]
        limit = item["reference_box_limit"]
        tail = 6*6800/(2*math.pi*limit)*math.exp(-math.pi*limit**2/6800)*(1+math.sqrt(6800))**2/partition
        assert tail < 1e-8 and math.isclose(tail, item["relative_omitted_mass_upper_bound"], rel_tol=1e-12)
        assert item["empirical_projection_outside_count"] == 0
        for weights, counts, recorded in zip(item["reference_projection_weights"], item["empirical_projection_counts"], item["ks_for_d1_d2_d3_d1_plus_d2"]):
            assert len(weights) == len(counts) == 4*limit+1
            assert all(w >= 0 and math.isfinite(w) for w in weights)
            assert sum(counts) == 10000
            assert math.isclose(sum(weights), partition, rel_tol=1e-7)
            cdf = empirical = difference = 0.0
            for weight, count in zip(weights, counts):
                cdf += weight/partition
                empirical += count/10000
                difference = max(difference, abs(cdf-empirical))
            assert difference <= 0.04 and math.isclose(difference, recorded, abs_tol=1e-12)
    assert report["engineering_screening_passed"] and not report["protocol_parameters_passed"]
    assert not report["stage_b_fully_passed"]
    print("Verified B2/B4 finite screening evidence and 3 public-solution counterexamples; B3 remains REJECTED.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    verify(json.loads(args.report.read_text(encoding="utf-8")))

"""Independently verify a configurable PLT JSON report using Python integer arithmetic.

Usage: python3 scripts/verify_plt_report.py target/plt_2_of_3.json
No third-party dependencies. Checks algebra and norm bounds, not distributions.
"""

import argparse
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path


def matrix(encoded, rows, columns, modulus=None):
    """Decode qfall's matrix serialization and check dimensions and modulus."""
    raw, separator, mod = encoded["matrix"].partition(" mod ")
    assert bool(separator) == (modulus is not None), "unexpected matrix type"
    if modulus is not None:
        assert int(mod) == modulus, "wrong modulus"
    result = json.loads(raw)
    assert len(result) == rows and all(len(row) == columns for row in result)
    assert all(type(value) is int for row in result for value in row)
    return result


def multiply(left, right):
    """Multiply matrices over the integers without reducing entries."""
    assert len(left[0]) == len(right)
    return [[sum(a * b for a, b in zip(row, column)) for column in zip(*right)]
            for row in left]


def congruent(left, right, modulus):
    """Assert equality of equally sized matrices modulo the given modulus."""
    assert len(left) == len(right)
    assert all(len(a) == len(b) for a, b in zip(left, right))
    assert all((a - b) % modulus == 0
               for row_a, row_b in zip(left, right) for a, b in zip(row_a, row_b))


def image(point, value, threshold, modulus):
    """Compute (1, point, ..., point^(t-1))^T tensor an n-row matrix."""
    return [[pow(point, block, modulus) * entry for entry in row]
            for block in range(threshold) for row in value]


def norm_squared(vector):
    """Compute the squared Euclidean norm of an integer column vector."""
    assert all(len(row) == 1 for row in vector)
    return sum(row[0] ** 2 for row in vector)


def verify_profile(profile, rows, columns, modulus):
    """Check the non-secret construction metadata when present in a new report."""
    ell = (modulus-1).bit_length()
    bar = columns-rows*ell
    assert profile["rows"] == rows and profile["columns"] == columns
    assert profile["gadget_length"] == ell and profile["m_bar"] == bar
    assert profile["rounding_parameter"] == 4
    assert profile["distribution"] == "MP12 PlusMinusOneZero: P(-1)=1/4, P(0)=1/2, P(1)=1/4"
    width = Fraction(profile["initial_width"]["value"])
    assert width.denominator == 1 and width.numerator % 4 == 0
    scale = width.numerator//4
    assert scale >= 64 and scale & (scale-1) == 0
    assert scale**2 > 5*rows*ell*(bar+1)+1


def verify(report):
    """Check all four congruences and exact integer products/sums in the report."""
    parameters = report["parameters"]
    q, n, m, t, k = (parameters[key] for key in ("q", "n", "m", "t", "k"))
    assert n > 0 and q > 1 and 1 <= t <= k < q
    assert m == n * ((q-1).bit_length()+1)
    rows, preimage_rows = n*t, 2*m*t
    trials = report.get("trials_per_set", 8)  # Also accept the historical 2-of-3 report.
    assert type(trials) is int and trials > 0
    a = matrix(report["a_mod_q"], rows, preimage_rows, q)
    if "master_sampling_profile" in report:
        verify_profile(report["master_sampling_profile"], rows, preimage_rows, q)
    width = Fraction(report["sigma_squared"]["value"])
    beta_squared = Fraction(report["beta_squared"]["value"])
    assert width > 0 and beta_squared == t*t*preimage_rows*width
    partials = {}
    for partial in report["partial_trapdoors"]:
        participant = partial["participant"]
        assert participant in range(k) and participant not in partials
        assert partial["point"] == participant + 1
        if "auxiliary_sampling_profile" in partial:
            verify_profile(partial["auxiliary_sampling_profile"], n, m, q)
        u = matrix(partial["u_integer"], preimage_rows, m)
        c = matrix(partial["c_mod_q"], n, m, q)
        congruent(multiply(a, u), image(participant + 1, c, t, q), q)
        partials[participant] = (u, c)
    assert set(partials) == set(range(k))
    coverage = {pair: set() for pair in combinations(range(k), t)}
    targets = {}
    for record in report["reconstructions"]:
        pair = tuple(record["recovery_set"])
        assert pair in coverage
        trial = record["trial"]
        assert trial in range(trials) and trial not in coverage[pair]
        coverage[pair].add(trial)
        y = matrix(record["target_mod_q"], rows, 1, q)
        assert targets.setdefault(trial, y) == y
        xs, ids, images = [], [], []
        for share in record["shares"]:
            participant = share["participant"]
            ids.append(participant)
            u, c = partials[participant]
            z = matrix(share["z_mod_q"], n, 1, q)
            d = matrix(share["d_integer"], m, 1)
            x = matrix(share["x_integer"], preimage_rows, 1)
            congruent(multiply(c, d), z, q)
            assert multiply(u, d) == x, "x_j must be the integer product U_j*d_j"
            vz = image(participant + 1, z, t, q)
            congruent(multiply(a, x), vz, q)
            norm = norm_squared(x)
            assert norm == int(share["norm_squared"]["value"]) and norm <= preimage_rows * width
            xs.append(x)
            images.append(vz)
        assert sorted(ids) == list(pair)
        congruent([[sum(value[i][0] for value in images)] for i in range(rows)], y, q)
        x = matrix(record["x_integer"], preimage_rows, 1)
        assert x == [[sum(value[i][0] for value in xs)] for i in range(preimage_rows)]
        congruent(multiply(a, x), y, q)
        norm = norm_squared(x)
        assert norm == int(record["norm_squared"]["value"]) and norm <= beta_squared
    assert all(observed == set(range(trials)) for observed in coverage.values())
    print(f"Verified {len(report['reconstructions'])} reconstructions for {t}-of-{k}, "
          "all 4 congruences, integer products/sums and norm bounds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    verify(json.loads(args.report.read_text(encoding="utf-8")))

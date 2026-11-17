"""Known-answer and tampering tests for public parameter evidence checks."""

import copy
import unittest
from verify_plt_parameter_exploration import exhaustive_pairs, verify_lattice_report


def encoded(rows, modulus=None):
    import json
    return {"matrix": json.dumps(rows) + (f" mod {modulus}" if modulus else "")}


def fixture():
    return {"case": "known", "status": "completed",
            "plan": {"target_rows": 2, "preimage_rows": 3, "q": 257},
            "probe": {"parameters": {"q": 257, "preimage_rows": 3},
                      "public_attack": {"a_mod_q": encoded([[1, 0, 2], [0, 1, 3]], 257),
                                        "target_mod_q": encoded([[1], [1]], 257)},
                      "assessment": {"accepted_full_norm_bound_squared": {"value": "2"}}}}


class Screens(unittest.TestCase):
    def test_pair_enumeration_known_answer(self):
        result = exhaustive_pairs(fixture())
        self.assertEqual((result["total_pairs"], result["invertible_pairs"]), (3, 3))
        self.assertEqual(result["accepted_pairs"], 1)
        self.assertEqual(int(result["minimum_norm_squared"]), 2)

    def test_witness_congruence_and_acceptance_are_recomputed(self):
        report = {"results": [fixture()]}
        lattice = {"security_bit_target": None, "results": [
            {"case": "known", "status": "completed", "target_source": "public_attack",
             "support": 3, "norm_squared": {"value": "2"},
             "x_integer": encoded([[1], [1], [0]]), "accepted_by_public_equation_and_norm": True}]}
        self.assertEqual(verify_lattice_report(report, lattice), 1)
        tampered = copy.deepcopy(lattice)
        tampered["results"][0]["accepted_by_public_equation_and_norm"] = False
        with self.assertRaises(AssertionError):
            verify_lattice_report(report, tampered)
        tampered = copy.deepcopy(lattice)
        tampered["results"][0]["x_integer"] = encoded([[0], [1], [0]])
        with self.assertRaises(AssertionError):
            verify_lattice_report(report, tampered)


if __name__ == "__main__":
    unittest.main()

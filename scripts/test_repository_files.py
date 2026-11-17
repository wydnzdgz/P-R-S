"""Credential detector and ignore-rule regressions; fixtures are synthetic."""

import unittest
from check_repository_files import content_findings, ignored_paths


class PublicationPolicy(unittest.TestCase):
    def test_private_keys_and_provider_credentials_are_detected_without_echoing_values(self):
        cases = [
            (b"-----BEGIN "+b"OPENSSH PRIVATE KEY-----", "private key"),
            (b"ghp_"+b"a"*36, "GitHub credential"),
            (b"github_pat_"+b"a"*60, "GitHub credential"),
            (b"AKIA"+b"A"*16, "AWS access key"),
            (b"https"+b"://user:"+b"private-value@github.com/repo", "credential in URL"),
        ]
        for fixture, category in cases:
            with self.subTest(category=category):
                self.assertIn(category, content_findings(fixture))
                self.assertNotIn(fixture.decode(), str(content_findings(fixture)))

    def test_public_experiment_json_and_ci_secret_references_are_allowed(self):
        self.assertEqual(content_findings(b'{"sigma_squared":{"value":"19200216832"}}'), [])
        self.assertEqual(content_findings(b"${{ secrets.GITHUB_TOKEN }}"), [])

    def test_ignore_rules_keep_source_lockfile_and_reviewed_evidence(self):
        protected = {".env", ".env.local", "nested/.env.production", ".git-credentials",
                     "credentials.json", "keys/id_ed25519", "keys/server.pem", "paper.pdf",
                     "archive.tar.gz", "target/debug/tool", "tmp/test.json", "trace.log",
                     ".migration-backups/checkpoint.tar.gz", "experiments/raw/capture.json"}
        retained = {"Cargo.lock", ".env.example", ".env.sample", ".env.template",
                    "docs/plt_sampling_B_results.json", "crates/plt/src/lib.rs",
                    "experiments/config.toml"}
        self.assertEqual(ignored_paths(sorted(protected | retained)), protected)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Check package ownership and one-way dependencies after the Workspace split.

--migration-snapshot additionally checks current files against the recorded
migration snapshot; later deliberate source changes require a new snapshot.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tomllib

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--migration-snapshot", action="store_true")
    args = parser.parse_args()
    manifest = tomllib.loads((ROOT / "Cargo.toml").read_text())
    assert "workspace" in manifest and "package" not in manifest
    metadata = json.loads(subprocess.check_output(
        ["cargo", "metadata", "--locked", "--offline", "--no-deps", "--format-version", "1"],
        cwd=ROOT, text=True))
    members = set(metadata["workspace_members"])
    packages = {p["name"]: p for p in metadata["packages"] if p["id"] in members}
    for name in ("qfall-tools", "plt"):
        assert name in packages
        assert Path(packages[name]["manifest_path"]) == ROOT / "crates" / name / "Cargo.toml"
    core_dependencies = {d["name"] for d in packages["qfall-tools"]["dependencies"]}
    assert not core_dependencies.intersection(packages), "baseline must not depend on an upper layer"
    plt_dependencies = {d["name"] for d in packages["plt"]["dependencies"] if d["kind"] is None}
    assert "qfall-tools" in plt_dependencies
    assert plt_dependencies.intersection(packages) == {"qfall-tools"}
    for old in ("src", "tests", "examples", "benches"):
        assert not (ROOT / old).exists(), f"legacy root directory remains: {old}"
    assert not (ROOT / "crates/qfall-tools/src/primitive/plt.rs").exists()
    for source in (ROOT / "crates").rglob("*.rs"):
        assert "qfall_tools::primitive::plt" not in source.read_text(), str(source)
    snapshot_path = ROOT / "docs/workspace_migration_manifest.json"
    if snapshot_path.exists():
        snapshot = json.loads(snapshot_path.read_text())
        for name, checksum in snapshot["historical_artifacts_sha256"].items():
            assert digest(ROOT / name) == checksum, f"historical evidence modified: {name}"
        if args.migration_snapshot:
            for entry in snapshot["file_mappings"]:
                assert digest(ROOT / entry["new_path"]) == entry["after_sha256"], entry["new_path"]
            for name, checksum in snapshot["migration_sources_sha256"].items():
                assert digest(ROOT / name) == checksum, name
            before = {tuple(item) for item in snapshot["external_dependency_versions_before"]}
            lock = tomllib.loads((ROOT / "Cargo.lock").read_text())
            after = {(p["name"], p["version"], p["source"]) for p in lock["package"] if "source" in p}
            assert before == after, "external dependency versions changed during migration"
    else:
        assert not args.migration_snapshot, "migration manifest not available"
    print("Workspace layout checked: plt -> qfall-tools; historical artifacts preserved.")


if __name__ == "__main__":
    main()

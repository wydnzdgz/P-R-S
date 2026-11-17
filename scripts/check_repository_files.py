#!/usr/bin/env python3
"""Reject ignored paths, blobs over 1 MiB and recognizable credentials.

Default: check the Git index (the proposed commit), rather than untracked local
files. --history checks all blobs reachable from HEAD, including earlier commits.
Never print matched credential contents. This is a finite screen, not a guarantee
that all forms of secrets can be recognized automatically.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MAX_BYTES = 1024 * 1024
PATTERNS = {
    "private key": rb"-----BEGIN (?:(?:RSA|EC|DSA|OPENSSH|ENCRYPTED) )?PRIVATE KEY-----",
    "GitHub credential": rb"(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,})",
    "AWS access key": rb"\b(?:AKIA|ASIA)[0-9A-Z]{16}\b",
    "Slack credential": rb"\bxox[baprs]-[A-Za-z0-9-]{20,}",
    "API key": rb"\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{32,}",
    "credential in URL": rb"[a-z][a-z0-9+.-]*://[^\s/@:]+:[^\s/@]+@",
}
COMPILED = {name: re.compile(pattern) for name, pattern in PATTERNS.items()}


@dataclass
class Blob:
    oid: str
    path: str
    size: int


def git(*args, input_data=None):
    return subprocess.check_output(["git", *args], cwd=ROOT, input=input_data)


def content_findings(content):
    return [name for name, pattern in COMPILED.items() if pattern.search(content)]


def blobs(history):
    if history:
        objects = git("rev-list", "--objects", "HEAD")
        info = git("cat-file", "--batch-check=%(objectname) %(objecttype) %(objectsize) %(rest)",
                   input_data=objects)
        result = []
        for line in info.decode(errors="surrogateescape").splitlines():
            fields = line.split(" ", 3)
            if len(fields) == 4 and fields[1] == "blob":
                result.append(Blob(fields[0], fields[3], int(fields[2])))
        return result
    entries = []
    for record in git("ls-files", "--stage", "-z").split(b"\0"):
        if not record:
            continue
        header, path = record.split(b"\t", 1)
        mode, oid, stage = header.decode().split()
        if stage != "0" or mode == "160000":
            raise ValueError("resolve index conflicts/submodule ownership before publication")
        entries.append((oid, path.decode(errors="surrogateescape")))
    info = git("cat-file", "--batch-check=%(objectname) %(objectsize)",
               input_data="".join(f"{oid}\n" for oid, _ in entries).encode())
    sizes = {line.split()[0]: int(line.split()[1]) for line in info.decode().splitlines()}
    return [Blob(oid, path, sizes[oid]) for oid, path in entries]


def ignored_paths(paths):
    if not paths:
        return set()
    result = subprocess.run(["git", "check-ignore", "--no-index", "-z", "--stdin"], cwd=ROOT,
                            input=b"\0".join(p.encode(errors="surrogateescape") for p in paths)+b"\0",
                            stdout=subprocess.PIPE, check=False)
    if result.returncode not in (0, 1):
        raise RuntimeError("git check-ignore failed")
    return {p.decode(errors="surrogateescape") for p in result.stdout.split(b"\0") if p}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history", action="store_true")
    parser.add_argument("--max-bytes", type=int, default=MAX_BYTES)
    args = parser.parse_args()
    if args.max_bytes <= 0:
        parser.error("max-bytes must be positive")
    entries = blobs(args.history)
    ignored = ignored_paths([entry.path for entry in entries])
    findings, content_cache = [], {}
    for entry in entries:
        if entry.path in ignored:
            findings.append((entry.path, "path excluded by .gitignore"))
        if entry.size > args.max_bytes:
            findings.append((entry.path, f"{entry.size} bytes exceeds {args.max_bytes}-byte limit"))
            continue
        if entry.oid not in content_cache:
            content_cache[entry.oid] = content_findings(git("cat-file", "blob", entry.oid))
        for category in content_cache[entry.oid]:
            findings.append((entry.path, f"possible {category}"))
    if findings:
        for path, reason in sorted(set(findings)):
            print(f"REJECT: {path}: {reason}")
        raise SystemExit(1)
    largest = max((entry.size for entry in entries), default=0)
    print(f"Checked {len(entries)} {'reachable historical' if args.history else 'indexed'} blobs; "
          f"largest {largest} bytes; no configured policy findings.")


if __name__ == "__main__":
    main()

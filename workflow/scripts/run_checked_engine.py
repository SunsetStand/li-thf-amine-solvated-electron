#!/usr/bin/env python3
"""Run an external engine and refuse to mark unconverged output successful."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

from solvelec.parsers import parse_output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", choices=["cp2k", "orca"], required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cwd", default=".")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("an engine command is required after --")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(command, cwd=args.cwd, stdout=handle, stderr=subprocess.STDOUT)
    result = parse_output(args.engine, output)
    if completed.returncode != 0 or not result.converged:
        print(f"Engine validation failed: {result.as_dict()}", file=sys.stderr)
        return completed.returncode or 4
    print(f"Validated: {shlex.join(command)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

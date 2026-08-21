from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path

from .engines import doctor_report


def sha256_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _git(root: Path, *args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *args], cwd=root, capture_output=True, text=True, check=True, timeout=10
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout.strip()


def git_state(root: Path) -> dict[str, object]:
    status = _git(root, "status", "--porcelain")
    return {
        "commit": _git(root, "rev-parse", "HEAD"),
        "branch": _git(root, "branch", "--show-current"),
        "dirty": None if status is None else bool(status),
        "remote": _git(root, "remote", "get-url", "origin"),
    }


def build_manifest(
    root: Path, inputs: Iterable[str | Path] = (), campaign: str | None = None
) -> dict[str, object]:
    input_records = []
    for path_value in inputs:
        path = Path(path_value).resolve()
        input_records.append(
            {"path": str(path), "size": path.stat().st_size, "sha256": sha256_file(path)}
        )
    return {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "campaign": campaign,
        "git": git_state(root),
        "runtime": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "hostname": platform.node(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        },
        "capabilities": doctor_report(),
        "inputs": input_records,
    }


def write_manifest(
    path: str | Path,
    root: Path,
    inputs: Iterable[str | Path] = (),
    campaign: str | None = None,
) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(build_manifest(root, inputs, campaign), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

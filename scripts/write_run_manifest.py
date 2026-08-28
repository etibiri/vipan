#!/usr/bin/env python3
"""Record input hashes and software versions for a completed workflow run."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Sequence


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def labeled(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected LABEL=VALUE")
    label, item = value.split("=", 1)
    if not label or not item:
        raise argparse.ArgumentTypeError("expected non-empty LABEL=VALUE")
    return label, item


def command_version(command: str) -> Dict[str, object]:
    try:
        result = subprocess.run(
            shlex.split(command),
            check=False,
            capture_output=True,
            text=True,
            timeout=20,
        )
        output = (result.stdout or result.stderr).strip().splitlines()
        return {
            "command": command,
            "exit_code": result.returncode,
            "first_output_line": output[0] if output else "",
        }
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"command": command, "error": str(error)}


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write a reproducibility manifest.")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--file", action="append", default=[], type=labeled, metavar="LABEL=PATH")
    parser.add_argument("--version-command", action="append", default=[], type=labeled, metavar="LABEL=COMMAND")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        files = {}
        for label, value in args.file:
            path = Path(value)
            if not path.is_file():
                raise FileNotFoundError(path)
            files[label] = {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        versions = {label: command_version(command) for label, command in args.version_command}
        manifest = {
            "run_id": args.run_id,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "python": sys.version.replace("\n", " "),
            "platform": platform.platform(),
            "files": files,
            "software": versions,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Run manifest: {args.output}")
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

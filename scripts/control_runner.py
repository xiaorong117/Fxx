#!/usr/bin/env python3
"""Hold a case flock and exec the fixed bubble solver without a shell."""
from __future__ import annotations

import argparse
import fcntl
import os
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--production-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.production_root.resolve(); output = args.output.resolve()
    if output.parent != root or not output.is_dir():
        raise SystemExit("output is not a direct production child")
    for path in (args.binary, args.config, args.checkpoint):
        if not path.resolve().is_file():
            raise SystemExit(f"required file missing: {path}")
    lock_fd = os.open(output / "run.lock", os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("case run.lock is already held")
    os.set_inheritable(lock_fd, True)
    log_fd = os.open(output / "driver.log", os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    os.dup2(log_fd, 1); os.dup2(log_fd, 2)
    os.chdir(Path(__file__).resolve().parents[2])
    binary = str(args.binary.resolve())
    os.execv(binary, [binary, "--config", str(args.config.resolve()),
                      "--restart", str(args.checkpoint.resolve())])


if __name__ == "__main__":
    main()

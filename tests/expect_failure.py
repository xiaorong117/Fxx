#!/usr/bin/env python3
import subprocess
import sys

result = subprocess.run(sys.argv[1:-1], text=True, capture_output=True)
combined = result.stdout + result.stderr
if result.returncode == 0 or sys.argv[-1] not in combined:
    print(combined)
    raise SystemExit("expected failure or diagnostic was not observed")
print(combined.strip())


import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOG_DIR = PROJECT_ROOT / "logs"

print(f"Root: {PROJECT_ROOT}")
print(f"Log Dir: {LOG_DIR}")
print(f"Log Dir Exists: {LOG_DIR.exists()}")

test_log = LOG_DIR / "test_debug.log"
try:
    with open(test_log, "w") as f:
        f.write("Test log content")
    print(f"Wrote to {test_log}")
except Exception as e:
    print(f"Error writing: {e}")

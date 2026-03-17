
import sys
from pathlib import Path
print(f"Python executable: {sys.executable}")
print(f"Current working directory: {Path.cwd()}")
try:
    PROJECT_ROOT = Path.cwd()
    sys.path.insert(0, str(PROJECT_ROOT))
    from scripts.utils.logger import print_status
    print_status("Logger imported successfully", "TITLE")
except Exception as e:
    print(f"Import failed: {e}")

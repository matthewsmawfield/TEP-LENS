import datetime
import sys


def print_status(message, level="INFO"):
    timestamp = datetime.datetime.now().strftime("%H:%M:%S")
    if level == "TITLE":
        print(f"\n[{timestamp}] [INFO] {'='*70}")
        print(f"[{timestamp}] [INFO] {message}")
        print(f"[{timestamp}] [INFO] {'='*70}")
    elif level == "WARNING":
        print(f"[{timestamp}] [WARN] {message}")
    elif level == "ERROR":
        print(f"[{timestamp}] [FAIL] {message}")
        sys.exit(1)
    else:
        print(f"[{timestamp}] [INFO] {message}")


def safe_json_default(obj):
    """JSON serializer for numpy types and other non-standard objects."""
    if hasattr(obj, 'item'):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

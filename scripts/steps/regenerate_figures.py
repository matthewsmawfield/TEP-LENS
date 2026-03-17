#!/usr/bin/env python3
"""
Regenerate all 9 manuscript figures with unified styling.
"""

import subprocess
import sys
from pathlib import Path

# The 9 figures we're keeping (in order)
FIGURES_TO_REGENERATE = [
    ("step_04_plot_closure.py", "Figure 1: Route-closure residuals"),
    ("step_07_observed_vs_predicted.py", "Figures 2-3: Observed vs predicted, residuals"),
    ("step_08_new_evidence.py", "Figures 4-5: Delay vs mu, alpha inference"),
    ("step_10_h0_tension.py", "Figure 8: H0 tension resolution"),
    ("step_11_model_dependence.py", "Figure 6: Model dependence"),
    ("step_16_independence_tier_significance.py", "Figure 7: Significance synthesis"),
]

def run_script(script_name):
    """Run a figure generation script."""
    script_path = Path(__file__).parent / script_name
    
    if not script_path.exists():
        print(f"✗ Script not found: {script_name}")
        return False
    
    print(f"\n{'='*60}")
    print(f"Running: {script_name}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=Path(__file__).parent.parent.parent,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0:
            print(f"✓ Success: {script_name}")
            if result.stdout:
                print(result.stdout)
            return True
        else:
            print(f"✗ Failed: {script_name}")
            if result.stderr:
                print(f"Error: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"✗ Timeout: {script_name}")
        return False
    except Exception as e:
        print(f"✗ Exception: {script_name} - {e}")
        return False

def main():
    """Regenerate all figures."""
    print("Regenerating manuscript figures with unified styling...")
    print(f"Total figures to regenerate: {len(FIGURES_TO_REGENERATE)}")
    
    results = []
    for script, description in FIGURES_TO_REGENERATE:
        print(f"\n{description}")
        success = run_script(script)
        results.append((script, description, success))
    
    # Summary
    print("\n" + "="*60)
    print("REGENERATION SUMMARY")
    print("="*60)
    
    success_count = sum(1 for _, _, success in results if success)
    
    for script, description, success in results:
        status = "✓" if success else "✗"
        print(f"{status} {description}")
    
    print(f"\nTotal: {success_count}/{len(results)} successful")
    
    if success_count == len(results):
        print("\n✓ All figures regenerated successfully!")
        return 0
    else:
        print(f"\n✗ {len(results) - success_count} figure(s) failed to regenerate")
        return 1

if __name__ == "__main__":
    sys.exit(main())

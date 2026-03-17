import os
import glob
import subprocess
from PIL import Image

def clean_and_run():
    print("Cleaning figures...")
    for p in glob.glob("results/figures/*.png") + glob.glob("site/dist/figures/*.png"):
        try:
            os.remove(p)
        except Exception as e:
            pass

    print("Running pipeline...")
    subprocess.run(["bash", "regenerate_all_figures.sh"], check=True)

    print("Inspecting figures...")
    figures = glob.glob("results/figures/*.png")
    for f in sorted(figures):
        try:
            with Image.open(f) as img:
                print(f"{os.path.basename(f)}: {img.size}")
        except Exception as e:
            print(f"Error reading {f}: {e}")

if __name__ == "__main__":
    clean_and_run()

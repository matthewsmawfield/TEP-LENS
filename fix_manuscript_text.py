
import sys
from pathlib import Path

def fix_file(file_path):
    path = Path(file_path)
    if not path.exists():
        print(f"File not found: {path}")
        return

    content = path.read_text(encoding='utf-8')
    
    # Target string segment to identify the sentence
    target_segment = "explanation for why both lensed supernova systems may be biased in the same direction"
    replacement_segment = "explanation for why both lensed supernova systems may be biased away from concordance"
    
    if target_segment in content:
        new_content = content.replace(target_segment, replacement_segment)
        path.write_text(new_content, encoding='utf-8')
        print(f"Successfully updated {file_path}")
    else:
        print(f"Target string not found in {file_path}")
        # Debug finding
        idx = content.find("explanation for why both lensed supernova systems")
        if idx != -1:
            print(f"Found context: ...{content[idx:idx+100]}...")

if __name__ == "__main__":
    fix_file("manuscripts/14manuscript-tep-lens.md")
    fix_file("site/components/4_results.html")

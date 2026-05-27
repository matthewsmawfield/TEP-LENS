#!/usr/bin/env python3
"""
Generate PDF from TEP-LENS site
================================

Generates a high-quality PDF from the built static HTML.
Uses html_to_pdf.py with settings aligned to other TEP papers (e.g. TEP-LLR).

Usage:
    python scripts/generate_site_pdf.py
    python scripts/generate_site_pdf.py --quality high --wait-time 5
"""

import asyncio
import argparse
import shutil
import subprocess
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / 'utils'))

from html_to_pdf import HTMLToPDFConverter, create_preset_configs

try:
    import yaml
except ImportError:
    yaml = None

PAPER_ID = '19'
PROJECT_CODE = 'TEP-LENS'
DEFAULT_CODENAME = 'Lisboa'


def load_citation_metadata():
    """Load version and codename from CITATION.cff."""
    base_dir = Path(__file__).parent.parent
    citation_file = base_dir / 'CITATION.cff'

    if not citation_file.exists():
        print("⚠️  CITATION.cff not found, using defaults")
        return {'version': '0.1', 'codename': DEFAULT_CODENAME, 'title': PROJECT_CODE}

    try:
        if yaml:
            with open(citation_file, 'r') as f:
                data = yaml.safe_load(f)
            version_str = data.get('version', 'v0.1')
        else:
            with open(citation_file, 'r') as f:
                content = f.read()
            version_match = re.search(r'version:\s*"?([^"\n]+)"?', content)
            version_str = version_match.group(1).strip() if version_match else 'v0.1'

        pattern = r'^(v?[\d.]+)(?:\s*\(([^)]+)\))?$'
        match = re.match(pattern, version_str.strip())

        if match:
            version = match.group(1).lstrip('v')
            codename = match.group(2) or DEFAULT_CODENAME
        else:
            version = version_str.lstrip('v')
            codename = DEFAULT_CODENAME

        return {'version': version, 'codename': codename, 'title': PROJECT_CODE}

    except Exception as e:
        print(f"⚠️  Error parsing CITATION.cff: {e}, using defaults")
        return {'version': '0.1', 'codename': DEFAULT_CODENAME, 'title': PROJECT_CODE}


def pdf_basename(metadata):
    version_str = f"v{metadata['version']}-{metadata['codename']}"
    return f"{PAPER_ID}-{PROJECT_CODE}-{version_str}"


def build_static_site():
    """Build the static site first."""
    print("🔨 Building static site...")
    site_dir = Path(__file__).parent.parent / 'site'
    build_script = site_dir / 'build.js'

    if not build_script.exists():
        print("❌ Build script not found:", build_script)
        return False

    try:
        result = subprocess.run(
            ['node', str(build_script)],
            cwd=site_dir,
            check=True,
            capture_output=True,
            text=True,
        )
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print("❌ Site build failed:")
        print(e.stderr)
        return False


def copy_pdf_to_docs(source_pdf: Path, docs_dir: Path):
    """Copy PDF to site/public/docs with project naming convention."""
    docs_dir.mkdir(parents=True, exist_ok=True)
    metadata = load_citation_metadata()
    target_name = f"{pdf_basename(metadata)}.pdf"
    target_path = docs_dir / target_name
    shutil.copy2(source_pdf, target_path)
    print(f"📄 Copied to: {target_path}")
    print(f"   Size: {target_path.stat().st_size / (1024 * 1024):.2f} MB")
    return target_path


def copy_pdf_to_root(source_pdf: Path, base_dir: Path):
    """Copy PDF to the project root directory."""
    metadata = load_citation_metadata()
    target_name = f"{pdf_basename(metadata)}.pdf"
    target_path = base_dir / target_name
    shutil.copy2(source_pdf, target_path)
    print(f"📄 Copied to root: {target_path}")
    print(f"   Size: {target_path.stat().st_size / (1024 * 1024):.2f} MB")
    return target_path


def process_pdf_with_metadata(pdf_path: Path):
    """Run the PDF processing script to add metadata and compress."""
    process_script = Path(__file__).parent / 'utils' / 'process_pdf.py'

    if not process_script.exists():
        print("⚠️  PDF processing script not found, skipping metadata embedding")
        return

    print("🔧 Processing PDF with metadata...")
    try:
        subprocess.run(
            [sys.executable, str(process_script), str(pdf_path), '--quality', 'ebook'],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"⚠️  PDF processing failed: {e}")


async def generate_pdf(quality: str = 'high', wait_time: float = 5.0, skip_build: bool = False):
    """Generate PDF from the TEP-LENS site."""
    if not skip_build:
        if not build_static_site():
            print("❌ Cannot generate PDF - site build failed")
            return False

    base_dir = Path(__file__).parent.parent
    dist_dir = base_dir / 'site' / 'dist'
    docs_dir = base_dir / 'site' / 'public' / 'docs'

    if not dist_dir.exists():
        print(f"❌ Dist directory not found: {dist_dir}")
        print("   Run site build first: node site/build.js")
        return False

    html_file = dist_dir / 'index.html'
    if not html_file.exists():
        print(f"❌ HTML file not found: {html_file}")
        return False

    output_pdf = dist_dir / 'manuscript.pdf'
    presets = create_preset_configs()

    if quality == 'maximum':
        options = presets['high_quality'].copy()
        options['scale'] = 0.72
        options['device_scale_factor'] = 3.0
        options['viewport'] = {'width': 1920, 'height': 1080}
        options['prefer_css_page_size'] = True
    elif quality == 'high':
        options = presets['high_quality'].copy()
        options['scale'] = 0.72
        options['device_scale_factor'] = 2.5
        options['viewport'] = {'width': 1920, 'height': 1080}
        options['prefer_css_page_size'] = True
    elif quality == 'print':
        options = presets['print_ready'].copy()
        options['scale'] = 0.72
        options['device_scale_factor'] = 2.0
    else:
        options = presets['web_optimized'].copy()

    options['wait_time'] = wait_time
    options['format'] = 'A4'
    options['margin_top'] = '1.2cm'
    options['margin_bottom'] = '1.5cm'
    options['margin_left'] = '1cm'
    options['margin_right'] = '1cm'
    options['display_header_footer'] = True
    options['header_template'] = '<div></div>'
    options['footer_template'] = '''
        <div style="font-size:9px; text-align:center; width:100%; color:#555555; font-family:system-ui,-apple-system,sans-serif; padding-bottom:5mm;">
            Page <span class="pageNumber"></span> of <span class="totalPages"></span>
        </div>
    '''

    print(f"\n📄 Generating PDF from: {html_file}")
    print(f"   Quality: {quality}")
    print(f"   Wait time: {wait_time}s (for MathJax rendering)")

    async with HTMLToPDFConverter() as converter:
        success = await converter.convert_file(
            str(html_file),
            str(output_pdf),
            options,
        )

        if not success:
            print("❌ PDF generation failed")
            return False

        print(f"✅ PDF generated: {output_pdf}")
        print(f"   Size: {output_pdf.stat().st_size / (1024 * 1024):.2f} MB")

        metadata = load_citation_metadata()
        print(f"   Version: v{metadata['version']}-{metadata['codename']}")

        final_pdf = copy_pdf_to_docs(output_pdf, docs_dir)
        process_pdf_with_metadata(final_pdf)
        copy_pdf_to_root(final_pdf, base_dir)

        print("\n✅ Complete! PDF available at:")
        print(f"   {final_pdf}")

        return True


def main():
    parser = argparse.ArgumentParser(
        description='Generate PDF from TEP-LENS site',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/generate_site_pdf.py
  python scripts/generate_site_pdf.py --quality maximum --wait-time 10
  python scripts/generate_site_pdf.py --skip-build
        """,
    )
    parser.add_argument(
        '--quality',
        choices=['maximum', 'high', 'print', 'web'],
        default='high',
        help='PDF quality preset (default: high)',
    )
    parser.add_argument(
        '--wait-time',
        type=float,
        default=5.0,
        help='Seconds to wait for MathJax rendering',
    )
    parser.add_argument(
        '--skip-build',
        action='store_true',
        help='Skip site build and use existing dist/ directory',
    )

    args = parser.parse_args()

    try:
        success = asyncio.run(generate_pdf(
            quality=args.quality,
            wait_time=args.wait_time,
            skip_build=args.skip_build,
        ))
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n⚠️  Cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()

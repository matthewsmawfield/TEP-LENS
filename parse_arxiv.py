import urllib.request
import tarfile
import io
import re
import sys

def parse_arxiv(arxiv_id, out_file):
    url = f'https://arxiv.org/e-print/{arxiv_id}'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        response = urllib.request.urlopen(req)
        tar = tarfile.open(fileobj=io.BytesIO(response.read()), mode='r:gz')
        with open(out_file, 'w') as out_f:
            for member in tar.getmembers():
                if member.name.endswith('.tex'):
                    out_f.write(f'Parsing {member.name}...\n')
                    f = tar.extractfile(member)
                    content = f.read().decode('utf-8', errors='ignore')
                    lines = content.split('\n')
                    for i, line in enumerate(lines):
                        if 'delay' in line.lower() and ('day' in line.lower() or 'dt' in line.lower() or r'\Delta t' in line):
                            out_f.write(f'L{i}: {line.strip()[:150]}\n')
                        if 'magnification' in line.lower() or r'\mu' in line:
                            out_f.write(f'M{i}: {line.strip()[:150]}\n')
        print(f"Successfully parsed {arxiv_id} to {out_file}")
    except Exception as e:
        print(f'Error parsing {arxiv_id}: {e}')

if __name__ == "__main__":
    parse_arxiv('2509.12301', 'encore_data.txt')
    parse_arxiv('2510.11719', 'h0pe_data.txt')
    parse_arxiv('2510.21694', 'winny_data.txt')

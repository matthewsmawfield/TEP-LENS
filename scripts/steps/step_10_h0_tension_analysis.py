import urllib.request
import tarfile
import io

urls = [
    'https://arxiv.org/e-print/2403.04873', # H0pe H0 measurement
    'https://arxiv.org/e-print/2403.08865', # Testing Lens Models
]

for url in urls:
    print(f"\n--- {url} ---")
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        response = urllib.request.urlopen(req)
        tar = tarfile.open(fileobj=io.BytesIO(response.read()), mode='r:gz')
        for member in tar.getmembers():
            if member.name.endswith('.tex'):
                f = tar.extractfile(member)
                content = f.read().decode('utf-8')
                for line in content.split('\n'):
                    if 'H_0' in line and ('km' in line or 'Mpc' in line):
                        pass
                        # print(line.strip())
                    if 'absolute' in line.lower() and 'delay' in line.lower():
                        print("ABSOLUTE DELAY:", line.strip()[:200])
                    if 'excess' in line.lower() and 'delay' in line.lower():
                        print("EXCESS DELAY:", line.strip()[:200])
    except Exception as e:
        print(f"Error: {e}")

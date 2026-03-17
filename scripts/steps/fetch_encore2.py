import urllib.request
import tarfile
import io

url = 'https://arxiv.org/e-print/2501.10286'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    response = urllib.request.urlopen(req)
    tar = tarfile.open(fileobj=io.BytesIO(response.read()), mode='r:gz')
    for member in tar.getmembers():
        if member.name.endswith('.tex'):
            f = tar.extractfile(member)
            content = f.read().decode('utf-8')
            for line in content.split('\n'):
                if 'Encore' in line and ('mu' in line or '\\mu' in line or 'magnification' in line.lower()):
                    print(line[:300])
except Exception as e:
    pass

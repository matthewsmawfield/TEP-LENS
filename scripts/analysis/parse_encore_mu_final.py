import urllib.request
import tarfile
import io

url = 'https://arxiv.org/e-print/2501.10286'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    response = urllib.request.urlopen(req)
    tar = tarfile.open(fileobj=io.BytesIO(response.read()), mode='r:gz')
    for member in tar.getmembers():
        if member.name.endswith('Results_3.tex'):
            f = tar.extractfile(member)
            content = f.read().decode('utf-8')
            for line in content.split('\n'):
                if '&' in line and 'Encore' in line:
                    print(line.strip())
except Exception as e:
    pass

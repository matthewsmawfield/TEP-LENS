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
            for i, line in enumerate(content.split('\n')):
                if 'Encore' in line and ('1a' in line.lower() or '1b' in line.lower() or '1c' in line.lower()):
                    if '&' in line and '\\' in line:
                        print(line)
except Exception as e:
    print(e)

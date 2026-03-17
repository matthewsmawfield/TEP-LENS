import urllib.request
import xml.etree.ElementTree as ET

urls = [
    "http://export.arxiv.org/api/query?search_query=id:2510.21694",
    "http://export.arxiv.org/api/query?search_query=id:2509.12301",
    "http://export.arxiv.org/api/query?search_query=id:2506.03023"
]

for url in urls:
    try:
        response = urllib.request.urlopen(url).read().decode("utf-8")
        root = ET.fromstring(response)
        ns = {"atom": "http://www.w3.org/2005/Atom"}
        entry = root.find("atom:entry", ns)
        if entry is not None:
            title = entry.find("atom:title", ns).text.strip()
            summary = entry.find("atom:summary", ns).text.strip()
            print(f"TITLE: {title}\nSUMMARY: {summary[:500]}...\n" + "-"*40)
    except Exception as e:
        print(f"Error fetching {url}: {e}")

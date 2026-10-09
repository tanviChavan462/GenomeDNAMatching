import re
with open("dashboard.html", "r", encoding="utf-8") as f:
    html = f.read()
matches = re.findall(r'alt="([^"]{0,80})"', html)
for m in matches:
    print(repr(m))

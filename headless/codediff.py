"""List every CODE line that differs between two sets of sources
(comment lines ignored).   python3 codediff.py DIR_A DIR_B"""
import sys, os, difflib
a, b = sys.argv[1:3]
def code(path):
    out = []
    for ln in open(path, encoding="latin-1").read().replace("\x1a", "").split("\n"):
        if ln[:1] in ("C", "c", "*") or not ln.strip():
            continue
        out.append(ln.rstrip())
    return out
tot = 0
for f in sorted(os.listdir(b)):
    pa, pb = os.path.join(a, f), os.path.join(b, f)
    if not os.path.exists(pa):
        print(f"== {f}: new file, {len(code(pb))} code lines"); continue
    d = [l for l in difflib.unified_diff(code(pa), code(pb), lineterm="", n=0) if l[:1] in "+-" and l[:3] not in ("+++", "---")]
    if d:
        print(f"== {f}: {len(d)} changed code lines"); print("\n".join(d)); tot += len(d)
print("total changed:", tot)

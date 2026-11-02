"""Tests for server.py public mode and the rate limiter.  Run: python3 test_eurosem_server.py
Uses a stand-in model, so no Fortran is needed."""
import os, sys, json, threading, urllib.request, urllib.error
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import server
from eurosem_form import DEFAULT_SPEC
fails = []
def check(n, c, d=""):
    print(("PASS  " if c else "FAIL  ") + n + ("" if c else "   " + str(d)[:300]))
    if not c: fails.append(n)

rl = server.RateLimiter(2, 100, 5)
check("limiter: two allowed then refused", rl.allow("a", 0)[0] and rl.allow("a", 1)[0] and not rl.allow("a", 2)[0])
check("limiter: another visitor unaffected", rl.allow("b", 3)[0])
ok, wait = rl.allow("a", 50)
check("limiter: wait time reported", (not ok) and 49 <= wait <= 52, wait)
check("limiter: window expires", rl.allow("a", 101)[0])
g = server.RateLimiter(100, 100, 2)
check("limiter: global cap", g.allow("x", 0)[0] and g.allow("y", 1)[0] and not g.allow("z", 2)[0] and g.allow("z", 61)[0])

def call(port, path, body=None, headers=None):
    req = urllib.request.Request("http://127.0.0.1:%d%s" % (port, path), data=None if body is None else json.dumps(body).encode(),
                                 headers=dict({"Content-Type": "application/json"}, **(headers or {})))
    try:
        with urllib.request.urlopen(req, timeout=60) as r: return r.status, r.read()
    except urllib.error.HTTPError as e: return e.code, e.read()

stub = os.environ.get("STUB_EXE", "/tmp/stubexe")
if not (os.path.isfile(stub) and os.access(stub, os.X_OK)):
    import tempfile
    stub = os.path.join(tempfile.mkdtemp(), "stub")
    ref = os.path.join(HERE, "reference")
    with open(stub, "w") as f:
        f.write("#!/bin/bash\nmapfile -t L < eurosem.ctl\ncp '%s/ref.dyn' \"${L[3]}\"; cp '%s/ref.sta' \"${L[5]}\"; touch \"${L[4]}\"; echo OK > eurosem.status\n" % (ref, ref))
    os.chmod(stub, 0o755)
try:
    srv = server.make_server(stub, "127.0.0.1", 0, public=True, rate_runs=2, rate_window=600, global_runs=50)
    port = srv.server_address[1]
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    check("healthz answers", call(port, "/healthz") == (200, b"ok"))
    check("config says public", json.loads(call(port, "/api/config")[1]) == {"public": True})
    check("sensitivity switched off", call(port, "/api/sensitivity", {"spec": DEFAULT_SPEC})[0] == 404 and call(port, "/api/sensitivity/plan", {"spec": DEFAULT_SPEC})[0] == 404)
    codes = [call(port, "/api/run", DEFAULT_SPEC, {"X-Forwarded-For": "1.2.3.4"})[0] for _ in range(3)]
    check("third run from one visitor is refused with 429", codes == [200, 200, 429], codes)
    check("a different visitor still runs", call(port, "/api/run", DEFAULT_SPEC, {"X-Forwarded-For": "5.6.7.8"})[0] in (200,))
    srv.shutdown()
    loc = server.make_server(stub, "127.0.0.1", 0, public=False)
    p2 = loc.server_address[1]; threading.Thread(target=loc.serve_forever, daemon=True).start()
    check("local mode keeps sensitivity route", call(p2, "/api/sensitivity", {"spec": {}})[0] != 404)
    check("local mode config", json.loads(call(p2, "/api/config")[1]) == {"public": False})
    loc.shutdown()
except Exception as ex:
    check("server tests ran", False, repr(ex))
print("\n%d failure(s)" % len(fails)); sys.exit(1 if fails else 0)

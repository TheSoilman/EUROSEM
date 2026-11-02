"""
server.py  --  local web service for EUROSEM (planes only).

    bash build_exe.sh             (once: compiles the model)
    python3 server.py             then open http://127.0.0.1:8080

Routes
  GET  /              the form and results page (app.html)
  GET  /api/defaults  the default form values (the two-plane reference case)
  POST /api/sensitivity[/plan]  one-at-a-time sensitivity (eurosem_sens.py)
  GET  /api/guide     soil-texture guide values from the user guide (eurosem_guide.py)
  POST /api/run       body: form values as JSON; returns the run result as JSON

The service uses only the Python standard library.  It listens on the
local machine only (127.0.0.1) unless you pass --host.  Each request runs
the model in its own temporary folder (see eurosem_run.py).
Options:  --exe PATH  --port N  --host ADDR  --timeout SECONDS
          --public   public deployment (see below); also set by EUROSEM_PUBLIC=1

Public mode (for a hosted copy, see ../DEPLOY.md): the sensitivity routes are switched
off, the page hides the sensitivity section, each visitor is limited to
--rate-runs model runs per --rate-window seconds (and all visitors together to
--global-runs per minute), the service listens on 0.0.0.0 and on the port named by
the PORT environment variable if --port is not given, and GET /healthz answers "ok".
Behind a proxy the visitor is taken from the first X-Forwarded-For address.
"""
import argparse, collections, json, os, sys, threading, time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import eurosem_io as E
from eurosem_form import build_run, DEFAULT_SPEC, FormError
from eurosem_run import run_case
from eurosem_export import make_results_csv, make_settings_csv
from eurosem_guide import guide_json
import eurosem_sens as SENS

MAX_BODY = 200_000
RUN_LOCK = threading.BoundedSemaphore(2)      # at most two model runs at once


class RateLimiter:
    """Sliding-window limit per visitor, plus one limit for everybody together."""
    def __init__(self, per_visitor, window, global_per_min):
        self.per_visitor, self.window, self.global_per_min = per_visitor, window, global_per_min
        self.hits = collections.defaultdict(collections.deque)
        self.all = collections.deque()
        self.lock = threading.Lock()

    def allow(self, who, now=None):
        """Return (True, 0) and record the run, or (False, seconds to wait)."""
        now = time.monotonic() if now is None else now
        with self.lock:
            q = self.hits[who]
            while q and q[0] <= now - self.window: q.popleft()
            while self.all and self.all[0] <= now - 60: self.all.popleft()
            if len(q) >= self.per_visitor:
                return False, max(1, int(q[0] + self.window - now) + 1)
            if len(self.all) >= self.global_per_min:
                return False, max(1, int(self.all[0] + 60 - now) + 1)
            q.append(now); self.all.append(now)
            if len(self.hits) > 5000:                      # forget visitors with nothing recent
                for k in [k for k, v in self.hits.items() if not v or v[-1] <= now - self.window]:
                    self.hits.pop(k, None)
            return True, 0


class Handler(BaseHTTPRequestHandler):
    exe = None
    timeout = 60.0
    public = False
    limiter = None

    def _who(self):
        fwd = self.headers.get("X-Forwarded-For", "")
        return fwd.split(",")[0].strip() or self.client_address[0] if self.public else self.client_address[0]

    def _send(self, code, body, ctype="application/json"):
        data = body if isinstance(body, bytes) else json.dumps(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype + "; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path in ("/", "/index.html"):
            with open(os.path.join(HERE, "app.html"), "rb") as f:
                self._send(200, f.read(), "text/html")
        elif self.path == "/api/defaults":
            self._send(200, DEFAULT_SPEC)
        elif self.path == "/api/config":
            self._send(200, {"public": self.public})
        elif self.path == "/healthz":
            self._send(200, b"ok", "text/plain")
        elif self.path == "/api/guide":
            self._send(200, guide_json())
        else:
            self._send(404, {"error": "not found"})

    def do_POST(self):
        if self.path not in ("/api/run", "/api/sensitivity", "/api/sensitivity/plan"):
            return self._send(404, {"error": "not found"})
        if self.public and self.path != "/api/run":
            return self._send(404, {"error": "not available on this site"})
        try:
            n = int(self.headers.get("Content-Length", "0"))
            if n <= 0 or n > MAX_BODY:
                return self._send(413, {"ok": False, "stage": "form", "error": "request too large or empty"})
            spec = json.loads(self.rfile.read(n))
            if not isinstance(spec, dict):
                raise ValueError("expected a JSON object")
        except (ValueError, json.JSONDecodeError) as exc:
            return self._send(400, {"ok": False, "stage": "form", "error": f"could not read the request: {exc}"})
        if self.path != "/api/run":
            return self._sensitivity(spec)
        try:
            run, accept = build_run(spec)
        except FormError as exc:
            return self._send(200, {"ok": False, "stage": "form", "error": str(exc), "problems": exc.problems})
        if self.limiter is not None:
            ok, wait = self.limiter.allow(self._who())
            if not ok:
                return self._send(429, {"ok": False, "stage": "run",
                                        "error": "too many runs from this connection; please wait about %d s and try again" % wait})
        if not RUN_LOCK.acquire(timeout=20 if self.public else 5):
            return self._send(503, {"ok": False, "stage": "run", "error": "the model is busy, try again in a moment"})
        try:
            out = run_case(run, self.exe, timeout=self.timeout, accept=accept)
        finally:
            RUN_LOCK.release()
        if out.get("ok"):
            try:
                out["csv"] = make_results_csv(spec, out["results"])
                out["settings_csv"] = make_settings_csv(spec)
            except Exception as exc:      # the run itself succeeded; only the export failed
                out["csv_error"] = str(exc)
        self._send(200, out)

    def _sensitivity(self, body):
        """body: {"spec": form values, "params": [keys], "overrides": {key: {low, high}}}"""
        form = body.get("spec")
        if not isinstance(form, dict):
            return self._send(400, {"ok": False, "error": "no form values supplied"})
        try:
            build_run(form)                  # the base case must be a valid form first
        except FormError as exc:
            return self._send(200, {"ok": False, "stage": "form", "error": str(exc), "problems": exc.problems})
        if self.path.endswith("/plan"):
            return self._send(200, {"ok": True, "plan": SENS.plan(form)})
        keys, ov = body.get("params") or [], body.get("overrides") or {}
        known = {p[0] for p in SENS.PARAMS}
        if not keys or not all(k in known for k in keys):
            return self._send(200, {"ok": False, "error": "choose at least one parameter"})
        for k, v in ov.items():
            try:
                for e in ("low", "high"):
                    if v.get(e) not in (None, ""):
                        float(v[e])
            except (TypeError, ValueError, AttributeError):
                return self._send(200, {"ok": False, "error": "the low and high values for %s must be numbers" % SENS.LABEL.get(k, (k,))[0]})
        exe, timeout = self.exe, self.timeout

        def run_one(sp):
            try:
                run, accept = build_run(sp)
            except FormError as exc:
                return False, "; ".join(exc.problems)
            o = run_case(run, exe, timeout=timeout, accept=accept)
            return (True, o["results"]) if o.get("ok") else (False, o.get("error", "model failed"))

        if not RUN_LOCK.acquire(timeout=5):
            return self._send(503, {"ok": False, "error": "the model is busy, try again in a moment"})
        try:
            out = SENS.run_sensitivity(form, keys, ov, run_one)
        except Exception as exc:
            return self._send(200, {"ok": False, "error": "sensitivity failed: %s" % exc})
        finally:
            RUN_LOCK.release()
        if out.get("ok"):
            try:
                out["csv"] = SENS.sensitivity_csv(form, out)
            except Exception as exc:
                out["csv_error"] = str(exc)
        self._send(200, out)

    def log_message(self, fmt, *a):
        sys.stderr.write("[%s] %s\n" % (self.log_date_time_string(), fmt % a))


def make_server(exe, host, port, timeout=60.0, public=False, rate_runs=20, rate_window=600, global_runs=40):
    """Build the HTTP server (used by main and the tests).  Settings are held on a
    per-server Handler subclass, so several servers can coexist in one process."""
    attrs = dict(exe=exe, timeout=timeout, public=public,
                 limiter=RateLimiter(rate_runs, rate_window, global_runs) if public else None)
    return ThreadingHTTPServer((host, port), type("ConfiguredHandler", (Handler,), attrs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=os.environ.get("EUROSEM_EXE", os.path.join(HERE, "eurosem")))
    ap.add_argument("--port", type=int, default=None)
    ap.add_argument("--host", default=None)
    ap.add_argument("--timeout", type=float, default=60.0)
    ap.add_argument("--public", action="store_true", default=os.environ.get("EUROSEM_PUBLIC") == "1")
    ap.add_argument("--rate-runs", type=int, default=20, help="public mode: runs per visitor per window")
    ap.add_argument("--rate-window", type=int, default=600, help="public mode: window in seconds")
    ap.add_argument("--global-runs", type=int, default=40, help="public mode: runs per minute, all visitors")
    a = ap.parse_args()
    if a.port is None:
        a.port = int(os.environ.get("PORT", "8080")) if a.public else 8080
    if a.host is None:
        a.host = "0.0.0.0" if a.public else "127.0.0.1"
    if not os.path.isfile(a.exe) or not os.access(a.exe, os.X_OK):
        sys.exit(f"Model program not found: {a.exe}\nBuild it first:  bash build_exe.sh")
    srv = make_server(a.exe, a.host, a.port, a.timeout, a.public, a.rate_runs, a.rate_window, a.global_runs)
    print(f"EUROSEM service on http://{a.host}:{a.port}  (model: {a.exe}; {'PUBLIC mode' if a.public else 'local mode'})  Ctrl-C to stop")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")

if __name__ == "__main__":
    main()

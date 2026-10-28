#!/usr/bin/env python3
"""Run one job file on the server and collect its small outputs into results/<name>/.

Usage:  python runner/run_job.py jobs/0001-test.toml

Settings come from environment variables (set in runner/config.sh):
  DEV_DIR       persistent clone of the code repository on the server
  MAX_CORES     upper limit on the cores a job may ask for
  GAIAMOCK_DIR  location of gaiamock (passed on to the job)

The job's command runs inside DEV_DIR, at the code version named in the job,
(nice 0 = normal priority by default; a job can set nice = 1..19 to yield to other users) and with
every parallel library limited to the number of cores asked for.
Only small result files are copied back; large outputs stay on the server.

While a job runs, its results/<name>/status.txt is rewritten and pushed every PROGRESS_MINUTES
(default 60): time so far, the last progress line of the output ("k/N" counts, as the project's
scripts print them), the share done and the estimated time left for the current step. Pushing is
best effort; a failed progress push never stops the job.
"""
from __future__ import annotations

import datetime as dt
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import tomllib
from pathlib import Path

RUNS = Path(__file__).resolve().parents[1]
KEEP_SUFFIXES = {".txt", ".csv", ".tsv", ".png", ".pdf", ".svg", ".json", ".md", ".dat", ".log"}
MAX_FILE_MB = 20      # larger files are skipped (they stay on the server)
MAX_TOTAL_MB = 100    # per job
LOG_TAIL_LINES = 3000
PROGRESS_MINUTES = float(os.environ.get("PROGRESS_MINUTES", "60"))
COUNT = re.compile(r"(\d[\d,]*)\s*/\s*(\d[\d,]*)")


def _hm(sec: float) -> str:
    sec = max(0.0, sec)
    return f"{int(sec // 3600)} h {int(sec % 3600 // 60):02d} min"


class Progress:
    """Follows "k/N" counts in the job output and estimates the time left of the current step."""

    def __init__(self):
        self.last = ""          # last line with a count
        self.lastline = ""      # last line of any kind
        self.total = None       # N of the current step
        self.first = None       # (time, k) when this N was first seen
        self.k = 0
        self.lock = threading.Lock()

    def feed(self, line: str) -> None:
        line = line.strip()
        if not line:
            return
        with self.lock:
            self.lastline = line[:200]
            m = COUNT.search(line)
            if not m:
                return
            k, n = (int(x.replace(",", "")) for x in m.groups())
            if n <= 0 or k > n:
                return
            if n != self.total or k < self.k:
                self.total, self.first = n, (dt.datetime.now(), k)
            self.k, self.last = k, line[:200]

    def report(self) -> list[str]:
        with self.lock:
            if self.total is None:
                return [f"last output: {self.lastline or '(none yet)'}"]
            out = [f"progress:   {self.k:,}/{self.total:,} ({100 * self.k / self.total:.0f}% of the current step)"]
            t1, k1 = self.first
            el = (dt.datetime.now() - t1).total_seconds()
            if self.k > k1 and el > 0:
                left = el / (self.k - k1) * (self.total - self.k)
                end = dt.datetime.now(dt.timezone.utc) + dt.timedelta(seconds=left)
                out.append(f"time left:  about {_hm(left)} for this step (ends about {end:%Y-%m-%d %H:%M} UTC; "
                           "a job with several steps takes longer)")
            out.append(f"last line:  {self.last}")
            return out


def push_status(path: Path, text: str, msg: str) -> None:
    """Write status.txt and push it, best effort (the job continues whatever happens)."""
    path.write_text(text)
    try:
        def g(*a):
            subprocess.run(["git", "-C", str(RUNS), *a], check=True, capture_output=True, timeout=120)
        g("add", str(path))
        g("commit", "-q", "-m", msg)
        for _ in range(3):
            try:
                g("pull", "-q", "--rebase")
                g("push", "-q")
                return
            except subprocess.SubprocessError:
                continue
    except (subprocess.SubprocessError, OSError):
        pass


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def git(dev: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(dev), *args], check=True, capture_output=True,
                          text=True, timeout=600).stdout.strip()


def load_job(path: Path) -> dict:
    job = tomllib.loads(path.read_text())
    for key in ("name", "command"):
        if key not in job:
            sys.exit(f"{path}: missing '{key}'")
    if not re.fullmatch(r"[A-Za-z0-9._-]+", job["name"]):
        sys.exit(f"{path}: name may only contain letters, digits, '.', '_' and '-'")
    job.setdefault("code_ref", "main")
    job.setdefault("cores", 8)
    job.setdefault("collect", [])
    job.setdefault("max_hours", 100)
    job.setdefault("nice", 0)
    if not 0 <= int(job["nice"]) <= 19:
        sys.exit(f"{path}: nice must be between 0 and 19")
    max_cores = int(os.environ.get("MAX_CORES", "96"))
    if not 1 <= int(job["cores"]) <= max_cores:
        sys.exit(f"{path}: cores must be between 1 and {max_cores}")
    return job


def checkout(dev: Path, ref: str) -> str:
    git(dev, "fetch", "--quiet", "--tags", "origin")
    try:
        git(dev, "checkout", "--quiet", "--detach", f"origin/{ref}")
    except subprocess.CalledProcessError:
        git(dev, "checkout", "--quiet", "--detach", ref)       # tag or commit
    return git(dev, "rev-parse", "HEAD")


def job_env(cores: int) -> dict:
    env = dict(os.environ)
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[var] = "1"                         # one thread per worker; joblib supplies the workers
    env["LOKY_MAX_CPU_COUNT"] = str(cores)     # joblib's n_jobs=-1 then means `cores`
    env["JOB_CORES"] = str(cores)
    return env


def collect(dev: Path, patterns: list[str], out: Path) -> list[str]:
    """Copy small result files; return a list of what was skipped and why."""
    skipped, total = [], 0.0
    for pat in patterns:
        matches = sorted(dev.glob(pat)) or [dev / pat]
        for src in matches:
            if not src.exists():
                skipped.append(f"{pat}: not found")
                continue
            files = [src] if src.is_file() else sorted(p for p in src.rglob("*") if p.is_file())
            for f in files:
                rel = f.relative_to(dev)
                mb = f.stat().st_size / 1e6
                if f.suffix.lower() not in KEEP_SUFFIXES:
                    skipped.append(f"{rel}: type {f.suffix or '(none)'} not copied")
                elif mb > MAX_FILE_MB:
                    skipped.append(f"{rel}: {mb:.0f} MB, too large")
                elif total + mb > MAX_TOTAL_MB:
                    skipped.append(f"{rel}: job total over {MAX_TOTAL_MB} MB")
                else:
                    dest = out / "files" / rel
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(f, dest)
                    total += mb
    return skipped


def main() -> int:
    job_file = Path(sys.argv[1])
    job = load_job(job_file)
    dev = Path(os.environ["DEV_DIR"])
    out = RUNS / "results" / job["name"]
    out.mkdir(parents=True, exist_ok=True)
    cores = int(job["cores"])
    start, t0 = now(), dt.datetime.now()

    commit = checkout(dev, job["code_ref"])
    command = job["command"].replace("{cores}", str(cores))
    print(f"job {job['name']}: code {commit[:10]} ({job['code_ref']}), {cores} cores on "
          f"{socket.gethostname()}\n$ {command}", flush=True)

    log_path = out / "log.txt"
    lines: list[str] = []
    timed_out = False
    prog = Progress()
    stop = threading.Event()

    def reporter():
        while not stop.wait(PROGRESS_MINUTES * 60):
            el = (dt.datetime.now() - t0).total_seconds()
            text = "\n".join([f"status:     running on {socket.gethostname()} since {start}",
                              f"running:    {_hm(el)} so far (updated {now()})", *prog.report()]) + "\n"
            push_status(out / "status.txt", text, f"Progress: {job['name']} ({_hm(el)})")

    th = threading.Thread(target=reporter, daemon=True)
    th.start()
    with subprocess.Popen(["nice", "-n", str(int(job["nice"])), "bash", "-c", command], cwd=dev, env=job_env(cores),
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
        deadline = t0 + dt.timedelta(hours=float(job["max_hours"]))
        for line in proc.stdout:
            print(line, end="", flush=True)    # live in the GitHub Actions log
            lines.append(line)
            prog.feed(line)
            if len(lines) > 2 * LOG_TAIL_LINES:
                lines = lines[-LOG_TAIL_LINES:]
            if dt.datetime.now() > deadline:
                timed_out = True
                proc.terminate()
                break
        code = proc.wait()
    stop.set()
    th.join(timeout=600)      # let a progress push in flight finish before the final status is written
    log_path.write_text("".join(lines[-LOG_TAIL_LINES:]))

    skipped = collect(dev, list(job["collect"]), out)
    hours = (dt.datetime.now() - t0).total_seconds() / 3600
    status = "timed out" if timed_out else ("finished" if code == 0 else f"failed (exit {code})")
    meta = [
        f"job:        {job_file.name}",
        f"status:     {status}",
        f"code:       {commit} ({job['code_ref']})",
        f"command:    {command}",
        f"host:       {socket.gethostname()}, {cores} cores",
        f"started:    {start}",
        f"ended:      {now()}  ({hours:.2f} h)",
        f"log:        last {LOG_TAIL_LINES} lines in log.txt",
    ]
    if skipped:
        meta += ["", "Not copied (still on the server):"] + [f"  {s}" for s in skipped]
    (out / "status.txt").write_text("\n".join(meta) + "\n")
    print("\n".join(meta))
    backup = os.environ.get("BACKUP_DIR")
    if backup:                      # kept on the server even if pushing to GitHub fails
        shutil.copytree(out, Path(backup) / job["name"], dirs_exist_ok=True)
    return 0 if code == 0 and not timed_out else 1


if __name__ == "__main__":
    sys.exit(main())

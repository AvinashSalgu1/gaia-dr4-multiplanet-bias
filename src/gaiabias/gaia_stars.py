"""Real Gaia DR3 stars for step 3.

Download (once; needs internet and `pip install astroquery`):
Gaia DR3 sources (gaiadr3.gaia_source) with

    parallax > 4 mas (closer than 250 pc)   and   parallax / error > 10

inside a RANDOM SET OF SMALL SKY CIRCLES (radius 3 deg, centres on an even
Fibonacci grid of 768 points, so circles never overlap), covering a chosen
fraction of the sky (default 25%). Why not the whole sky: any all-sky query on
parallax has to visit millions of scattered rows of the 1.8-billion-row table
and took 15-25+ minutes per piece on both archives. A circle query uses the
archives' sky-position index and reads only that patch, so each piece takes
seconds to a minute, and progress is visible right away. The stars of nearby
space look the same in every direction, so a random sky fraction gives a fair
sample of stars (all planet rates in step 3 are per star); step 3 only uses a
random subset of ~300k stars anyway. Circles are fetched in a fixed random
order, so stopping early still leaves a random sky sample.

The table is served by two archives (TAP services); the code can use either:
  * "vizier" - CDS VizieR, catalogue I/355/gaiadr3 (identical DR3 values;
    VizieR column names are renamed to the Gaia ones in the query);
  * "esa"    - the ESA Gaia archive (gea.esac.esa.int).
With source="auto" (default) VizieR is tried first; a piece that does not
finish within STALL_S on one archive is requested from the other. Optional ESA
login (faster queue): user name and password on two lines in ~/.gaia_credentials.

The FLAME stellar masses (Creevey et al. 2023) are fetched afterwards from the
ESA archive ONLY for the stars actually used (uploading their source_ids), and
cached; if that is slow or fails, the mass-magnitude relation is used instead.

Selection, following Lammers & Winn (2026, AJ 171:18) Section 4.1:
  * absolute magnitude M_G = G + 5 log10(parallax / 100 mas), no extinction;
  * main sequence: 0.4 < M_G - 2.9 (BP - RP) < 4.5            (their Eq. 38);
  * distance below their r_max(M_G)                          (their Eqs. 36-37);
plus, for this project:
  * distance below 250 pc (they find > 90% of DR4 detections within 250 pc);
  * RUWE <= 1.4 by default (higher RUWE usually means an unresolved binary;
    every simulated star here is single). Stars above are kept in the file
    with a flag, so the cut can be changed.
Masses: the FLAME mass where Gaia provides one, otherwise from M_G with the
main-sequence mass-magnitude relation (as Lammers & Winn); stars fainter
than M_G = 12 (below ~0.2 M_sun, beyond the relation) are left out.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from .population import _MS_MASS, _MS_MG, StarSample, r_max_lw

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "gaia_dr3_250pc"
MERGED = DATA_DIR / "stars.npz"
N_CONES = 768              # Fibonacci grid on the sky; nearest neighbours are >= 6.4 deg apart
CONE_R_DEG = 3.0           # so 3-deg circles never overlap; each covers 28.3 deg^2 (0.069% of the sky)
SKY_FRACTION = 0.25        # default sky fraction downloaded (max 768 circles = 52.6%)
CONE_SEED = 2017
POLL_S = 30
STALL_S = 15 * 60          # give a piece this long on one archive before trying the other
DOWNLOAD_JOBS = 6          # circles requested at the same time (archive jobs run in parallel)
FLAME_TIMEOUT_S = 10 * 60
D_MAX_PC = 250.0
MG_FAINT = 12.0
COLS = ["source_id", "ra", "dec", "parallax", "parallax_error", "pmra", "pmdec", "phot_g_mean_mag",
        "bp_rp", "ruwe"]
FLAME_CACHE = DATA_DIR / "flame_masses.npz"
FLAME_BATCH = 50_000
CREDENTIALS = Path.home() / ".gaia_credentials"

ESA_WHERE = """
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra0}, {dec0}, {r}))
  AND parallax >= {plx_min}
  AND parallax_over_error > 10
  AND phot_g_mean_mag IS NOT NULL
  AND bp_rp IS NOT NULL
"""
ESA_QUERY = ("SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag, bp_rp, ruwe"
             + ESA_WHERE)
# VizieR I/355/gaiadr3 holds the same DR3 values under VizieR names (RA_ICRS/DE_ICRS are the
# epoch-2016.0 positions, i.e. gaia_source ra/dec); the aliases give the Gaia names back.
VIZIER_WHERE = """
FROM "I/355/gaiadr3"
WHERE 1 = CONTAINS(POINT('ICRS', "RA_ICRS", "DE_ICRS"), CIRCLE('ICRS', {ra0}, {dec0}, {r}))
  AND "Plx" >= {plx_min}
  AND "Plx" > 10 * "e_Plx"
  AND "Gmag" IS NOT NULL
  AND "BP-RP" IS NOT NULL
"""
VIZIER_QUERY = ("""SELECT "Source" AS source_id, "RA_ICRS" AS ra, "DE_ICRS" AS dec, "Plx" AS parallax,
       "e_Plx" AS parallax_error, "pmRA" AS pmra, "pmDE" AS pmdec, "Gmag" AS phot_g_mean_mag,
       "BP-RP" AS bp_rp, "RUWE" AS ruwe""" + VIZIER_WHERE)
MIRRORS = {
    "vizier": dict(label="CDS VizieR (I/355/gaiadr3)",
                   url="https://tapvizier.cds.unistra.fr/TAPVizieR/tap",
                   query=VIZIER_QUERY, where=VIZIER_WHERE,
                   check='SELECT TOP 1 "Source" FROM "I/355/gaiadr3"'),
    "esa": dict(label="ESA Gaia archive (gaiadr3.gaia_source)",
                url=None,
                query=ESA_QUERY, where=ESA_WHERE,
                check="SELECT TOP 1 source_id FROM gaiadr3.gaia_source"),
}
# sizes at which TAP services commonly truncate results: getting exactly one of these is checked
_SUSPICIOUS_ROWS = {1000, 2000, 3000, 5000, 10_000, 20_000, 50_000, 100_000}


def cone_centres(n: int = N_CONES):
    """(ra, dec) in degrees of an even Fibonacci grid on the sphere."""
    i = np.arange(n)
    dec = np.degrees(np.arcsin(1 - (2 * i + 1) / n))
    ra = np.degrees((i * np.pi * (3 - np.sqrt(5))) % (2 * np.pi))
    return ra, dec


def cone_area_deg2(r_deg: float = CONE_R_DEG):
    return 2 * np.pi * (1 - np.cos(np.radians(r_deg))) * (180 / np.pi) ** 2


def cone_order(seed: int = CONE_SEED):
    """Fixed random order of the circles (so any first k of them are a random sky sample)."""
    return np.random.default_rng(seed).permutation(N_CONES)


def n_cones_for(sky_fraction: float):
    n = int(np.ceil(sky_fraction * 41252.96 / cone_area_deg2()))
    if n > N_CONES:
        raise SystemExit(f"sky fraction {sky_fraction} is above the maximum "
                         f"{N_CONES * cone_area_deg2() / 41252.96:.3f} for non-overlapping circles")
    return max(n, 1)


FLAME_QUERY = """
SELECT u.source_id, ap.mass_flame
FROM tap_upload.ids AS u
JOIN gaiadr3.astrophysical_parameters AS ap ON ap.source_id = u.source_id
WHERE ap.mass_flame IS NOT NULL
"""


def _piece_path(k):
    return DATA_DIR / f"cone_{k:03d}.npz"


def _gaia(login=True):
    try:
        from astroquery.gaia import Gaia
    except ImportError as exc:
        raise SystemExit("astroquery is needed to download the Gaia stars: "
                         "pip install astroquery   (then rerun)") from exc
    Gaia.ROW_LIMIT = -1
    if login and CREDENTIALS.exists():
        try:
            Gaia.login(credentials_file=str(CREDENTIALS))
            print(f"  logged in to the Gaia archive (credentials from {CREDENTIALS})", flush=True)
        except Exception as exc:
            print(f"  login failed ({type(exc).__name__}: {exc}); continuing anonymously", flush=True)
    return Gaia


def _service(name, fresh=False):
    """A TAP client for one mirror (astroquery TapPlus; the Gaia class for ESA).
    fresh=True gives a separate client (one per download thread)."""
    if name == "esa":
        if not fresh:
            return _gaia()
        from astroquery.gaia import GaiaClass
        try:
            g = GaiaClass(show_server_messages=False)
        except TypeError:
            g = GaiaClass()
        if CREDENTIALS.exists():
            try:
                g.login(credentials_file=str(CREDENTIALS))
            except Exception:
                pass
        return g
    try:
        from astroquery.utils.tap.core import TapPlus
    except ImportError as exc:
        raise SystemExit("astroquery is needed to download the Gaia stars: "
                         "pip install astroquery   (then rerun)") from exc
    return TapPlus(url=MIRRORS[name]["url"], verbose=False)


def _table_to_arrays(tab):
    names = {c.lower(): c for c in tab.colnames}
    arr = {}
    for c in COLS:
        if c not in names:
            raise KeyError(f"column {c} missing from the archive answer (got {list(tab.colnames)})")
        col = tab[names[c]]
        if c == "source_id":
            arr[c] = np.asarray(col, dtype=np.int64)
        else:
            col = col.astype(float)
            arr[c] = np.asarray(col.filled(np.nan) if hasattr(col, "filled") else col, float)
    return arr


def _abort(job):
    """Best effort: tell the archive to stop a job we gave up on (UWS PHASE=ABORT)."""
    try:
        import requests
        loc = getattr(job, "remoteLocation", None)
        if loc:
            requests.post(loc.rstrip("/") + "/phase", data={"PHASE": "ABORT"}, timeout=20)
    except Exception:
        pass


def _fmt(sec):
    sec = int(sec)
    return f"{sec // 3600}h{(sec % 3600) // 60:02d}m" if sec >= 3600 else f"{sec // 60}m{sec % 60:02d}s"


def _run_job(tap, q, label, t_start, done_before, rows_before, n_total, timeout_s=STALL_S):
    """Launch one asynchronous archive job and report its progress until the rows are in."""
    t0 = time.perf_counter()
    try:
        job = tap.launch_job_async(q, background=True)
    except TypeError:                                   # very old astroquery: no background option
        print(f"  {label}: running (no progress available with this astroquery version) ...", flush=True)
        return tap.launch_job_async(q).get_results()
    last = 0.0
    while True:
        try:
            phase = job.get_phase(update=True)
        except Exception:                              # status not available: just wait for the rows
            phase = "UNKNOWN"
        if phase in ("COMPLETED", "ERROR", "ABORTED", "UNKNOWN"):
            break
        now = time.perf_counter()
        if now - t0 >= 60 and now - last >= POLL_S:
            el = now - t_start
            eta = (f", about {_fmt(el / done_before * (n_total - done_before) - (now - t0))} left in total"
                   if done_before else "")
            print(f"  {label}: archive status {phase}, waiting {_fmt(now - t0)} (gives up at {_fmt(timeout_s)})"
                  f"  [{done_before}/{n_total} new pieces done, {rows_before:,} stars so far{eta}]", flush=True)
            last = now
        if now - t0 > timeout_s:
            _abort(job)
            raise TimeoutError(f"no answer after {_fmt(timeout_s)}")
        time.sleep(5)
    if phase in ("ERROR", "ABORTED"):
        raise RuntimeError(f"archive job {phase}")
    return job.get_results()


def _check(name):
    """Connect to a mirror and run a one-row query. Returns the client, or None if it does not answer."""
    t0 = time.perf_counter()
    try:
        tap = _service(name)
        tap.launch_job(MIRRORS[name]["check"]).get_results()
        print(f"  {MIRRORS[name]['label']}: reachable ({time.perf_counter() - t0:.1f} s)", flush=True)
        return tap
    except SystemExit:
        raise
    except Exception as exc:
        print(f"  {MIRRORS[name]['label']}: not answering ({type(exc).__name__}: {str(exc)[:200]})", flush=True)
        return None


def _count(tap, name, fmt):
    q = "SELECT COUNT(*) AS n" + MIRRORS[name]["where"].format(**fmt)
    tab = tap.launch_job_async(q).get_results()
    return int(tab[tab.colnames[0]][0])


def download(sky_fraction: float = SKY_FRACTION, source: str = "auto", stall_s: float = STALL_S,
             plx_min: float = 1000.0 / D_MAX_PC, jobs: int = DOWNLOAD_JOBS):
    """Download the first n circles (in the fixed random order) not yet on disk, then merge.

    source: "auto" (VizieR first, ESA as fallback per piece; whichever answers keeps being
    used first), "vizier" or "esa" (only that one). jobs: circles requested in parallel."""
    import os
    import threading
    from concurrent.futures import ThreadPoolExecutor, as_completed

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    order = ["vizier", "esa"] if source == "auto" else [source]
    if any(o not in MIRRORS for o in order):
        raise SystemExit(f"unknown source {source}; choose auto, vizier or esa")
    want = cone_order()[:n_cones_for(sky_fraction)]
    todo = [k for k in want if not _piece_path(k).exists()]
    n_have = len(want) - len(todo)
    rows0 = sum(len(np.load(_piece_path(k))["ra"]) for k in want if _piece_path(k).exists())
    print(f"  sky sample: {len(want)} circles of {CONE_R_DEG:g} deg = {100 * sky_fraction:.1f}% of the sky; "
          f"{n_have} already on disk ({rows0:,} stars), {len(todo)} to fetch, {jobs} at a time", flush=True)
    if todo:
        print("  checking which archives answer ...", flush=True)
        order = [name for name in order if _check(name) is not None]
        if not order:
            raise SystemExit("no Gaia archive answers. Check the internet connection, or try later "
                             "(https://www.cosmos.esa.int/web/gaia/news lists ESA archive downtime).")
        timeout = stall_s if len(order) > 1 else 3600
        ra_c, dec_c = cone_centres()
        lock = threading.Lock()
        local = threading.local()
        st = dict(order=list(order), n_new=0, rows=rows0, t_start=time.perf_counter())

        def client(name):
            if not hasattr(local, "taps"):
                local.taps = {}
            if name not in local.taps:
                local.taps[name] = _service(name, fresh=True)
            return local.taps[name]

        def fetch(k):
            fmt = dict(ra0=f"{ra_c[k]:.5f}", dec0=f"{dec_c[k]:.5f}", r=f"{CONE_R_DEG:g}", plx_min=f"{plx_min:.6f}")
            where = f"circle {k:03d} (RA {ra_c[k]:5.1f}, Dec {dec_c[k]:+5.1f})"
            with lock:
                pref = list(st["order"])
            attempts = [(name, timeout) for name in pref] + [(name, 4 * timeout) for name in pref]
            for i, (name, t_lim) in enumerate(attempts):
                label = f"{where} from {name}"
                try:
                    tp = time.perf_counter()
                    with lock:
                        snap = (st["n_new"], st["rows"])
                    tab = _run_job(client(name), MIRRORS[name]["query"].format(**fmt), label, st["t_start"],
                                   snap[0], snap[1], len(todo), timeout_s=t_lim)
                    arr = _table_to_arrays(tab)
                    n_rows = len(arr["ra"])
                    if n_rows in _SUSPICIOUS_ROWS and _count(client(name), name, fmt) != n_rows:
                        raise RuntimeError(f"exactly {n_rows:,} rows - result cut off by a row limit")
                    np.savez_compressed(_piece_path(k), source=np.array(name), **arr)
                    with lock:
                        st["n_new"] += 1
                        st["rows"] += n_rows
                        n, el = st["n_new"], time.perf_counter() - st["t_start"]
                        left = el / n * (len(todo) - n)
                        print(f"  piece {n_have + n}/{len(want)}: {label}: {n_rows:,} stars in "
                              f"{_fmt(time.perf_counter() - tp)}  ->  {st['rows']:,} stars, "
                              f"about {_fmt(left)} left", flush=True)
                        if name != st["order"][0]:     # this archive works better now: use it first
                            st["order"] = [name] + [o for o in st["order"] if o != name]
                    return True
                except Exception as exc:      # stalled, network hiccup or archive busy: try the next option
                    with lock:                         # this archive struggles: try the others first
                        if len(st["order"]) > 1 and st["order"][0] == name:
                            st["order"] = st["order"][1:] + [name]
                    nxt = attempts[i + 1][0] if i + 1 < len(attempts) else None
                    print(f"  {label} failed ({type(exc).__name__}: {str(exc)[:200]}); "
                          + (f"trying {nxt}" if nxt else "giving up on this circle"), flush=True)
                    time.sleep(10)
            return False

        pool = ThreadPoolExecutor(max_workers=max(1, jobs))
        try:
            results = [f.result() for f in as_completed([pool.submit(fetch, k) for k in todo])]
        except KeyboardInterrupt:
            print("\n  stopped - finished circles are kept; rerun the same command to continue", flush=True)
            os._exit(1)
        pool.shutdown()
        if not all(results):
            raise SystemExit(f"{results.count(False)} circles failed on every archive; rerun later to fetch "
                             "them (finished ones are kept), or use a smaller --sky-fraction.")
    merge(want)


def merge(cones=None):
    """Join the downloaded circles into one file (the circles in `cones`, default: all on disk)."""
    if cones is None:
        cones = [k for k in range(N_CONES) if _piece_path(k).exists()]
    parts = [dict(np.load(_piece_path(k))) for k in cones]
    out = {c: np.concatenate([p[c] for p in parts]) for c in COLS}
    frac = len(cones) * cone_area_deg2() / 41252.96
    np.savez_compressed(MERGED, sky_fraction=np.array(frac), n_cones=np.array(len(cones)), **out)
    print(f"merged {len(out['ra']):,} Gaia DR3 sources from {len(cones)} circles "
          f"({100 * frac:.1f}% of the sky) -> {MERGED}")


def have_stars() -> bool:
    return MERGED.exists()


def _flame_cache() -> dict:
    if FLAME_CACHE.exists():
        f = np.load(FLAME_CACHE)
        return dict(zip(f["source_id"].tolist(), f["mass_flame"].tolist())), set(f["asked"].tolist())
    return {}, set()


def fetch_flame(source_ids) -> dict:
    """FLAME masses for the given stars (only those Gaia has), fetched in batches and cached.

    Returns {source_id: mass}. Stars Gaia has no FLAME mass for are simply absent.
    If the archive cannot be reached, returns what is cached (the mass-magnitude
    relation is then used for the rest).
    """
    masses, asked = _flame_cache()
    todo = np.array(sorted(set(int(x) for x in source_ids) - asked), dtype=np.int64)
    if len(todo) == 0:
        return masses
    try:
        from astropy.table import Table
        Gaia = _gaia()
    except (ImportError, SystemExit) as exc:
        print(f"  FLAME masses not fetched ({exc}); using the mass-magnitude relation", flush=True)
        return masses
    tmp = DATA_DIR / "upload_ids.xml"
    n_b = int(np.ceil(len(todo) / FLAME_BATCH))
    print(f"  fetching Gaia FLAME masses for {len(todo):,} stars in {n_b} batches ...", flush=True)
    for b in range(n_b):
        ids = todo[b * FLAME_BATCH:(b + 1) * FLAME_BATCH]
        try:
            Table({"source_id": ids}).write(tmp, format="votable", overwrite=True)
            t0 = time.perf_counter()
            try:
                job = Gaia.launch_job_async(FLAME_QUERY, upload_resource=str(tmp), upload_table_name="ids",
                                            background=True)
                while job.get_phase(update=True) not in ("COMPLETED", "ERROR", "ABORTED"):
                    if time.perf_counter() - t0 > FLAME_TIMEOUT_S:
                        _abort(job)
                        raise TimeoutError(f"no answer after {_fmt(FLAME_TIMEOUT_S)}")
                    time.sleep(5)
                tab = job.get_results()
            except TypeError:                      # very old astroquery: no background option
                tab = Gaia.launch_job_async(FLAME_QUERY, upload_resource=str(tmp),
                                            upload_table_name="ids").get_results()
            sid = np.asarray(tab["source_id"], np.int64)
            mf = np.asarray(tab["mass_flame"].astype(float).filled(np.nan) if hasattr(tab["mass_flame"], "filled")
                            else tab["mass_flame"], float)
            masses.update({int(a): float(m) for a, m in zip(sid, mf) if np.isfinite(m)})
            asked.update(int(x) for x in ids)
            np.savez_compressed(FLAME_CACHE, source_id=np.array(list(masses), np.int64),
                                mass_flame=np.array(list(masses.values()), float),
                                asked=np.array(sorted(asked), np.int64))
            print(f"    batch {b + 1}/{n_b}: {len(sid):,} FLAME masses ({time.perf_counter() - t0:.0f} s)", flush=True)
        except Exception as exc:
            print(f"    batch {b + 1}/{n_b} failed ({type(exc).__name__}: {str(exc)[:200]}); the remaining "
                  "stars use the mass-magnitude relation (rerun later to retry)", flush=True)
            break                              # archive struggling: do not wait on every batch
    return masses


def mass_from_MG(M_G):
    """Main-sequence mass [M_sun] from absolute G (inverse of population.abs_G);
    linear extrapolation of the last segment down to M_G = 12."""
    M_G = np.asarray(M_G, float)
    order = np.argsort(_MS_MG)
    m = np.interp(M_G, _MS_MG[order], _MS_MASS[order])
    slope = (_MS_MASS[0] - _MS_MASS[1]) / (_MS_MG[0] - _MS_MG[1])
    faint = M_G > _MS_MG.max()
    m[faint] = _MS_MASS[0] + slope * (M_G[faint] - _MS_MG[0])
    return m


def select(raw: dict, ruwe_max: float | None = 1.4) -> tuple[dict, dict]:
    """Apply the selection. Returns (selected columns incl. mass, counts after each cut)."""
    plx, G, bprp = raw["parallax"], raw["phot_g_mean_mag"], raw["bp_rp"]
    d = 1000.0 / plx
    M_G = G + 5 * np.log10(plx / 100.0)
    sky = f", {100 * float(raw['sky_fraction']):.1f}% of the sky" if "sky_fraction" in raw else ""
    counts = {f"downloaded (closer than 250 pc, parallax/error > 10{sky})": len(plx)}
    keep = np.isfinite(M_G) & np.isfinite(bprp)
    x = M_G - 2.9 * bprp
    keep &= (x > 0.4) & (x < 4.5)
    counts["main sequence (L&W Eq. 38)"] = int(keep.sum())
    keep &= d < np.minimum(r_max_lw(M_G), D_MAX_PC)
    counts["within L&W r_max(M_G) and 250 pc"] = int(keep.sum())
    keep &= M_G <= MG_FAINT
    counts[f"M_G <= {MG_FAINT:g} (mass relation valid)"] = int(keep.sum())
    ruwe_ok = raw["ruwe"] <= (ruwe_max if ruwe_max is not None else np.inf)
    n_high = int((keep & ~ruwe_ok).sum())
    keep &= ruwe_ok
    counts[f"RUWE <= {ruwe_max}" if ruwe_max is not None else "no RUWE cut"] = int(keep.sum())
    flame = raw.get("mass_flame", np.full(len(plx), np.nan))
    mass = np.where(np.isfinite(flame), flame, mass_from_MG(M_G))
    sel = {k: v[keep] for k, v in raw.items() if np.ndim(v) == 1 and len(v) == len(plx)}
    sel.update(M_G=M_G[keep], mass=mass[keep], distance_pc=d[keep], mass_from_flame=np.isfinite(flame[keep]))
    counts["(removed by the RUWE cut)"] = n_high
    return sel, counts


def load_sample(n: int | None, seed: int, ruwe_max: float | None = 1.4, flame: bool = True):
    """A random subset of n selected stars as a StarSample, plus selection info.
    FLAME masses are fetched (and cached) for the chosen stars only."""
    raw = dict(np.load(MERGED))
    sel, counts = select(raw, ruwe_max)
    idx = np.arange(len(sel["ra"]))
    if n is not None and n < len(idx):
        idx = np.sort(np.random.default_rng([seed, 77]).choice(idx, n, replace=False))
    if flame and "mass_flame" not in raw:
        fm = fetch_flame(sel["source_id"][idx])
        f_sel = np.array([fm.get(int(x), np.nan) for x in sel["source_id"][idx]])
        has = np.isfinite(f_sel)
        sel["mass"][idx[has]] = f_sel[has]
        sel["mass_from_flame"][idx] = has
    stars = StarSample(m_star=sel["mass"][idx], G=sel["phot_g_mean_mag"][idx], parallax=sel["parallax"][idx],
                       ra=sel["ra"][idx], dec=sel["dec"][idx], pmra=np.nan_to_num(sel["pmra"][idx]),
                       pmdec=np.nan_to_num(sel["pmdec"][idx]))
    info = dict(counts=counts, n_selected=len(sel["ra"]), n_used=len(idx),
                frac_flame_mass=float(sel["mass_from_flame"][idx].mean()),
                source_id=sel["source_id"][idx], ruwe=sel["ruwe"][idx])
    return stars, info

"""Read Gaia DR4 epoch astrometry (VOTable, BINARY2 serialisation) without astropy.

The DR4 pre-release (June 2026) gives epoch astrometry for 12 sources as one
VOTable. Each row is one field-of-view transit with per-CCD arrays (up to 10
CCDs: SM + AF1-AF9 for 2D/1D windows). See the DR4 draft data model, Sec. 3.1.

Weights: the data model says that to reproduce AGIS, the per-CCD variance is
    sigma_tot^2 = centroid_pos_error_al^2 + agis_source_excess_noise^2
(centroid_pos_error_al already contains the IPD and calibration terms). SM
errors are "significantly underestimated", so SM is excluded (as AGIS does via
used_by_agis_al).
"""
from __future__ import annotations

import base64
import re
import struct

import numpy as np

_SCALAR = {"long": ("q", 8), "double": ("d", 8), "float": ("f", 4), "short": ("h", 2), "int": ("i", 4),
           "unsignedByte": ("B", 1), "boolean": ("c", 1)}


def read_votable_binary2(path):
    """Return {field name: list of values (arrays for variable-length fields)}."""
    txt = open(path, encoding="utf-8").read()
    fields = [(m.group(2), m.group(1) is not None)
              for m in re.finditer(r'<FIELD\s+(arraysize="\*"\s+)?datatype="(\w+)"', txt)]
    names = re.findall(r'<FIELD[^>]*\sname="([^"]+)"', txt)
    stream = re.search(r"<STREAM encoding='base64'>(.*?)</STREAM>", txt, re.S).group(1)
    raw = base64.b64decode("".join(stream.split()))
    nf = len(names)
    nb = (nf + 7) // 8
    out = {n: [] for n in names}
    pos = 0
    while pos < len(raw):
        nulls = raw[pos:pos + nb]
        pos += nb
        for k, (name, (dtype, var)) in enumerate(zip(names, fields)):
            fmt, size = _SCALAR[dtype]
            isnull = (nulls[k // 8] >> (7 - k % 8)) & 1
            if var:
                (n,) = struct.unpack(">i", raw[pos:pos + 4])
                pos += 4
                vals = raw[pos:pos + n * size]
                pos += n * size
                if dtype == "boolean":
                    v = np.array([c in (ord("T"), ord("t"), ord("1")) for c in vals])
                else:
                    v = np.array(struct.unpack(f">{n}{fmt}", vals), float if fmt in "df" else np.int64)
            else:
                vals = raw[pos:pos + size]
                pos += size
                if dtype == "boolean":
                    v = vals in (b"T", b"t", b"1")
                else:
                    (v,) = struct.unpack(f">{fmt}", vals)
            out[name].append(np.nan if (isnull and not var and fmt in "df") else v)
    return out


def ccd_table(path):
    """Flatten to one entry per CCD measurement: dict of numpy arrays."""
    t = read_votable_binary2(path)
    rows = {k: [] for k in ("source_id", "transit", "ccd", "t_ns", "psi_deg", "pf_al", "w", "err", "ipd_err",
                            "excess_src", "used", "g_mag", "ra0", "dec0")}
    for i in range(len(t["source_id"])):
        w = np.atleast_1d(t["centroid_pos_al"][i])
        n = len(w)
        for c in range(n):
            rows["source_id"].append(t["source_id"][i])
            rows["transit"].append(t["transit_id"][i])
            rows["ccd"].append(c)
            rows["t_ns"].append(np.atleast_1d(t["obs_time_tcb"][i])[c] + t["obs_time_bary_corr"][i])
            rows["psi_deg"].append(np.atleast_1d(t["scan_pos_angle"][i])[c])
            rows["pf_al"].append(t["parallax_factor_al"][i])
            rows["w"].append(w[c])
            rows["err"].append(np.atleast_1d(t["centroid_pos_error_al"][i])[c])
            rows["ipd_err"].append(np.atleast_1d(t["ipd_error_al"][i])[c])
            rows["excess_src"].append(t["agis_source_excess_noise"][i])
            ua = np.atleast_1d(t["used_by_agis_al"][i])
            rows["used"].append(bool(ua[c]) if c < len(ua) else False)
            rows["g_mag"].append(t["g_mag"][i])
            rows["ra0"].append(t["ra0"][i])
            rows["dec0"].append(t["dec0"][i])
    return {k: np.array(v) for k, v in rows.items()}

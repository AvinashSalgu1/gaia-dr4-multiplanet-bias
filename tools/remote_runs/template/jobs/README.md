# Jobs

One file per run, named with a number so they run in order: `0001-test.toml`, `0002-step3-all-stars.toml`, ...
Pushing a new file here starts it on the server. A job runs once: when `results/<name>/` exists it
is skipped. To run it again, delete that folder or add a new job file with a new name.

```toml
name      = "0002-step3-all-stars"     # folder name under results/
code_ref  = "main"                     # branch, tag or commit of gedi
cores     = 64                         # at most MAX_CORES (runner/config.sh)
max_hours = 100                        # stopped after this (optional, default 100; GitHub's own limit is 5 days)
command   = "python scripts/step3.py --jobs {cores}"   # runs in the code folder; {cores} is replaced
collect   = ["results/step3_gaia/*.txt", "results/step3_gaia/*.png"]   # files or folders to copy back
```

Only small files come back (txt, csv, png, pdf, svg, json, md, dat, log; at most 20 MB each and
100 MB per job). Anything else stays on the server and is listed in `results/<name>/status.txt`.
Several commands can be chained in `command` with `&&`.

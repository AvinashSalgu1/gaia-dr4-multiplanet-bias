# Running the code on a remote server through GitHub

For long runs on a shared server without a queue system (a department workstation, for example):
you push a small job file to a **private** GitHub repository, a GitHub Actions runner on the server
picks it up within seconds, runs it, and pushes back the summaries, tables and figures. Large outputs
stay on the server. You do not need to log in to the server to start or follow a run; the live output
is on the repository's Actions tab, and GitHub e-mails you when a run ends.

Not for clusters with a queue system (SLURM, PBS): their rules usually forbid background services,
and jobs go through the queue instead (see ../slurm/).

## What is in template/

```
.github/workflows/run-jobs.yml  starts a sweep when a job file is pushed (or with "Run workflow")
runner/config.sh                paths, Python environment, code repository, runner name, core limit
runner/setup_server.sh          one-time setup on the server
runner/sweep.sh                 runs every job without a results folder, oldest first, and pushes results
runner/run_job.py               runs one job and copies back its small output files
jobs/README.md                  job file format
examples/0001-test.toml         first test job (tests + gaiamock check)
```

## Set-up (about 20 minutes)

1. Create a new **private** repository on GitHub (empty) and copy the contents of `template/` into it.
   Edit `runner/config.sh` (paths on the server, your code repository, the runner name and label) and
   push. Keep it private: the runner executes whatever is pushed to `jobs/`.
2. On the server, install Miniforge in your space and create the environment
   (`conda create -n gaia python=3.12 numpy scipy matplotlib astropy healpy joblib astroquery pytest rebound dynesty gsl`).
3. Copy the new repository to the server (it is private: clone it on a computer that has access and
   `rsync` it over), then on the server:

   ```bash
   bash runner/setup_server.sh key       # only if the code repository is private: prints a key to add as a
                                         # read-only deploy key of the code repository
   bash runner/setup_server.sh code      # clones the code and gaiamock and builds gaiamock against conda's GSL
   # gaiamock's scan files (49,152 files, ~1 GB) go into $GAIAMOCK_DIR/healpix_scans/
   bash runner/setup_server.sh runner TOKEN   # TOKEN: Settings > Actions > Runners > New self-hosted runner
   bash runner/setup_server.sh status
   ```

   The runner keeps going in a `tmux` session and restarts after a reboot (`crontab @reboot`).
4. Copy `examples/0001-test.toml` into `jobs/` and push. Within a minute the run appears on the
   Actions tab; its results arrive in `results/0001-test/`.

## Jobs

A job is a TOML file in `jobs/` (format in `jobs/README.md`): a name, the code version (branch, tag or
commit), the number of cores, a time limit, the command (run in the code folder; `{cores}` is replaced)
and the files to bring back. Jobs run one at a time, oldest file name first; a job runs once (delete its
results folder to rerun). The command runs at normal priority (`nice 0`; set `nice = 1..19` in a job
file to yield to other users), with every parallel library limited to the cores asked for. A copy of every results folder is also
kept on the server (`BACKUP_DIR`), and a run whose results cannot be pushed is marked as failed.

While a job runs, its `results/<name>/status.txt` is updated and pushed every hour
(`PROGRESS_MINUTES`): time so far, the last `k/N` progress line of the output, the share done and the
estimated time left of the current step. Scripts run this way should print such counts (`k/N`) with
`flush=True`. On a shared machine without a queue, run times depend on other users' load; plan with a
margin.

Push several job files in one commit: each push starts a workflow run, and GitHub keeps only one
waiting run per queue, replacing (and reporting as "cancelled") any earlier waiting one. Nothing is
lost, because the running sweep picks up every job file pushed while it was busy.

## Limits

- One GitHub job runs for at most 5 days; split longer runs (most scripts have `--part k/n`) or resume.
- Only files under 20 MB (100 MB per job) are copied back.
- Ask the server's administrator before installing a background service on a shared machine.

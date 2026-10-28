# Settings for the runner on your server. Copy this repository template, then edit these values.
export BASE=/path/to/your/space                              # a disk with room for data and results
export CONDA_DIR=$BASE/miniforge3                            # Python environment manager (see README)
export CONDA_ENV=gaia                                        # environment with the packages of the code
export CODE_REPO=OWNER/gedi                                 # GitHub repository with the code (yours or a fork)
export CODE_PRIVATE=no                                       # yes: clone it with a read-only deploy key
export RUNS_REPO_URL=https://github.com/OWNER/gaia-dr4-runs  # this repository (keep it PRIVATE)
export DEV_DIR=$BASE/gedi                                   # clone of the code on the server
export GAIAMOCK_DIR=$BASE/gaiamock
export BACKUP_DIR=$BASE/runs-results                         # server copy of every results/<name>/ folder
export RUNNER_DIR=$BASE/actions-runner
export RUNNER_NAME=my-server
export RUNNER_LABEL=gaia-runner                              # must match runs-on in .github/workflows/run-jobs.yml
export MAX_CORES=64                                          # largest "cores" a job may ask for
export GIT_NAME="Your Name"                                  # author of the results commits
export GIT_EMAIL="you@example.org"

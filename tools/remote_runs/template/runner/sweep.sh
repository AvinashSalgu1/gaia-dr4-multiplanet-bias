#!/usr/bin/env bash
# Run every job in jobs/*.toml that has no results/<name>/ folder yet, oldest file name first.
# After each job the results are committed and pushed, so they appear on GitHub one by one.
# To run a job again, delete its results/<name>/ folder (or copy the job under a new name).
set -uo pipefail
{   # the whole script is read before it runs, so a `git pull` that updates this file cannot disturb it
cd "$(dirname "$0")/.."
source runner/config.sh
source "$CONDA_DIR/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
# identity for every git operation here (a rebase in `git pull --rebase` rewrites commits and needs it too)
git config user.name "$GIT_NAME"
git config user.email "$GIT_EMAIL"
# never wait for ever on the network: no password prompts, abort transfers that stall for 60 s,
# and give every pull or push at most 5 minutes (a hung push once blocked the queue for a day)
export GIT_TERMINAL_PROMPT=0 GIT_SSH_COMMAND="ssh -o BatchMode=yes -o ServerAliveInterval=30 -o ServerAliveCountMax=4"
git config http.lowSpeedLimit 1000
git config http.lowSpeedTime 60
g() { timeout 300 git "$@"; }

push() {
    git add -A results
    git commit -q -m "$1" || return 0
    for _ in 1 2 3 4 5; do
        g pull -q --rebase && g push -q && return 0
        git rebase --abort 2>/dev/null
        sleep 30
    done
    echo "::error::could not push '$1' to GitHub; the results are kept on the server in $BACKUP_DIR"
    failed=1                                  # the run then shows as failed, and GitHub emails
}

failed=0
shopt -s nullglob
declare -A bad=()
next_job() {        # first job file (by name) without a results folder, or nothing
    for job in jobs/*.toml; do
        [ -n "${bad[$job]:-}" ] && continue
        name=$(python -c "import tomllib,sys; print(tomllib.load(open(sys.argv[1],'rb'))['name'])" "$job") || {
            echo "skipping $job: cannot read it" >&2; bad[$job]=1; failed=1; continue; }
        [ -e "results/$name" ] || { echo "$job $name"; return; }
    done
}
# After every job, fetch the repository again: jobs pushed while one was running are run
# by this same sweep, so nothing waits for (or is lost with) a queued workflow run.
while true; do
    g pull -q --rebase || { git rebase --abort 2>/dev/null; echo "warning: could not fetch new jobs"; }
    read -r job name < <(next_job) || true
    [ -z "${job:-}" ] && break
    mkdir -p "results/$name"
    printf 'status:     running on %s since %s\n' "$(hostname)" "$(date -u '+%F %T UTC')" \
        > "results/$name/status.txt"
    push "Start: $name"
    echo "::group::$name"
    python runner/run_job.py "$job" || failed=1
    echo "::endgroup::"
    push "Results: $name ($(sed -n 's/^status: *//p' "results/$name/status.txt" | head -1))"
    job=""
done
exit $failed
}

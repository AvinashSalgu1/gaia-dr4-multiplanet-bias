#!/usr/bin/env bash
# One-time setup on the server. Edit runner/config.sh first, then run the parts in order:
#   bash runner/setup_server.sh key        # 1. (private code repository only) makes a read-only key
#   bash runner/setup_server.sh code       # 2. clones the code and gaiamock, builds gaiamock
#   bash runner/setup_server.sh runner TOKEN   # 3. installs and starts the GitHub runner
#   bash runner/setup_server.sh status     # check that everything is in place
set -euo pipefail
cd "$(dirname "$0")/.."
RUNS_DIR=$(pwd)
source runner/config.sh
RUNNER_VERSION=2.338.1          # https://github.com/actions/runner/releases
KEY=~/.ssh/gaia_code_readonly

case ${1:-} in
key)
    mkdir -p ~/.ssh && chmod 700 ~/.ssh
    [ -f "$KEY" ] || ssh-keygen -q -t ed25519 -N "" -C "$RUNNER_NAME read-only" -f "$KEY"
    if ! grep -q "Host github-code" ~/.ssh/config 2>/dev/null; then
        cat >> ~/.ssh/config <<EOF
Host github-code
    HostName ssh.github.com
    Port 443
    User git
    IdentityFile $KEY
    IdentitiesOnly yes
EOF
        chmod 600 ~/.ssh/config
    fi
    echo
    echo "Add this line as a read-only deploy key of $CODE_REPO"
    echo "(GitHub > that repository > Settings > Deploy keys > Add deploy key; leave 'Allow write access' OFF):"
    echo
    cat "$KEY.pub"
    ;;
code)
    source "$CONDA_DIR/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"
    if [ "$CODE_PRIVATE" = yes ]; then
        # ssh -T always exits with 1 on GitHub (no shell), so check its message instead of its exit code
        reply=$(ssh -o StrictHostKeyChecking=accept-new -T github-code 2>&1 || true)
        case $reply in
            *"successfully authenticated"*) ;;
            *) echo "The deploy key is not working yet (step 1). ssh said:"; echo "$reply"; exit 1 ;;
        esac
        url=git@github-code:$CODE_REPO.git
    else
        url=https://github.com/$CODE_REPO.git
    fi
    [ -d "$DEV_DIR/.git" ] || git clone "$url" "$DEV_DIR"
    pip install -q --no-deps -e "$DEV_DIR"      # makes `import gedi` work (editable: follows checkouts)
    [ -d "$GAIAMOCK_DIR/.git" ] || git clone https://github.com/kareemelbadry/gaiamock "$GAIAMOCK_DIR"
    (cd "$GAIAMOCK_DIR" && gcc -O2 -shared -fPIC -o kepler_solve_astrometry.so kepler_solve_astrometry.c \
        -I"$CONDA_PREFIX/include" -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" -lgsl -lgslcblas -lm)
    echo "Code and gaiamock ready. gaiamock still needs its scan files in $GAIAMOCK_DIR/healpix_scans/"
    echo "(download healpix_scans.zip as described in gaiamock's README, or copy them with rsync)."
    ;;
runner)
    TOKEN=${2:?give the registration token from GitHub > this repository > Settings > Actions > Runners > New self-hosted runner}
    mkdir -p "$RUNNER_DIR" && cd "$RUNNER_DIR"
    if [ ! -f config.sh ]; then
        curl -sSL -o runner.tar.gz \
            "https://github.com/actions/runner/releases/download/v$RUNNER_VERSION/actions-runner-linux-x64-$RUNNER_VERSION.tar.gz"
        tar xzf runner.tar.gz && rm runner.tar.gz
    fi
    ./config.sh --unattended --replace --url "$RUNS_REPO_URL" \
        --token "$TOKEN" --name "$RUNNER_NAME" --labels "$RUNNER_LABEL" --work _work
    # keep it running: start now in tmux, and again after every reboot
    tmux has-session -t runner 2>/dev/null || tmux new-session -d -s runner "cd $RUNNER_DIR && ./run.sh"
    line="@reboot tmux new-session -d -s runner 'cd $RUNNER_DIR && ./run.sh'"
    { crontab -l 2>/dev/null | grep -v "actions-runner" || true; echo "$line"; } | crontab -
    echo "Runner started (tmux session 'runner'). It should now show as Idle on GitHub."
    ;;
status)
    echo -n "deploy key:  "; [ -f "$KEY" ] && echo ok || echo missing
    echo -n "code:        "; [ -d "$DEV_DIR/.git" ] && git -C "$DEV_DIR" log --oneline -1 || echo missing
    echo -n "gaiamock:    "; [ -f "$GAIAMOCK_DIR/kepler_solve_astrometry.so" ] && echo built || echo "not built"
    echo -n "scan files:  "; ls "$GAIAMOCK_DIR/healpix_scans" 2>/dev/null | wc -l; echo "             (should be 49152)"
    echo -n "runner:      "; tmux has-session -t runner 2>/dev/null && echo running || echo "not running"
    echo -n "autostart:   "; crontab -l 2>/dev/null | grep -q actions-runner && echo "on (after reboot)" || echo "off"
    ;;
*)
    sed -n '2,7p' "$0"
    ;;
esac

#!/usr/bin/env bash
# Extended analyses that make new fits or integrations (the g2_* scripts). Every step resumes where
# it stopped, so the script can simply be started again after an interruption.
# Logs go to results/extended/logs/.
#
#   bash scripts/extended/run_extended.sh all                  # every step, in order
#   bash scripts/extended/run_extended.sh census evidence ...  # chosen steps
#
# Steps: table (needs results/step2_population.csv), noise, evidence, evidence3, evidence03, truth
# (after the evidence runs), census (needs internet), stability, stabilityall, observability,
# mixed (needs results/step3_gaia). Run from anywhere with the virtual environment active.
set -euo pipefail
cd "$(dirname "$0")/../.."
mkdir -p results/extended/logs
run() {
    local name=$1; shift
    echo "=== $name started $(date) on $(hostname) ===" | tee -a results/extended/logs/progress.log
    python "$@" 2>&1 | tee "results/extended/logs/$name.log"
    echo "=== $name finished $(date) ===" | tee -a results/extended/logs/progress.log
}
step() {
    case $1 in
        table)         run table         scripts/extended/g2_step2_exclusive.py ;;
        noise)         run noise         scripts/extended/g2_noise_robustness.py ;;
        evidence)      run evidence      scripts/extended/g2_evidence.py ;;
        evidence3)     run evidence3     scripts/extended/g2_evidence.py --prior-scale 3 ;;
        evidence03)    run evidence03    scripts/extended/g2_evidence.py --prior-scale 0.3 ;;
        truth)         run truth         scripts/extended/g2_evidence_truth.py ;;
        census)        run census        scripts/extended/g2_census.py ;;
        stability)     run stability     scripts/extended/g2_stability_nbody.py ;;
        stabilityall)  run stabilityall  scripts/extended/g2_stability_nbody.py --all ;;
        observability) run observability scripts/extended/g2_observability.py ;;
        mixed)         run mixed         scripts/extended/g2_mixed_flags.py ;;
        *) echo "unknown step: $1"; exit 1 ;;
    esac
}
case ${1:-all} in
    all) steps=(table noise evidence evidence3 evidence03 truth census stability stabilityall observability mixed) ;;
    *)   steps=("$@") ;;
esac
for s in "${steps[@]}"; do step "$s"; done
echo "Done: ${steps[*]}."

#!/bin/sh
# Reproduce every figure and number of the paper from the data in data/.
# Printed results go to outputs/*.txt and figures to outputs/figures/.
set -e
cd "$(dirname "$0")"
mkdir -p outputs/figures
for s in energy_fits rms_residuals waic_indistinguishability convention_refit hashrate_series \
         nb_gate_count landauer_estimate intermediate_forms make_errorbar_figures \
         make_fig3_deltaP make_fig4_ksn; do
    echo "== $s.py"
    case "$s" in
        waic_indistinguishability) python3 "$s.py" --both > "outputs/$s.txt" ;;
        nb_gate_count) python3 "$s.py" 16384 > "outputs/$s.txt" ;;  # the paper's N_b
        *) python3 "$s.py" > "outputs/$s.txt" ;;
    esac
done
echo "Done: see outputs/."

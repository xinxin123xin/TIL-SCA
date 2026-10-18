#!/usr/bin/env python3
"""Recompute the additional validation evidence released with the TIL-SCA manuscript.

The script uses only the Python standard library and the CSV/JSON files under
``data/additional_checks``.  It re-derives the numbers quoted in the manuscript:
the 896-position family-wise null, the single-event injection sensitivity,
the per-entry sec_add statistic, the CP/plaintext-HW relation, the Verilator
baseline ratios, the split-storage PPA proxy, and the VCD crossing check.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'additional_checks'


def read(name):
    with open(DATA / name, newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def percentile(values, q):
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def check(label, condition, detail):
    print(f"{'PASS' if condition else 'FAIL'}  {label}: {detail}")
    return condition


def main():
    ok = True

    # --- 896-position family-wise null -------------------------------------
    null = [float(row['max_abs_t']) for row in read('null_maxabs.csv')]
    crossings = sum(value >= 4.5 for value in null)
    p95 = percentile(null, 0.95)
    p99 = percentile(null, 0.99)
    ok &= check('null draws=200', len(null) == 200, f'{len(null)} draws')
    ok &= check('null max<4.5', max(null) < 4.5, f'max={max(null):.2f}')
    ok &= check('null p95~3.92', abs(p95 - 3.92) < 0.05, f'p95={p95:.2f}')
    ok &= check('null p99~4.17', abs(p99 - 4.17) < 0.05, f'p99={p99:.2f}')
    ok &= check('null no crossing', crossings == 0, f'{crossings}/200 above 4.5')

    # --- single-event injection -------------------------------------------
    sweep = read('injection_sweep.csv')
    deltas = [float(row['delta_bits']) for row in sweep]
    detection = [float(row['detection_prob']) for row in sweep]
    index = next(i for i in range(len(detection) - 1)
                 if detection[i] < 0.8 <= detection[i + 1])
    delta80 = deltas[index] + (0.8 - detection[index]) * (
        deltas[index + 1] - deltas[index]) / (detection[index + 1] - detection[index])
    at_half = next(float(row['detection_prob']) for row in sweep
                   if abs(float(row['delta_bits']) - 0.5) < 1e-9)
    at_half_loc = next(float(row['localization_prob']) for row in sweep
                       if abs(float(row['delta_bits']) - 0.5) < 1e-9)
    ok &= check('injection delta80~0.48', abs(delta80 - 0.48) < 0.02,
                f'delta80={delta80:.3f} bit')
    ok &= check('injection 0.5-bit detection', abs(at_half - 0.85) < 0.01,
                f'{at_half:.2f}')
    ok &= check('injection 0.5-bit localization', abs(at_half_loc - 0.975) < 0.01,
                f'{at_half_loc:.3f}')

    # --- per-entry sec_add statistic --------------------------------------
    entries = read('sec_add_per_entry_t.csv')
    seed_max = {}
    for row in entries:
        if row['entry'] != 'entry3':
            continue
        seed = int(row['seed'])
        seed_max[seed] = max(seed_max.get(seed, 0.0), float(row['abs_t']))
    low, high = min(seed_max.values()), max(seed_max.values())
    averages = [float(row['abs_t']) for row in entries
                if row['entry'] == 'all_mean_of_5' and row['register'] == 'g0_r']
    ok &= check('per-entry range', abs(low - 97.61) < 0.1 and abs(high - 101.30) < 0.1,
                f'entry3={low:.2f}-{high:.2f}')
    ok &= check('five-entry average lower', 48.0 < min(averages) < 51.0,
                f'average={min(averages):.2f}-{max(averages):.2f}')

    # --- CP versus plaintext Hamming weight -------------------------------
    cp_rows = read('cp_plaintext_hw.csv')
    cp_ok = True
    for row in cp_rows:
        predicted = float(row['predicted_mean_fixed'])
        observed = float(row['CP_mean_fixed'])
        cp_ok &= abs(predicted - observed) < 0.4
    ok &= check('CP follows (W-2HW)/4', cp_ok,
                'all rows within 0.4 of the predicted mean')
    a2b = next(row for row in cp_rows if row['case'].startswith('a2b'))
    ok &= check('a2b CP near zero', abs(float(a2b['CP_t'])) < 1.0,
                f"t={float(a2b['CP_t']):.2f}")

    # --- Verilator baseline ------------------------------------------------
    baseline = read('verilator_baseline.csv')
    for row in baseline:
        if row['design'] == 'AES-128':
            continue
        ratio = float(row['oim_s']) / float(row['verilator_s'])
        ok &= check(f"{row['design']} ratio", abs(ratio - float(row['ratio'])) / ratio < 0.05,
                    f'{ratio:.1f}x')

    # --- split-storage PPA proxy ------------------------------------------
    ppa = {row['metric']: row for row in read('ppa_split_storage.csv')}
    ok &= check('PPA cells +12', int(ppa['cells']['delta']) == 12,
                f"delta={ppa['cells']['delta']}")
    ok &= check('PPA FF +66', int(ppa['ff_bits']['delta']) == 66,
                f"delta={ppa['ff_bits']['delta']}")

    # --- VCD crossing check ------------------------------------------------
    summary = json.loads((DATA / 'vcd_check_summary.json').read_text(encoding='utf-8-sig'))
    ok &= check('VCD events', summary['events'] == 896, f"{summary['events']} events")
    ok &= check('VCD HD mismatches', summary['hd_mismatch'] == 0,
                f"{summary['hd_mismatch']} of 896")
    ok &= check('VCD result mismatches', summary['result_mismatch'] == 0,
                f"{summary['result_mismatch']} of 896")

    print('\nRESULT:', 'PASS' if ok else 'FAIL')
    raise SystemExit(0 if ok else 1)


if __name__ == '__main__':
    main()




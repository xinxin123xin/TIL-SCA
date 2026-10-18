# Additional validation

This directory documents the additional checks added in a later revision of
**TIL-SCA: In-Loop Power Side-Channel Analysis with RTL Operation Context for
PQC Hardware**. The released CSVs are the exact inputs used by
`scripts/reproduce_additional_checks.py`; the script re-derives every number
quoted in the manuscript and exits non-zero on any mismatch.

## Files

| File | Contents | Manuscript statement |
|---|---|---|
| `data/additional_checks/null_maxabs.csv` | 200 label-permutation draws of `max|t|` over the 896 NTT positions | Family-wise null: max 4.43, p99 4.17, 0/200 above 4.5 |
| `data/additional_checks/injection_sweep.csv` | Detection/localization probability versus injected mean shift (0.1--1.5 bit) | `delta_80` = 0.48 bit; 85% detection and 97.5% localization at 0.5 bit |
| `data/additional_checks/sec_add_per_entry_t.csv` | Three seeds, five FSM entries, both `g0_r`/`g1_r` | Entry 3 = 97.61--101.30; five-entry average 48.2--50.8 |
| `data/additional_checks/cp_plaintext_hw.csv` | Predicted `(W-2*HW)/4` versus observed CP mean | CP mean follows the fixed plaintext Hamming weight; `a2b_core` CP near zero because `x=1234` is near-balanced |
| `data/additional_checks/verilator_baseline.csv` | OIM/tsim versus Verilator wall time, N=100 | 10x CBD2-256, 14x NTT-256, ~980x short-transaction AES boundary case |
| `data/additional_checks/ppa_split_storage.csv` | Generic-cell and flip-flop counts before/after the two-bank split | +66 FF bits, +12 cells (+3.3%), no timing-closure claim |
| `data/additional_checks/vcd_check_summary.json` | Verilator VCD check over one NTT trace | 0/896 HD, result, and zeta mismatches |
| `data/additional_checks/vcd_crossing_check.csv` | Five sampled crossings with operands, result, and HD | Crossings confirmed on the RTL waveform |

## Methods in brief

**Null distribution.** The NTT-256 feature matrix (5000 traces x 896 positions)
is relabelled 200 times (2500/2500) and the Welch statistic is recomputed at
every position. The family-wise maximum never exceeds 4.5, so the screening
threshold has no observed false positive at this scale.

**Injection sensitivity.** A mean shift is added to one randomly chosen event of
the permuted feature stream and the detection/locus decision is recomputed over
the remaining 895 positions. The 80% point is 0.48 bit; at 0.5 bit the injected
event is detected in 85% of trials and localized in 97.5%. This is a
feature-stream injection, not an RTL-level leak injection.

**Per-entry statistic.** The `sec_add` profile stores the five `state[3]` entry
HD values per trace. Reporting the third entry instead of the trace-level mean
raises the statistic from 48.2--50.8 to 97.61--101.30 and shows that the leak is
concentrated in one entry.

**CP and plaintext Hamming weight.** For Boolean two-share encoding, the per-bit
centered product has expectation `(W-2*HW(x))/4`; the observed CP means match
this prediction (within 0.4) for `sec_and`, `sec_add`, and `a2b_core`. A balanced
fixed plaintext removes the inherent term, so only a residual crossing then
indicates a defect beyond the second-order dependence of a two-share design.

**Verilator baseline.** The unchanged RTL is run through the existing direct-WPI
SystemVerilog testbenches with Verilator 5.050 (`--binary --timing`) on the same
frozen N=100 stimuli. OIM/tsim is 10x (CBD2-256) and 14x (NTT-256) slower; the
AES-128 number is the short-transaction boundary case from the earlier benchmark.
Verilator testbenches use `#1` delays and a 10 ns clock, so the ratios are a
lower bound on the gap.

**PPA proxy.** Yosys 0.65 `read_verilog; hierarchy; flatten; proc; opt; stat` on
the baseline and split-storage coprocessor tops. Values are generic cells and FF
bits, not library-mapped area or `F_max`.

**VCD crossing check.** One NTT trace is run under Verilator with `--trace`; for
each of the 896 `result_valid` pulses the script reads the operands, the return
value, and the previous return, reconstructs
`HD = popcount(result XOR previous)`, and compares it with the OIM in-loop
feature matrix row 0. All 896 HD values, results, and twiddle values match. Five
crossings (events 0, 200, 400, 513, 700) are released with their operands and
records. The waveform confirms the operands, return value, and observation; the
stage/j labels remain schedule metadata because the coprocessor executes a flat
CT_MUL stream and holds no NTT-stage state.

## Reproduce

```bash
python scripts/reproduce_additional_checks.py
```

## Provenance

The source runs are archived in the private experiment tree: the null and
injection campaign, the per-entry tables, the CP re-analysis, the Yosys PPA logs,
the CBD2/NTT Verilator baselines, and the VCD dump with its parser. The released
CSVs are derived from those runs and are sufficient to reproduce every reported
number.

The full VCD is 4.47 MB and is not bundled; the released summary and the five
sampled crossings are sufficient to recheck the reported claims.

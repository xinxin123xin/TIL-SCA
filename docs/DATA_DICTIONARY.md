# Released data and their connection to the paper

All trace-level values below are digital RTL activity features. Counts and HW/HD
values are in bits. Every campaign has 5,000 traces: 2,500 fixed and 2,500 random.
Seed campaigns are analyzed separately. `data/protocol.json` supplies input
conditions, centering, threshold, widths, snapshot windows, and coverage.

| File | Rows | Evidence |
|---|---:|---|
| `masked_outputs.csv` | 75,000 | Output HW/CP in Table II; nine adder seed/operand campaigns |
| `sec_add_entry_profiles.csv` | 30,000 | Adder internal HD in Table II and signed mean-HD profiles in Fig. 2 |
| `a2b_internal_hd.csv` | 15,000 | A2B internal HD in Table II |
| `intervention_traces.csv` | 30,000 | Six matched baseline/split-storage campaigns |
| `intervention_statistics.csv` | 60 | Original in-loop/offline statistics and statuses for those six campaigns |
| `pqc_paired_comparisons.csv` | 1,152 | Table I: 256 CBD2 and 896 NTT paired event summaries |
| `sec_add_paired_checks.csv` | 9 | Archived second-order in-loop/offline discrepancies across adder conditions |
| `masked_paired_checks.csv` | 6 | Archived CP discrepancies for three AND and three A2B campaigns |
| `campaigns.csv` | 15 | Module, seed, condition, and class sizes for output records |

## Masked outputs

`campaign` names combine module, operand condition, and seed. `trace_id` is
the trace index within that campaign; `group` is `fixed` or `random`.
`L0` and `L1` are HW of the two predefined output shares, serialized as decimal
numbers. AND uses 32-bit shares; adder and A2B output HW use 13 bits.
`internal_hd` is populated only for AND, where it is the sum of HD over four
pipeline-register banks between cycles 1 and 2; every released value is zero.
Blank values for other modules mean that this column is not their internal
observation.

The CP feature is `(L0-pooled_mean_L0)*(L1-pooled_mean_L1)`. Each pooled mean
uses both fixed and random samples from that campaign. Both first-order share
tests and the second-order CP test use one value per trace. Table II's HW
number for a seed is the maximum absolute statistic of its two shares.
Its range is then taken over the three seeds, using the all-ones adder condition.

## Adder entry profiles

`g0_e1` through `g0_e5` and `g1_e1` through `g1_e5` hold the original integer
HD observations of `g0_r` and `g1_r` at five entries to `state[3]`.
The first compares the preceding `state[2]` entry; subsequent entries compare
the preceding `state[6]` entry. The profiles are pivoted from long-form source
records into one row per trace. No values are averaged in this released file.

For Table II, average the five entries within each trace and channel, compute
Welch statistics across traces, and take the maximum over `g0_r` and `g1_r`.
For Fig. 2, compute the fixed-minus-random mean of `g0_r` at each entry separately.
The file includes all-ones and all-zero conditions, three seeds each. The five
entries are repeated observations within a trace, rather than five independent
sampling units. These FSM records were retained for post-run diagnosis.

## A2B internal boundary

`raw0_hd` and `raw1_hd` compare the respective register snapshots between
`ST_ADD1_RUN` and entry to `ST_FILL_RND2` (state 4). There is one row per trace;
Table II reports `raw0_r`. The evaluated A2B RTL exposes an external
`rnd_bypass_i` interface; its adaptation patch is included. Internal-HD records
are from the post-run FSM diagnosis.

## Matched interface records

`variant` is `baseline` or `split_storage`; `seed` selects an independent pool.
`input_hash`, `output_hash`, and `protocol_hash` are experiment matching
digests, not credentials. `x`, `y`, `share0`, and `share1` are simulated
32-bit words in hexadecimal. `cycles` is the total transaction cycle count;
the three `*_ok` flags record functionality, protocol, and endpoint-identity
checks.

Feature names describe their boundaries:

- `result_reg_between_returns`: HD of the original shared register between
  the two share-return snapshots;
- `share_return0_write`, `share_return1_write`: HD on each new bank's write;
- `share_return0_between_returns`, `share_return1_between_returns`: bank
  snapshot HD across the share-return interval;
- `return_bus_valid_endpoint`: HD of the two valid returned words;
- `return_bus_adjacent_sum`: sum of sampled adjacent-cycle bus HD over that
  interval;
- `share0_hw`, `share1_hw`: HW of returned shares, also the predefined CP pair.

Each trace resets the DUT. A new bank therefore writes from zero to one share.
Returns occur at cycles 22 and 25, with zero bus values between them, and each
transaction lasts 26 cycles. Full cycle-by-cycle logs are omitted; the released
trace audits retain the resulting features and matching evidence. Rebuilding
the interface experiment generates the detailed cycle logs again.

`nan` means that a feature is absent or that Welch's denominator is zero.
`intervention_statistics.csv` distinguishes `not_present`, `undefined`, and
`finite`. A constant observation produces an undefined statistic, not `t=0`.
The maximum bank-write statistic 1.397 is a different observation from the
baseline shared-register statistic near 494.

## PQC paired summaries and numerical precision

`t_offline` and `t_inloop` are the originally serialized paired statistics.
`reported_abs_delta` preserves the original discrepancy column. In particular,
NTT's printed statistics have fewer digits than the accumulator used to
compute its discrepancy. Re-subtracting the printed `t` columns does not
recover that discrepancy. `scripts/reproduce.py` checks the archived residual
bound and coverage; it does not reconstruct PQC kernel statistics from a
trace-wise matrix.

The compact release omits the PQC trace-wise matrices and NTT input/event pools.
Their original drivers are included for inspection and integration. CBD2's
15,000 dedicated traces feed one accumulated reducer in three batches. NTT's
5,000 by 896 observations total 4,480,000; the first-order reducer retains
5,376 logical statistical scalars, excluding other application state.

The archive's masked-pair residual files also contain originally reported
discrepancies. The separate C replay check performs a fresh comparison against
two-pass Python with an absolute tolerance of `1e-9` to accommodate accumulation
order, while preserving the paper's archived numerical results.

## Provenance and reuse

`provenance/sources.csv` maps every copied or transformed source to its original
workspace-relative path and raw SHA-256. Transformations select columns, add
campaign labels, concatenate campaigns, or pivot entry rows. `SHA256SUMS`
anchors the released files. The generated experiment data follow the root MIT
license. A release tag should be used when referencing this artifact.

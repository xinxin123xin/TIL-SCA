# Release validation

The compact package was checked on Windows on 2026-10-11 using Python 3.14.7,
GCC 15.2.0 (w64devkit), and Yosys 0.66. Python reanalysis uses only the standard
library. The archived paper environment is recorded separately in
`data/protocol.json`.

## Released data and extracted reducers

`python scripts/reproduce.py --check-core` passed for:

- 15 masked-output campaigns, comprising 75,000 trace-level feature rows;
- 30,000 five-entry adder profiles and 15,000 A2B internal-boundary records;
- six matched interface campaigns, comprising 30,000 trace-level records;
- archived PQC paired residuals and observation coverage.

The extracted C first-order and centered-product reducers agreed with
independent two-pass Python calculations to a maximum absolute difference
of `8.853362487570848e-12` across the 15 masked-output campaigns. This replay
comparison is distinct from the paper's archived in-loop/offline comparisons.
All three masked modules had an indication beyond share-local output HW.
The interface checks verified returned-share identities, functional results,
cycle counts, and matching input/output/protocol hashes across variants.

PQC verification checks the supplied per-event summaries and archived
residuals; the compact package does not recompute PQC statistics from full
trace-wise feature matrices.

## Actual RTL smoke tests

All seven configurations built from the separately obtained, hash-checked WPI
RTL and ran successfully with ten samples per class, seed 0:

| Configuration | Result |
|---|---|
| `sec_and` | PASS |
| `sec_add` | PASS |
| `sec_add_entries` (`allones`) | PASS |
| `a2b` | PASS |
| `a2b_entries` | PASS |
| `interface_baseline` | PASS |
| `interface_split` | PASS |

The interface variants each passed 20/20 functional, protocol, and endpoint
identity checks; their operands, returned shares, cycle counts, and audit
hashes matched trace by trace. These small runs check the build and execution
path; paper-size screening conclusions come from the released 5,000-trace
campaigns. Full PQC execution is outside this compact build entry point.

## Windows compiler selection

The compiler was selected through `CC`, with another GCC installation earlier
on the parent process PATH. Both Python entry points successfully compiled
using the selected compiler's own `bin` directory at the front of their
compiler subprocess PATH. See `REPRODUCING_RTL.md` for the PowerShell setup.

# TIL-SCA: Companion Code and Data

This repository contains the core source code and selected experimental data
supporting the manuscript **TIL-SCA: In-Loop Power Side-Channel Analysis with
RTL Operation Context for PQC Hardware**, by Xin Zhou, Yijun Cui, Chenghua Wang,
and Weiqiang Liu. Publication metadata will be added when available.

Tensor-Based In-Loop Side-Channel Analysis (TIL-SCA) embeds first-order and
preselected-pair second-order TVLA reduction into tensor-based RTL execution.
Event hooks retain RTL signal and event identities for examining output shares,
internal-register transitions, and interface storage.

## What is included

| Directory | Contents |
|---|---|
| `src/` | Welford and nine-moment centered-product routines extracted verbatim from the experiment driver |
| `drivers/` | Actual CBD2, NTT, masked-output, FSM-entry, and storage-intervention drivers |
| `wrappers/`, `recipes/`, `patches/` | Observation wrappers, per-design Yosys recipes, evaluated RTL adaptations, and split-return-storage intervention |
| `third_party/wpi/` | Unchanged ISC-licensed bit-level `tsim.h`, OIM generator, and original license |
| `data/` | Selected trace-level RTL features, paired statistics, and experiment configuration |
| `data/additional_checks/` | Null distribution, injection sensitivity, per-entry statistics, CP relation, Verilator baseline, PPA proxy, and VCD crossing check (see `docs/ADDITIONAL_CHECKS.md`) |
| `scripts/` | Feature-level result reproduction and optional RTL rebuilding |
| `provenance/` | Original workspace-relative source names, SHA-256 hashes, transformations, and upstream file identities |

The release keeps the main paper evidence: 75,000 masked-output records,
30,000 five-entry adder profiles, 15,000 A2B boundary records, and 30,000 matched
interface records. PQC data contain 1,152 per-event paired summaries (256 CBD2
positions and 896 NTT positions). This is a compact companion artifact; the
PQC trace-wise matrices and complete project history are outside this release.

## Quick start: recompute the released results

Python 3.10 or newer and its standard library are sufficient:

```bash
python scripts/reproduce.py
```

On systems where Python is named `python3`, substitute that command. Results
are written to `results/`:

- `table1.csv`: archived PQC paired residuals and event coverage;
- `table2.csv`: recomputed output-HW, output-CP, and internal-HD ranges;
- `figure2_points.csv`: signed mean-HD contrasts for all five adder entries;
- `intervention_recomputed.csv`: recomputed interface statistics;
- `verification.json`: numerical, functional-identity, and matching checks.

## Additional validation

`data/additional_checks/` and `scripts/reproduce_additional_checks.py` collect
the additional checks added in a later revision of the manuscript: the
896-position family-wise null (max `|t|` = 4.43, 0/200 above 4.5), the
single-event injection sensitivity (`delta_80` = 0.48 bit; 85% detection and
97.5% localization at 0.5 bit), the per-entry `sec_add` statistic (entry 3 =
97.61--101.30), the CP/plaintext-HW relation, the Verilator baseline (10x CBD2,
14x NTT, ~980x AES boundary case), the split-storage PPA proxy (+66 FF bits,
+12 cells), and the VCD crossing check (0/896 HD, result, and zeta mismatches).
Run:

```bash
python scripts/reproduce_additional_checks.py
```

The script uses only the standard library and exits non-zero if any released
number does not reproduce. See `docs/ADDITIONAL_CHECKS.md` for provenance and
limitations.

With a GCC/Clang C compiler, also check the extracted C accumulators against
independent two-pass Python statistics:

```bash
python scripts/reproduce.py --check-core --cc gcc
```

This replays archived features through the actual accumulator routines.
The original execution-integrated implementation is in `drivers/`.
Both scripts accept `--cc` and the `CC` environment variable. They prepend the
selected compiler's `bin` directory to its subprocess PATH; see the Windows
example in [RTL reproduction](docs/REPRODUCING_RTL.md).

`SHA256SUMS` covers source, data, provenance, and license files. After an
intentional change to these files, run `python scripts/update_checksums.py`
and rerun reproduction before releasing an update. README, citation metadata,
and documentation can be edited without changing this manifest.

## Run an RTL experiment

Obtain the [WPI reference source](https://github.com/Secure-Embedded-Systems/open-source-reference-for-reproducible-side-channel-assessment-of-Post-Quantum-Standards)
separately. The source tree must match the file identities recorded in
`provenance/rtl_sources.json`. With Yosys and GCC/Clang available:

```bash
python scripts/build_rtl.py --case sec_and --upstream /path/to/wpi-reference
```

The script prints the build directory and a 20-trace smoke command. Other
cases cover adder/A2B outputs and FSM entries, and both interface variants.
See [RTL reproduction](docs/REPRODUCING_RTL.md) for commands, full campaign
sizes, and the separate requirements of the archived PQC drivers.

## Interpreting the evidence

Across the three evaluated masked modules, share-local output HW remained
below the screening threshold, while output CP or internal HD crossed in
every module. Statistical order and RTL observation boundary therefore
provide complementary assessment dimensions.

The matched, reset-per-trace storage experiment removes the monitored
shared-register overwrite (baseline `|t|=488.615–494.412`). The original
register then has an undefined statistic because it is constant. The new
bank-write tests have maximum `|t|=1.397`; output-pair CP persists at
`|t|=41.220–42.775`. These are different observations, and the valid-return
endpoint relation also persists.

HW/HD are RTL activity proxies. The screening threshold is `|t| >= 4.5`, with
preselected pairs and no multiple-testing correction. A trace is the sampling
unit, including when it supplies five FSM-entry observations. Detailed fields,
centering, snapshot windows, and undefined values are defined in the
[data dictionary](docs/DATA_DICTIONARY.md).

The [release validation record](docs/VALIDATION.md) describes the feature-level
reanalysis, extracted-reducer checks, and seven successful RTL smoke tests.

## License and citation

TIL-SCA-authored materials use the MIT license; the bundled WPI runtime and
generator retain their ISC license. See [third-party notices](THIRD_PARTY_NOTICES.md).
The source data are RTL simulation features rather than measured power traces.

Please cite the manuscript and the release used. `CITATION.cff` supplies the
current title and authors without assuming an acceptance status or DOI.

## 中文说明

这是支撑上述 TIL-SCA 论文的精简代码与数据包。先运行
`python scripts/reproduce.py` 可重算论文主要 masked-RTL 结果、图 2 数据点
及接口干预结果；增加 `--check-core` 可验证 C 归约器。PQC 部分提供已归档的
逐 event 配对统计，完整 PQC 特征矩阵未随包提供。实际 RTL 小规模复现另需
Yosys、GCC/Clang 和单独取得的 WPI RTL，详见 `docs/REPRODUCING_RTL.md`。

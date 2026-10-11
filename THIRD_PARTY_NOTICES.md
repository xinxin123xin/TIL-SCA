# Third-party sources and licenses

The bit-level RTL execution path comes from the WPI reference artifact. The
tensor-based execution formulation is described in RTeAAL Sim; its Scala,
Java/SBT, FIRRTL, and benchmark sources are not required by this release.

## Bundled WPI components

Source repository:
https://github.com/Secure-Embedded-Systems/open-source-reference-for-reproducible-side-channel-assessment-of-Post-Quantum-Standards

| Released file | Upstream path | SHA-256 |
|---|---|---|
| `third_party/wpi/tsim.h` | `pre-silicon-mlkem-tvla/src/tsim_kernel/tsim.h` | `fe945d84b70689b2cb32f9e5c49e80126cff744569f9a82858fe914a41863ad4` |
| `third_party/wpi/gen_oim_fullchip.py` | `pre-silicon-mlkem-tvla/src/oim/gen_oim_fullchip.py` | `ad000cee48d2ff757cec0f7d442e05a008c2e5f42bab626c000077442bc4f0b8` |

These files are copied unchanged. The complete original ISC license,
including its scope statement and copyright attribution to Dillibabu
Shanmugam, Patrick Schaumont, and Worcester Polytechnic Institute, is
preserved at `third_party/wpi/LICENSE`.

## Separately obtained RTL

The WPI `rtl/obi/`, `rtl/ise/`, and coprocessor sources are upstream materials.
The archived root README describes Apache-2.0 licensing, while the artifact's
subdirectory license explicitly leaves the PQC RTL under the parent project's
terms. This release includes source identities and adaptation patches, and
obtains the original RTL from the user's upstream source tree at build time.
Upstream context in patches retains its original ownership and terms; the
root MIT license covers the TIL-SCA-authored additions.
The complete Apache-2.0 license text is included at
[`third_party/Apache-2.0.txt`](third_party/Apache-2.0.txt) for that upstream
context, following the parent README's license statement. The original
upstream source tree supplies any file-specific notices.

`provenance/rtl_sources.json` records normalized-text hashes of the evaluated
upstream files. The recorded archive base commit is
`46dbacf0975c85877db49885c2f13a8bb7767297`; the archive also records local edits,
so per-file checksums identify the actual inputs used here.

## Statistical methods

Welford accumulation, Welch's test, and moment-based higher-order TVLA are
established methods. The manuscript provides their references. TIL-SCA's
contribution is their integration with the RTL execution path and observation
identity, together with the reported boundary diagnosis and matched intervention.

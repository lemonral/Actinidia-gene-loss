# Frozen production phylogeny and gene-family toolchain

## Scope and decision

The production stack is frozen in
`config/phylogeny/toolchain.tsv`. It is deliberately separate from the
submitted-method reproduction stack. The revised workflow uses OrthoFinder
3.1.5, DIAMOND 2.1.10, MAFFT 7.526, IQ-TREE 2.4.0, the current
ASTRAL-Pro implementation distributed as ASTRAL-Pro3 1.25.3.8 in ASTER
v1.25, PAML/MCMCTree 4.10.10, CAFE 5.1.0, NCBI BLAST+ 2.16.0+, and R
4.4.0.

IQ-TREE 2.4.0 is intentionally the last pinned release in the IQ-TREE 2
series used by this protocol; the existence of IQ-TREE 3 does not silently
change the method. Likewise, `ASTRAL-Pro3` is the current C++ reimplementation
of the ASTRAL-Pro objective, not a switch to the single-copy ASTRAL method.
The official ASTER tutorial describes ASTRAL-Pro3 as faster and lower-memory
while retaining the ASTRAL-Pro objective and accepting multi-copy gene-family
trees.

This file freezes software only. It does not authorize a run. Cohort, assembly,
annotation, locus, calibration, and family-count gates in the other project
documents still apply.

## Final production audit status

The production executable audit is complete. The path-free V4 table and its
validation are curated in `results/phylogeny/toolchain_audit`; each accepted
analysis also has a run-specific executable/input checksum binding. Text below
that describes an executable as awaiting a formal audit records the earlier
discovery phase and should be read together with that final V4 artifact. CAFE5
is the sole special case: it has no version banner, so the accepted 5.1.0 build
is bound by its executable checksum and frozen source/build-file checksums
rather than by help text.

## Why this installation is not a single unconstrained Conda environment

Several programs publish self-contained binaries, while ASTRAL-Pro3 and CAFE
are best tied to an exact source commit or release archive. A single solver
transaction can change helper programs without changing the top-level package
request. Production therefore uses isolated, versioned prefixes outside Git
and an explicit PATH assembled for each run. Conda or mamba may be used to
provide compilers and libraries, but a solver environment is not accepted as
the software lock by itself.

The recommended external layout is:

```text
$DATA/tools/
  orthofinder/3.1.5/
  iqtree2/2.4.0/
  aster/v1.25-db2b3e95/
  paml/4.10.10/
  cafe5/5.1.0/
  mafft/7.526/
  diamond/2.1.10/
  ncbi-blast/2.16.0+/
  R/4.4.0/
```

Downloaded archives, compilers, binaries, and scientific outputs stay outside
Git. The repository keeps only versions, upstream URLs, checksums, commands,
and path-free run manifests.

## Fail-closed acquisition procedure

Perform acquisition in a new staging directory, never directly in a final
tool prefix.

1. Read the row for the tool from `config/phylogeny/toolchain.tsv`.
2. Download only the exact `artifact_url` over HTTPS. A transport proxy may be
   used, but proxy URLs, credentials, and node names must never be written to a
   project manifest or log.
3. Verify the archive before extraction. For rows with a hexadecimal
   `artifact_sha256`, require exact equality. For ASTER, require the exact Git
   commit `db2b3e95da5bb0318b933afe1a144eb943ef7cbf`; a branch name alone is
   insufficient.
4. For a future BLAST+ reacquisition, NCBI publishes the pinned MD5
   `48f66c9e01ea5136e381b2bf6fc62036`; verify it, calculate the downloaded
   archive SHA-256, and add that value to a dated acquisition lock before a new
   installation. The completed production searches instead bind all four
   already installed 2.16.0+ executables by exact V4 and per-run checksums; the
   acquisition placeholder must not be misread as an unbound completed run.
5. Extract into staging with no superuser privileges. Inspect archive member
   paths before extraction and reject absolute paths or `..` traversal.
6. Run the declared probe. `exact_version_banner` requires the pinned version
   token. `program_identity_smoke_only` proves only that the executable starts
   as the expected program; it cannot establish the pinned release and must be
   combined with archive/build provenance. A similarly named executable on
   PATH is not evidence. The probe must also return a code in the reviewed
   `allowed_probe_exit_codes_json` list. The default contract is `[0]`; any
   nonzero exception must be tool-specific and documented.
7. Calculate SHA-256 for every executable actually invoked and validate the
   file identity and hash both before and after its probe, then atomically
   rename staging to the final versioned prefix.
8. Write a path-free `software_versions.tsv` and full command log into the run
   directory. Include the tool ID, expected version, observed banner, release
   or commit, archive checksum, executable checksum, compiler and linked BLAS
   metadata where relevant, command, UTC time, and status.

For a downloaded archive, the minimum check is:

```bash
printf '%s  %s\n' EXPECTED_SHA256 ARCHIVE > expected.sha256
sha256sum --check expected.sha256
```

Do not pipe an unchecked archive directly into an extraction command.

## Tool-specific installation and verification

### OrthoFinder 3.1.5 and its helpers

Use the official release page and SHA-256-pinned Linux bundle, but install it
in an isolated Python 3.12 virtual environment as described by the official
OrthoFinder 3.1.5 README. An older completed run records OrthoFinder 3.0.1b1,
but the 2026-07-17 read-only inventory also found a live 3.1.5 executable. The
live executable is a reuse candidate only after its launcher **and installed
package code** hashes, release-archive binding, configuration, and helper
resolution have been recorded. The observed launcher is only a 251-byte
wrapper, so its hash cannot identify the installed package by itself. The old
3.0.1b1 run cannot satisfy this production pin.

OrthoFinder includes helper binaries, including DIAMOND. The production run
must not select an embedded helper by accident. Construct a minimal PATH in
which the pinned DIAMOND 2.1.10, MAFFT 7.526, and IQ-TREE 2.4.0 directories
precede the OrthoFinder bundle, retain a run-local copy of the release
`user_config.json`, and resolve each command with `command -v` before the run.
Record those resolved executable hashes. Use the explicit MSA workflow; the
OrthoFinder internal species tree is diagnostic, not the final dated topology.

## Display topology, species covariance tree, and optional unit QC

The primary publication topology contains one representative haploid
complement per nonhybrid biological species. The confirmed F1 hybrid
*A. x zhejiangensis* is the sole exception: parental haplomes A and B are
retained as separately labelled genome-unit tips so their two placements remain
visible. This mixed genome-unit display is not called a strictly
biological-species tree and is prohibited as input to dating, CAFE, or PGLS.

The primary parental-lineage analysis retains *A. x zhejiangensis* A and B as
separate terminal lineages and separate PGLS rows. Each receives its own
callable non-shared loss numerator and denominator. The project does not add a
shared-individual correlation constraint and does not require an A/B-pruned
sensitivity analysis.

Every accepted haplotype or subgenome still receives a complete, separate
JCVI, SynOrths, translated-genome search, callable-state, and gene-loss result.
Complete loss, partial/homeolog loss, and callable copy-opportunity loss are
additional aggregation summaries; they never replace the per-unit calls.

An all-unit tree is not a mandatory result. It is generated only when QC,
orthogroup occupancy, contamination checks, or an unexpected haplome placement
needs diagnosis. Representative-swap trees remain required sensitivities for
polyploid or phased species used in the biological-species analysis.

The production search/alignment/tree settings must be explicit, for example
`-S diamond -M msa -A mafft`, plus the reviewed OrthoFinder tree-method key
from the exact 3.1.5 configuration. Never assume that a default is unchanged
between OrthoFinder releases. Total project scientific workers remain at or
below fifteen, with server load and available memory checked before another
heavy process is launched.

Official sources: [OrthoFinder 3.1.5 release](https://github.com/OrthoFinder/OrthoFinder/releases/tag/v3.1.5)
and [OrthoFinder documentation](https://orthofinder.github.io/OrthoFinder/).

### IQ-TREE 2.4.0

Use `iqtree2 --version` and require version 2.4.0. Per-locus and concatenated
runs must record the complete ModelFinder, partition, support, seed, and thread
arguments. `-T AUTO` is not allowed under the project-wide ten-worker cap;
pass a bounded integer.

Official sources: [IQ-TREE 2.4.0 release](https://github.com/iqtree/iqtree2/releases/tag/v2.4.0)
and [IQ-TREE manual](https://www.iqtree.org/doc/).

### ASTRAL-Pro3 1.25.3.8

Checkout the exact ASTER commit and build only `make astral-pro`. The TAPER
submodule is not needed for that target and must not be fetched implicitly.
The built program prints `Version: v1.25.3.8` to stderr. Record the C++ compiler
version, CPU architecture, complete build command, and resulting executable
SHA-256. Provide the exact gene-copy-to-biological-species mapping and retain
the rooted output plus stderr. The output is a topology/branch-length analysis,
not a time tree.

Official sources: [ASTER v1.25 release](https://github.com/chaoszhang/ASTER/releases/tag/v1.25)
and [ASTRAL-Pro3 tutorial](https://github.com/chaoszhang/ASTER/blob/v1.25/tutorial/astral-pro3.md).

### PAML/MCMCTree 4.10.10

The release archive SHA-256 is published upstream. The official 4.10.10 binary
prints `MCMCTREE in paml version 4.10.10, 27 Jan 2026`, so its probe is an
`exact_version_banner`, not an identity-only smoke test. Run it in an empty
temporary directory so no stray `mcmctree.ctl` can trigger an analysis. In
that deliberately empty directory it exits 255 after the banner because the
control file is absent. Exit 255 is therefore the sole reviewed nonzero probe
code in this manifest; any other MCMCTree exit code fails. Preserve the release
README, archive binding, and executable hash even though the exact banner now
provides direct version evidence.

Every dating run needs an immutable control file, accepted rooted topology,
alignment checksum, calibration table checksum, random seed, burn-in/sample
settings, stdout/stderr, and convergence diagnostics. An executable passing a
smoke test does not activate a fossil calibration.

Official sources: [PAML 4.10.10 release](https://github.com/abacus-gene/paml/releases/tag/v4.10.10)
and [PAML repository](https://github.com/abacus-gene/paml).

### CAFE 5.1.0

The observed CAFE5 executable starts successfully and `--help` prints its
usage, but neither `--help`, `--version`, `-v`, nor `version` reports `5.1.0`.
The help probe is therefore declared `program_identity_smoke_only`; it must
never be reported as a version match. Reuse is possible only when the
executable SHA-256 is bound to the pinned source archive or documented build,
compiler, and BLAS/LAPACK linkage. Otherwise compile the checksum-pinned
upstream release. CAFE must receive the matching rooted, binary, ultrametric
MCMCTree tree and exactly one reviewed representative-complement count per
biological species.

Official source: [CAFE 5.1 release](https://github.com/hahnlab/CAFE5/releases/tag/v5.1).

The completed production bundle binds the executable, source/build metadata,
dated-tree checksum, and family-count checksum independently of that limited
banner. Base Poisson passed exact output closure with lambda
0.0085698614157905 and 15,066 analyzed families. Gamma3 failed initialization
and is frozen as unavailable; it is neither retried through further family
filtering nor presented as a supported result.

### MAFFT 7.526

The 2026-07-17 read-only inventory found MAFFT 7.526 and hashed its launcher.
Reuse it after the formal path-free audit records both the launcher and the
implementation executable it dispatches. If it is replaced, use the exact
upstream 7.526 package and manifest checksum. MAFFT core is BSD licensed;
the full package contains separately licensed extension code, so retain the
upstream license notices even though this project uses protein alignment only.

Official sources: [MAFFT 7.526 downloads](https://mafft.cbrc.jp/alignment/software/linux.html)
and [MAFFT license notices](https://mafft.cbrc.jp/alignment/software/license66.txt).

### DIAMOND 2.1.10

The read-only inventory found DIAMOND 2.1.10 and reported its executable hash.
Reuse is allowed only after that observation is captured in the formal
path-free audit and its archive provenance closes; otherwise install the
pinned binary. Its directory must precede the OrthoFinder bundle on PATH.

Official source: [DIAMOND 2.1.10 release](https://github.com/bbuchfink/diamond/releases/tag/v2.1.10).

### NCBI BLAST+ 2.16.0+

This pin covers the standalone suite used by SynOrths and the translated
gene-loss evidence stage. It is not substituted by DIAMOND. The 2026-07-17
inventory found `blastp`, `tblastn`, `tblastx`, and `makeblastdb` reporting
2.16.0+ and reported their executable hashes, but this does not yet prove that
the full suite came from the same archive. Require
`blastp -version`, `tblastn -version`, `tblastx -version`, and
`makeblastdb -version` to report the same 2.16.0+ package, and hash all four
executables. Do not install until the SHA-256 acquisition lock described above
is completed.

Official sources: [NCBI 2.16.0 release directory](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/),
[BLAST+ release notes](https://www.ncbi.nlm.nih.gov/books/NBK131777/), and
[NCBI public-domain statement](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html).

### R 4.4.0

The source release is SHA-256 pinned, and the 2026-07-17 inventory found R
4.4.0 and reported its launcher hash. Reuse still requires that hash in the
formal path-free audit plus build provenance. Base R is not a lock for an R
analysis: create and review a project `renv.lock` for every package used, record
`sessionInfo()`, and never allow a later package restore to change a completed
run. The primary PGLS implementation may remain the reviewed Python workflow;
R is retained for declared companion validation, plotting, or a separately
specified R PGLS sensitivity, not as an undocumented second implementation.

Official sources: [CRAN R 4 source archive](https://cran.r-project.org/src/base/R-4/)
and [R licensing](https://www.r-project.org/Licenses/).

## Server inventory interpretation

The probe and inventory columns in the manifest distinguish what a command can
actually prove from what still needs installation or provenance audit:

- `exact_version_observed_checksum_pending` means the version was found, but
  it is not production-ready until its executable hash and provenance pass;
- `different_version_observed` means a legacy executable or run exists but
  cannot satisfy the production pin;
- `present_version_unverified` means legacy use is evident but no acceptable
  version/hash record exists; and
- `program_identity_observed_version_provenance_pending` means the program
  starts but its own probe contains no exact version, so build/archive
  provenance remains the version lock; and
- `not_observed` means the local audit contains no evidence of an installation.

The `server_inventory_status` columns retain the deliberately conservative
initial-discovery state and are not rewritten after a run. The final 2026-07-17
V4 audit is curated at `results/phylogeny/toolchain_audit/toolchain_audit.tsv`
and is checksum-bound to the unchanged manifest. It records stable executable
hashes and exact-version matches for OrthoFinder 3.1.5, IQ-TREE 2.4.0,
ASTRAL-Pro3 1.25.3.8, MCMCTree 4.10.10, MAFFT 7.526, DIAMOND 2.1.10, all four
declared BLAST+ executables at 2.16.0+, and R 4.4.0. CAFE5 exposes program
identity but no version banner; its accepted 5.1.0 identity is therefore bound
separately by the executable checksum and frozen source/build file hashes in
`results/phylogeny/toolchain_audit/cafe5_build_provenance.json`. Production run
manifests further bind the exact executable and inputs used by each analysis.

## Read-only, path-free executable audit

`scripts/phylogeny/inventory_toolchain.py` closes the executable-hash and
declared-probe part of the gate without putting server paths into Git. It
accepts a path-bearing TSV registry that **must remain outside this
repository**:

```text
tool_id	executable_id	executable_path
orthofinder	orthofinder	/absolute/tool-prefix/bin/orthofinder
blast_plus	blastp	/absolute/tool-prefix/bin/blastp
blast_plus	tblastn	/absolute/tool-prefix/bin/tblastn
blast_plus	tblastx	/absolute/tool-prefix/bin/tblastx
blast_plus	makeblastdb	/absolute/tool-prefix/bin/makeblastdb
```

Run the audit outside any scientific analysis directory:

```bash
python3 scripts/phylogeny/inventory_toolchain.py \
  --registry /external/private/tool-registry.tsv \
  --output /external/run-manifests/toolchain-audit.tsv \
  --strict
```

The script executes only the manifest-declared probe arguments, without a
shell, in an empty temporary directory. Its child environment is rebuilt from
a narrow runtime allowlist; proxy, API-token, credential, user-profile, and
unrelated application variables are not inherited. `HOME` and temporary-file
locations point to the disposable directory, user Python/R startup files are
disabled, `TERM=dumb` is fixed for deterministic noninteractive output, and
numerical-library thread counts are fixed to one.

It hashes and records the resolved executable before and after the probe and
fails if the symlink target, file identity, mode, size, modification time, or
content changes. It also refuses a registry stored inside the repository and
refuses to overwrite an audit output. The output contains tool/executable IDs,
expected versions, declared evidence level, executable sizes and SHA-256
values, the exact allowed-exit-code JSON, observed probe exit, raw probe-output
byte counts/SHA-256 values, normalized-view SHA-256, stripped-SGR count, and a
restricted matched probe token. It contains neither resolved executable paths
nor raw banners.

Version matching never runs on unchecked terminal output. The raw bytes are
hashed first, decoded as strict UTF-8, and normalized only by converting CRLF
or CR line endings to LF and removing standard ANSI SGR color/style sequences
of the form `ESC [ parameters m`. Any other escape sequence, C0/C1 control,
Unicode control/format character, or invalid UTF-8 fails the audit. The exact
version regex itself is not broadened. This permits the colored OrthoFinder
3.1.5 banner to match after its SGR codes are removed while retaining the raw
28-byte observation and its checksum as evidence.

`PASS_EXACT_VERSION_MATCH` means both the declared exact-version pattern and a
reviewed exit code matched. An undeclared nonzero code fails even if the banner
contains the expected version; conversely, an allowed code cannot rescue a
wrong version. All tools use `[0]` except the empty-directory MCMCTree probe,
which uses `[255]` for its documented missing-control-file exit.
`PASS_PROGRAM_IDENTITY_SMOKE` means only the program-name smoke pattern
matched with a reviewed exit code; its `version_match` value is `not_tested`.
Thus a strict CAFE inventory can confirm that the executable starts without
falsely claiming the pinned release. Archive/build provenance remains
mandatory.

## Per-run preflight gate

Before OrthoFinder, any gene-tree job, dating, CAFE, or PGLS starts, require all
of the following:

1. no unresolved checksum placeholder for a tool used by that stage;
2. all resolved executables lie under reviewed versioned prefixes;
3. version banners and executable SHA-256 values match the run manifest;
4. the exact input cohort and every input checksum are frozen;
5. the aggregate requested scientific-worker count is at most 15;
6. commands use explicit output directories, seeds where supported, and
   captured stdout/stderr; and
7. the run cannot overwrite or resume a directory from another cohort or tool
   version.

A failed item stops the run. It is never converted into a warning merely to
keep a batch moving.

# Upstream Blockers

Workarounds in this repo that exist solely because of an upstream constraint
(broken release, hardcoded behavior, naming policy, missing feature). Each
entry should be cheap to revisit: it lists what to remove, where it lives,
and a one-line check that will start failing once the upstream catches up.

When you fix or revisit one of these, **delete the entry** — this file is
not a changelog. If a check turns out to be ambiguous or expensive, simplify
it before merging the entry; otherwise the file will rot.

## How to use

1. When you upgrade a pinned package or notice a related upstream release,
   skim this file and run the **Check** for any entry that touches what
   you changed.
2. If the check passes, lift the workaround and remove the entry.
3. If you add a new workaround, append an entry below following the same
   shape (Constraint / Workaround / Upstream / Check / Remove when).

---

## `setuptools<71` in `epibutton_chip.yaml`

- **Constraint:** `epic2` (broad-peak caller for CUT&RUN) imports
  `pkg_resources` at startup. setuptools 71 unbundled `pkg_resources` into
  a separate PyPI distribution (`pkg_resources` package), and conda-forge
  has not yet shipped it, so any setuptools >=71 in the env breaks `epic2`.
- **Workaround:** `setuptools<71` pin in `workflow/envs/epibutton_chip.yaml`.
- **Upstream:**
  - <https://github.com/pypa/setuptools/issues/4519> (deprecation/removal of
    bundled `pkg_resources`)
  - <https://github.com/biocore-ntnu/epic2> (no release that drops the
    `pkg_resources` import yet)
- **Check:** in a fresh chip env without the pin,
  `python -c "import epic2"` succeeds.
- **Remove when:** epic2 cuts a release that imports `importlib.resources`
  / `importlib.metadata` instead, OR conda-forge ships the standalone
  `pkg_resources` package.

## `numpy<=1.19` in `epibutton_idr.yaml`

- **Constraint:** the `idr` package (Irreproducible Discovery Rate) uses
  `np.int` / `np.float` / `np.bool` aliases that numpy 1.20 deprecated and
  1.24 removed. Without the pin, `idr` imports fail with
  `AttributeError: module 'numpy' has no attribute 'int'`.
- **Workaround:** `numpy<=1.19` pin in `workflow/envs/epibutton_idr.yaml`.
  The env is single-purpose (IDR + bedtools + meme) so the old numpy is
  only loaded for IDR-calling rules, not the rest of the pipeline.
- **Upstream:** <https://github.com/nboley/idr> — repo is dormant; the
  scalar-alias deprecations have been raised in issues but no release.
- **Check:** in a fresh idr env without the pin,
  `python -c "import idr"` succeeds.
- **Remove when:** a maintained fork of IDR ships a numpy-1.24-compatible
  release, or we replace `idr` with an alternative caller.

## DMRcaller `parallel=TRUE, BPPARAM=...` plumbing

- **Constraint:** DMRcaller >= 1.42 added a new `parallel=FALSE` default
  to `computeDMRs` that hardcodes `BPPARAM = SerialParam()` regardless
  of `cores=` or any globally-registered BiocParallel backend. The old
  `cores=` parameter is silently ignored on those versions, so a job that
  requests 16 CPUs runs single-threaded.
- **Workaround:** `workflow/scripts/R_call_DMRs.R` registers
  `MulticoreParam(workers=threads)` AND passes
  `parallel=TRUE, BPPARAM=MulticoreParam(workers=threads)` explicitly when
  those formals are present (introspected via `formals()` so the script
  still works on DMRcaller 1.38 and earlier, which used `parallel::mclapply`
  directly via `cores=`).
- **Upstream:**
  - GitHub issue tracking the regression in this repo:
    <https://github.com/cahnlab/EPICC/issues/23>
  - DMRcaller upstream:
    <https://bioconductor.org/packages/DMRcaller/> (no flag has been added
    to make registered BPPARAM the default again).
- **Check:** with DMRcaller pinned to a future release, run a `dmC` test
  with `--cores 8` and confirm `call_DMRs_pairwise` consumes >1 CPU
  (e.g. `top` shows ~700% during the rule, or `time` shows
  `user >> real`).
- **Remove when:** DMRcaller defaults `BPPARAM` to `bpparam()` again (i.e.
  honors the registered backend), at which point the explicit pair can
  collapse back to a one-line `register(MulticoreParam(...))` call.

## Per-genome `.libPaths()` for AnnotationForge GO packages

- **Constraint:** `AnnotationForge::makeOrgPackage()` enforces a strict
  package name format `org.<G><species>.eg.db` (one initial of the genus
  + species). For projects that use multiple reference genomes of the
  *same species* (e.g. ColCEN and TAIR10 — both `org.Athaliana.eg.db`)
  this collides in the conda env's shared R library: the second build
  silently overwrites the first.
- **Workaround:** `workflow/scripts/R_build_GO_database.R` and
  `workflow/scripts/R_GO_analysis.R` install/load each genome's GO
  package into its own per-genome library directory under
  `genomes/<refgenome>/GO/`, with `.libPaths(c(getwd(), .libPaths()))`
  prepended so `library(dbname, character.only=TRUE)` resolves to the
  right one.
- **Upstream:** AnnotationForge has not relaxed the dbname constraint;
  the request would be a `name=` argument to `makeOrgPackage()` that
  bypasses the auto-derived `org.<G><species>.eg.db`.
- **Check:** in a fresh dev env, build GO databases for two same-species
  genomes back-to-back without the per-genome lib dance and confirm both
  `library(org.Athaliana.eg.db)` invocations return the genome-specific
  GAF table.
- **Remove when:** `makeOrgPackage()` exposes a `package_name=` (or
  similar) argument that lets us namespace by reference genome, OR we
  switch off AnnotationForge to a tool with no name-collision policy.

## Bismark 3.x aligner `--nucleotide_coverage` not implemented

- **Constraint:** Bismark 3.x (the Rust rewrite of the suite) recognises the
  aligner's `--nucleotide_coverage` flag but no-ops it — the log prints
  "these options are recognised but not yet active in this build (wired in a
  later phase)" and no `*.nucleotide_stats.txt` is written. Passing the
  resulting (nonexistent) file to `bismark2report --nucleotide_report` then
  crashes it ("os error 2"). Perl Bismark ≤0.25 produced this file at align
  time.
- **Workaround:** `workflow/rules/mC.smk` drops `--nucleotide_coverage` from
  the `bismark` align call and instead runs a separate `bam2nuc` step on the
  deduplicated BAM (bam2nuc *is* implemented in 3.x). Output is named
  `{PE,SE}__<sample>.deduplicated.nucleotide_stats.txt` (after the BAM, not the
  align-time `..._bismark_bt2_pe.nucleotide_stats.txt`), declared as a rule
  output and fed to `bismark2report --nucleotide_report`.
- **Upstream:** <https://github.com/FelixKrueger/Bismark> (Rust suite; the flag
  is staged "for a later phase").
- **Check:** with a future Bismark, run `bismark ... --nucleotide_coverage`
  on a small mC sample and confirm it writes `*_bismark_bt2_pe.nucleotide_stats.txt`
  at align time (no "not yet active" note in the log).
- **Remove when:** the 3.x aligner implements `--nucleotide_coverage`, at which
  point the separate `bam2nuc` step can be dropped and the align-time flag +
  original nucleotide-report path restored.

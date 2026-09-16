# Reviewer reproduction boundary

This repository is a sanitized source derivative of staged release v3.5.2.
The frozen server release was not changed when this repository was prepared.
The per-file hashes in provenance.json distinguish byte-identical files from
files with replaced deployment paths and pilot identifiers.

## Available now

- Caller/QC source and command construction.
- Stage ordering, receipt validation, incremental exports and regression tests.
- Single-pass pilot recount source and legacy analysis modules for inspection.
- Synthetic test results for this exact source snapshot.

## Required for independent end-to-end reproduction

1. Archive the exact production container digest and make the base image or a
   fully resolved, tested rebuild recipe available under appropriate licenses.
2. Supply exact reference, capture BED and resource manifests with versions and
   hashes; public download locations alone do not identify the executed bytes.
3. Obtain controlled access to required study inputs or provide a separately
   validated synthetic/public end-to-end example.
4. Regenerate configuration and hash-bound receipts for the authorized inputs.
   Sanitized pilot identifiers intentionally do not match private receipts.
5. Freeze and validate the post-gate cohort recount scope implementation before
   reporting it as executed. Legacy global-union code is not the cohort plan.
6. Attach completed production receipts, calibration results and final analysis
   provenance when available, after reviewing them for disclosure risk.

The included shell helpers are study deployment templates. They require
untracked private inputs, including the analysis-plan documents referenced by
prepare_final_run.sh. Their absence is intentional. Do not disable hash checks
or replace absent evidence with synthetic success receipts to launch a study.

Synthetic test success demonstrates covered software behavior only. It does
not measure low-VAF sensitivity, validate biological conclusions or establish
that the cohort run has finished.

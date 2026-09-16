# Genomic Aging: staged WES SNV/indel pipeline

Research source accompanying the Genomic Aging in Li-Fraumeni Syndrome study.
This repository records the implementation derived from production release
**v3.5.2**, with a staged SAGE, Mutect2 and LoFreq workflow and separately
developed pileup evidence recount.

## Status and validation boundary

As of the latest operator report in September 2026, production QC for the
1,329-sample study was in progress. This is not a receipt of completed cohort
calling, somatic-origin classification, or biological validation.

The production release passed 101 synthetic/regression tests locally. A
three-sample pilot recount matched supported baseline evidence across 228,384
sample–variant pairs. Neither result establishes sensitivity or specificity
for the full study. See `validation/local_tests.json` for verification of this
sanitized repository snapshot.

## Workflow

1. **QC:** BAM validation, compact interval and chromosome-compartment coverage,
   global coverage histogram, GetPileupSummaries and CalculateContamination.
2. **SAGE:** primary discovery, with frozen `-skip_bqr -skip_msi_jitter`.
3. **Mutect2:** later sensitivity stage, reusing QC contamination evidence.
4. **LoFreq:** later sensitivity stage; serial caller implementation.
5. **CALLERS_COMPLETE:** issued after verification of all four stage receipts.

The `primary` invocation runs QC followed by SAGE and stops. Individual stages
can be invoked separately. Hash-bound receipts support resume and verify
predecessor outputs; an interrupted sample is not treated as complete.

The runnable staged image does not expose union/recount entry points. Candidate
union, calibration, frozen evaluation scopes and subsequent evidence recount
are a separate workflow. The pilot recount uses
`pysam.AlignmentFile.pileup`, including a single-coordinate-pass engine in
`aging/single_pass.py`; this is not the `samtools mpileup` command.
The full-cohort scope policy must be reconciled with the frozen analysis plan
and validated before that separate workflow is released. Presence of
`scopes.py` or legacy recount code does not imply production authorization.
Not-evaluated evidence must remain missing rather than become zero.

## Repository layout

- `aging/staged_production.py`: production stage entry point.
- `aging/runtime.py`: QC and caller commands.
- `aging/compact_qc.py`: coverage summaries, with zero-depth territory retained.
- `aging/stages.py`, `stage_contract.py`: orchestration and receipts.
- `aging/stage_export.py`: incremental export of completed, verified outputs.
- `aging/single_pass.py`, `window_recount.py`, `evidence.py`: pilot/legacy recount.
- `tests/`: synthetic fixtures and regression tests; no patient BAMs.
- `resources/compartments.GRCh38.json`: coverage compartment definitions.
- `provenance.json`: per-file upstream and repository hashes.

Additional historical analysis modules are retained for source audit and
regression testing; they are not stages of the current production image.

## Run the synthetic tests

Python 3.11 or 3.12 is recommended:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-test.txt
python -m unittest discover -s tests -q
```

Tests generate their own temporary data. They do not run the full production
cohort or replace validation with the actual caller binaries.

## Production reproduction

The supplied Dockerfile extends the study's internal
`genomic-aging-lfs:3.0-callers-permfix2` image. That image is **not distributed
here**, so the repository alone cannot yet rebuild the full caller environment.
Its toolchain includes SAGE 5.0.2, Mutect2/GATK, LoFreq, samtools and GRCh38
resources. Review reproduction requires the exact base image/digest and resource
manifests from the study team, subject to their respective licenses.

With the base image available, `bash build_image.sh` tests the source before
building the staged image. Production preparation additionally requires:

- authorized workbook and exact sample linkage;
- BAMs and indexes, reference FASTA and appropriate resource indexes;
- capture territory, population/contamination resources and SAGE resources;
- frozen base configuration, analysis plan and evidence receipts.

The shell launch/preparation helpers document the study deployment layout.
Paths have been replaced with `/srv/genomic-aging` and pilot IDs with
`PILOT_A/B/C`. They are templates requiring authorized local configuration,
not a turnkey demo. Private inputs and plans are intentionally absent. Receipt
hashes must be regenerated from authorized local evidence; the original study
receipts cannot be reused with the sanitized identifiers. See
`docs/REPRODUCIBILITY.md`.

## SAGE protocol and limitations

The observed SAGE binary required a sample-specific REDUX BQR table that was
absent from the study inputs. The study therefore explicitly selected
`skip_bqr` and `skip_msi_jitter`. The failed enabled-BQR benchmark does not
measure the scientific benefit or performance of BQR. REDUX base-quality
recalibration is not applied, and default jitter behavior is retained.
Sensitivity calibration, strand/orientation evidence and spectrum checks are
planned analytical safeguards, not substitutes for demonstrated validation.

Caller-native evidence and germline/PoN evidence must be interpreted with
callability and calibration. Tumor-only candidate calling does not by itself
establish somatic origin. The pipeline is research software.

## Data, citation and release

No workbook, patient name, date of birth, clinical linkage, BAM/VCF, private
results or operational logs are included. Example identifiers in source are
placeholders. Do not commit study data to this repository.

For review, cite the exact Git commit and `provenance.json`. A publication DOI,
archival code release and licensing decision remain to be added by the authors.
No third-party tools or reference datasets are redistributed here.

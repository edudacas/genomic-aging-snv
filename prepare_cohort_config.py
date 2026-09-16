"""Prepare the final frozen v3.5.2 config and protocol receipt; never launches."""
import argparse
import copy
import json
from pathlib import Path
from aging.bqr_decision import validate_skip_bqr_decision
from aging.common import sha, write_json
from aging.compact_qc import DEFAULT_THRESHOLDS
from aging.protocol import fingerprint as protocol_fingerprint


def prepare(base,decision,evidence_root,release,out,workers=12,cpus=16,memory_gib=96,
            source_base_sha256=None):
    c=copy.deepcopy(base);out=Path(out);release=Path(release)
    validate_skip_bqr_decision(decision,evidence_root)
    if c.get('build')!='GRCh38':raise ValueError('GRCh38 required')
    if c.get('threads')!=2 or c.get('java_heap_gb')!=4 or c.get('keep_derived_bam') is not False:
        raise ValueError('base caller resource protocol changed')
    if c.get('callers')!=['mutect','lofreq','sage']:
        raise ValueError('base caller list changed')
    if c.get('sage_extra_args')!=['-skip_bqr','-skip_msi_jitter']:
        raise ValueError('base SAGE flags are not the confirmed protocol')
    if not 1<=workers<=12 or not 1.25*workers<=cpus<=16 or not 4*workers<memory_gib<=96:
        raise ValueError('invalid final resource budget')
    for key in ('recount_scope','single_pass','evidence_engine','export_parquet','study_input'):
        c.pop(key,None)
    c.update(execution_mode='STAGED_CALLERS_ONLY',profile='exploratory',label_origins=False,
             sample_workers=workers,cpus=cpus,memory_gib=memory_gib,
             source_base_config_sha256=source_base_sha256,
             analysis_plan_path='/work/ANALYSIS_PLAN_GENOMIC_AGING_LFS.v3.1.md',
             analysis_plan_sha256=sha(release/'docs/ANALYSIS_PLAN_GENOMIC_AGING_LFS.v3.1.md'),
             compact_qc=dict(mask_path='/opt/aging/resources/compartments.GRCh38.json',
                 mask_sha256=sha(release/'resources/compartments.GRCh38.json'),thresholds=DEFAULT_THRESHOLDS),
             runtime_sha256=sha(release/'aging/runtime.py'),
             bqr_decision=dict(selected='skip_bqr',skip_msi_jitter=True,confirmed_by_user=True,
                 evidence_receipt='/work/bqr_decision.json',evidence_sha256=sha(out/'bqr_decision.json'),
                 copied_evidence_root='/work/bqr_evidence'),
             caller_protocol_receipt='/work/caller_protocol_unchanged.json')
    receipt=dict(status='CALLER_PROTOCOL_FROZEN',sage_variant='skip_bqr_and_skip_msi_jitter',
      sage_extra_args=['-skip_bqr','-skip_msi_jitter'],
      sage_other_arguments_unchanged_from_v3_4=True,
      mutect_arguments_unchanged_from_v3_4=True,
      mutect_contamination_reused_from_qc=True,
      lofreq_arguments_unchanged_from_v3_4=True,
      lofreq_parallelized=False,
      bqr_evidence_sha256=c['bqr_decision']['evidence_sha256'],
      analysis_plan_sha256=c['analysis_plan_sha256'],
      no_union=True,no_recount=True,protocol_projection_sha256=protocol_fingerprint(c))
    write_json(out/'caller_protocol_unchanged.json',receipt)
    c['caller_protocol_receipt_sha256']=sha(out/'caller_protocol_unchanged.json')
    write_json(out/'config.json',c)
    return c,receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('base','decision','evidence-root','release','out'):p.add_argument('--'+name,required=True)
    p.add_argument('--workers',type=int,default=12);p.add_argument('--cpus',type=int,default=16)
    p.add_argument('--memory-gib',type=int,default=96);a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    prepare(json.loads(Path(a.base).read_text()),json.loads(Path(a.decision).read_text()),
            a.evidence_root,a.release,out,a.workers,a.cpus,a.memory_gib,sha(a.base))
    print(json.dumps(dict(status='COHORT_CONFIG_PREPARED_NOT_LAUNCHED',config=str(out/'config.json'),
                          config_sha256=sha(out/'config.json'),protocol_receipt_sha256=sha(out/'caller_protocol_unchanged.json'))))


if __name__=='__main__':main()

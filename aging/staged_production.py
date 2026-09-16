"""FULL1329 stage-major entry point. No automatic union or recount."""
import argparse
import fcntl
from pathlib import Path
from . import runtime
from .common import load_config, sha, write_json, write_tsv
from .staged_inputs import fingerprint, freeze, selected
from .stages import STAGES, execute_stage, resource_budget
from .stage_contract import validate_outputs
from .bqr_decision import validate_skip_bqr_decision
from .protocol import fingerprint as protocol_fingerprint, verify_implementation
from .workbook_input import load_workbook_samples


def validate_config(c):
    workers=resource_budget(c)
    if c.get('execution_mode')!='STAGED_CALLERS_ONLY' or c.get('keep_derived_bam') is not False:
        raise ValueError('require staged callers-only with derived BAM disposal')
    if sha(c['analysis_plan_path'])!=c['analysis_plan_sha256']:
        raise ValueError('analysis plan changed')
    if len(c.get('source_base_config_sha256') or '')!=64:
        raise ValueError('source base config hash missing')
    if c.get('profile')!='exploratory' or c.get('label_origins') is not False:
        raise ValueError('discovery only; no somatic origin classification')
    if c.get('callers')!=['mutect','lofreq','sage']:
        raise ValueError('preserve all three caller protocols')
    if not c.get('compact_qc'):
        raise ValueError('compact coverage required')
    if sha(c['compact_qc']['mask_path'])!=c['compact_qc']['mask_sha256']:
        raise ValueError('coverage mask changed')
    for resource in ('gnomad','contamination_sites'):
        if resource not in c['resources']:
            raise ValueError('missing protocol resource '+resource)
    if c['resources']['sage_jar']['path']!='/opt/sage/sage.jar':
        raise ValueError('require frozen image SAGE JAR')
    if any(c.get('secondary',{}).get(k) for k in ('cnv','telseq')):
        raise ValueError('secondary modules are outside this run')
    decision=c.get('bqr_decision',{})
    if decision.get('selected')!='skip_bqr' or decision.get('skip_msi_jitter') is not True or decision.get('confirmed_by_user') is not True:
        raise ValueError('cohort config requires confirmed skip_bqr and skip_msi_jitter')
    if sha(decision['evidence_receipt'])!=decision['evidence_sha256']:
        raise ValueError('BQR decision receipt changed')
    validate_skip_bqr_decision(load_config(decision['evidence_receipt']),decision['copied_evidence_root'])
    if c.get('sage_extra_args')!=['-skip_bqr','-skip_msi_jitter']:
        raise ValueError('SAGE arguments disagree with BQR decision')
    protocol=load_config(c['caller_protocol_receipt'])
    if sha(c['caller_protocol_receipt'])!=c['caller_protocol_receipt_sha256']:
        raise ValueError('caller protocol receipt changed')
    verify_implementation(c)
    if protocol.get('protocol_projection_sha256')!=protocol_fingerprint(c) or protocol.get('sage_variant')!='skip_bqr_and_skip_msi_jitter':
        raise ValueError('caller protocol receipt/config projection mismatch')
    return workers


def stage_callback(stage,c,out):
    def process(row,d):
        row=dict(row,germline_vcf='',carrier_variant='')
        if stage=='qc':
            runtime.qc(row,c,d)
            runtime.contamination(row,c,d)
        else:
            options=dict(c,callers=[stage])
            runtime.call_sample(row,options,d,c['pon_by_assay'][row['assay_id']],
                                contamination_directory=out/'qc'/row['sample_id'])
            if stage=='sage':
                runtime.run(['bcftools','index','-f','-t',d/'sage.vcf.gz'],d/'commands.log')
        validate_outputs(d,stage)
    return process


def execute(a):
    c=load_config(a.config);workers=validate_config(c)
    rows,input_receipt=load_workbook_samples(a.workbook,a.linkage,expected=1329,check_files=True)
    rows=selected(rows)
    if len(rows)!=1329 or any(r['explore_status']!='SELECTED_NOT_SERVER_VALIDATED' for r in rows):
        raise ValueError('FULL1329 cannot silently exclude technically pending participants')
    for key in ('person_id','analysis_id'):
        if not all(r.get(key) for r in rows) or len({r[key] for r in rows})!=len(rows):
            raise ValueError('missing or duplicated '+key)
    for row in rows:
        if not c.get('pon_by_assay',{}).get(row['assay_id']):
            raise ValueError('assay PoN missing')
        row['germline_vcf_deferred_source']=row.get('germline_vcf','')
        row['germline_vcf']=''
        row['carrier_variant']=''
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    with (out/'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        resources=freeze(c,None)
        fp=fingerprint(rows,c,dict(resources=resources,workbook=input_receipt))
        marker=out/'RUN_INPUTS.json'
        if marker.exists():
            if load_config(marker)['fingerprint']!=fp:
                raise ValueError('frozen inputs changed: use a new run directory')
        else:
            write_json(marker,dict(fingerprint=fp,config_sha256=sha(a.config),
                       selected_samples=1329,pilot_cache_imported=False,
                       caller_protocol_unchanged=dict(mutect=True,lofreq=True,
                           sage_except_confirmed_flags=True,bqr_variant='skip_bqr',
                           msi_jitter_variant='skip_msi_jitter',
                           receipt_sha256=c['caller_protocol_receipt_sha256']),
                       workbook_input=input_receipt))
            write_json(out/'resources.lock.json',resources)
            write_tsv(out/'samples.private.tsv',rows)
        stages=('qc','sage') if a.stage=='primary' else (a.stage,)
        for stage in stages:
            execute_stage(out,stage,rows,fp,workers,stage_callback(stage,c,out))
    return 0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=STAGES+('primary',))
    for name in ('workbook','linkage','config','out'):
        parser.add_argument('--'+name,required=True)
    return execute(parser.parse_args())


if __name__=='__main__':
    raise SystemExit(main())

"""Validation contract for the confirmed skip-BQR/skip-jitter decision."""
from pathlib import Path
from .common import sha


def validate_skip_bqr_decision(receipt,evidence_root=None):
    expected={
      'status':'BQR_NOT_FEASIBLE_WITH_CURRENT_INPUTS',
      'reason':'REQUIRED_SAMPLE_SPECIFIC_REDUX_BQR_TABLES_ABSENT',
      'expected_bqr_table_filename':'PILOT_A.redux.bqr.tsv',
      'sage_readme_describes_bqr_as_internal':True,
      'failed_arm_is_not_scientific_evidence_against_bqr':True,
      'bqr_performance_or_call_effect_measured':False,
      'redux_jar_available_in_image':False,
      'redux_program_in_bam_header':False,
      'redux_resources_available':False,
      'user_decision':'skip_bqr_and_skip_msi_jitter',
      'decision_confirmed':True,
      'requires_user_decision':False,
    }
    for key,value in expected.items():
        if receipt.get(key)!=value:raise ValueError('invalid BQR decision evidence: '+key)
    note=receipt.get('sage_readme_bqr_note','')
    if '5.0.2' not in note or 'REDUX' not in note or 'binário' not in note:
        raise ValueError('BQR README/binary discrepancy note missing')
    for key in ('failure_sha256','log_sha256'):
        digest=receipt.get(key,'')
        if len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
            raise ValueError('invalid benchmark evidence hash: '+key)
    if evidence_root is not None:
        root=Path(evidence_root)
        for key,path in {'failure_sha256':root/'FAILURE.json','log_sha256':root/'commands.log'}.items():
            if not path.is_file() or sha(path)!=receipt[key]:raise ValueError('copied BQR evidence does not match '+key)
    if receipt.get('limitations')!=[
      'base qualities not recalibrated by REDUX BQR',
      'context-dependent error not modeled by the caller BQR score',
      'microsatellite indels use default jitter parameters']:
        raise ValueError('declared BQR limitations changed')
    if receipt.get('compensatory_controls')!=[
      'calibration battery','strand and orientation recount counts',
      'mutational-spectrum deviation']:
        raise ValueError('declared compensatory controls changed')
    return receipt

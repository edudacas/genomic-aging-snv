"""Finalize the observed server receipt without modifying benchmark artifacts."""
import argparse
import json
from pathlib import Path
from aging.bqr_decision import validate_skip_bqr_decision
from aging.common import sha, write_json


def finalize(root):
    root=Path(root).resolve();path=root/'bqr_feasibility.json'
    if not path.is_file():raise ValueError('bqr_feasibility.json missing')
    receipt=json.loads(path.read_text())
    required={
      'status':'BQR_NOT_FEASIBLE_WITH_CURRENT_INPUTS',
      'reason':'REQUIRED_SAMPLE_SPECIFIC_REDUX_BQR_TABLES_ABSENT',
      'failure_sha256':'5b942eee5f2692e47798fa794d8579b64463401436ec3bee0bee06e68d0ead3d',
      'log_sha256':'d95fd00ea8bfbb54ceada26dd818e4c5186a28cd8b6bb8c44c3647d442a1b561',
      'failed_arm_is_not_scientific_evidence_against_bqr':True,
    }
    for key,value in required.items():
        if receipt.get(key)!=value:raise ValueError('unexpected original feasibility field: '+key)
    failure=root/'enabled/PILOT_A/FAILURE.json'
    log=root/'enabled/PILOT_A/commands.log'
    if sha(failure)!=receipt['failure_sha256'] or sha(log)!=receipt['log_sha256']:
        raise ValueError('benchmark artifacts do not match preserved hashes')
    backup=root/'bqr_feasibility.pre_decision.json'
    if backup.exists():
        if json.loads(backup.read_text()).get('status')!=required['status']:
            raise ValueError('unexpected existing backup')
    else:write_json(backup,receipt)
    receipt.update(
      expected_bqr_table_filename='PILOT_A.redux.bqr.tsv',
      sage_readme_describes_bqr_as_internal=True,
      sage_readme_bqr_note=('A documentação do SAGE 5.0.2 descreve o BQR como interno, mas está '
        'desatualizada em relação ao binário observado, que consome a tabela específica produzida pelo REDUX.'),
      user_decision='skip_bqr_and_skip_msi_jitter',
      decision_confirmed=True,
      requires_user_decision=False,
      limitations=['base qualities not recalibrated by REDUX BQR',
                   'context-dependent error not modeled by the caller BQR score',
                   'microsatellite indels use default jitter parameters'],
      compensatory_controls=['calibration battery','strand and orientation recount counts',
                             'mutational-spectrum deviation'],
    )
    validate_skip_bqr_decision(receipt)
    write_json(path,receipt)
    print(json.dumps(dict(status='BQR_DECISION_RECORDED',file=str(path),sha256=sha(path))))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--benchmark',required=True)
    finalize(p.parse_args().benchmark)


if __name__=='__main__':main()

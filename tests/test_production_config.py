import json
import tempfile
import unittest
from pathlib import Path
from aging.common import sha
from aging.staged_production import validate_config
from aging.protocol import fingerprint as protocol_fingerprint


class ProductionConfigTest(unittest.TestCase):
    def test_explicit_capacity_and_protocol(self):
        with tempfile.TemporaryDirectory() as temp:
            mask=Path(temp)/'mask'; mask.write_text('{}')
            evidence=Path(temp)/'evidence';evidence.mkdir()
            (evidence/'FAILURE.json').write_text('failure');(evidence/'commands.log').write_text('log')
            report=Path(temp)/'decision.json'
            decision=dict(status='BQR_NOT_FEASIBLE_WITH_CURRENT_INPUTS',
                reason='REQUIRED_SAMPLE_SPECIFIC_REDUX_BQR_TABLES_ABSENT',
                expected_bqr_table_filename='PILOT_A.redux.bqr.tsv',
                sage_readme_describes_bqr_as_internal=True,
                sage_readme_bqr_note='SAGE 5.0.2 README desatualizado; binário consome tabela REDUX.',
                failed_arm_is_not_scientific_evidence_against_bqr=True,
                bqr_performance_or_call_effect_measured=False,redux_jar_available_in_image=False,
                redux_program_in_bam_header=False,redux_resources_available=False,
                user_decision='skip_bqr_and_skip_msi_jitter',decision_confirmed=True,
                requires_user_decision=False,
                failure_sha256=sha(evidence/'FAILURE.json'),
                log_sha256=sha(evidence/'commands.log'),limitations=[
                  'base qualities not recalibrated by REDUX BQR',
                  'context-dependent error not modeled by the caller BQR score',
                  'microsatellite indels use default jitter parameters'],compensatory_controls=[
                  'calibration battery','strand and orientation recount counts','mutational-spectrum deviation'])
            report.write_text(json.dumps(decision))
            runtime=Path(__file__).parents[1]/'aging/runtime.py'
            plan=Path(temp)/'plan.md';plan.write_text('plan')
            protocol=Path(temp)/'protocol.json'
            base=dict(execution_mode='STAGED_CALLERS_ONLY',profile='exploratory',label_origins=False,
                      sample_workers=12,cpus=16,memory_gib=96,
                      build='GRCh38',threads=2,java_heap_gb=4,keep_derived_bam=False,
                      callers=['mutect','lofreq','sage'],sage_extra_args=['-skip_bqr','-skip_msi_jitter'],
                      resources={'gnomad':{},'contamination_sites':{},'sage_jar':{'path':'/opt/sage/sage.jar'}},
                      compact_qc={'mask_path':str(mask),'mask_sha256':sha(mask)},
                      runtime_sha256=sha(runtime),
                      source_base_config_sha256='0'*64,
                      analysis_plan_path=str(plan),analysis_plan_sha256=sha(plan),
                      bqr_decision={'selected':'skip_bqr','skip_msi_jitter':True,'confirmed_by_user':True,
                                    'evidence_receipt':str(report),'evidence_sha256':sha(report),
                                    'copied_evidence_root':str(evidence)},
                      caller_protocol_receipt=str(protocol))
            base['caller_protocol_receipt_sha256']='pending'
            protocol.write_text(json.dumps(dict(sage_variant='skip_bqr_and_skip_msi_jitter',
                                                protocol_projection_sha256=protocol_fingerprint(base))))
            base['caller_protocol_receipt_sha256']=sha(protocol)
            self.assertEqual(validate_config(base),12)
            for changes in [dict(sample_workers=13),dict(cpus=14),dict(memory_gib=97),
                            dict(bqr_decision={}),dict(sage_extra_args=[]),dict(execution_mode='CALLERS_ONLY')]:
                with self.assertRaises(ValueError):validate_config(base|changes)

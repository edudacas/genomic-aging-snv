import json
import tempfile
import unittest
from pathlib import Path
from aging.common import sha
from aging.protocol import fingerprint as protocol_fingerprint
from prepare_cohort_config import prepare


class FinalConfigTests(unittest.TestCase):
    def test_final_config_and_receipt_are_bound_to_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);release=root/'release';out=root/'run';evidence=root/'evidence'
            (release/'resources').mkdir(parents=True);(release/'aging').mkdir()
            (release/'docs').mkdir()
            (release/'resources/compartments.GRCh38.json').write_text('{}')
            (release/'aging/runtime.py').write_text('runtime')
            (release/'docs/ANALYSIS_PLAN_GENOMIC_AGING_LFS.v3.1.md').write_text('plan')
            out.mkdir();evidence.mkdir()
            (evidence/'FAILURE.json').write_text('failure');(evidence/'commands.log').write_text('log')
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
            (out/'bqr_decision.json').write_text(json.dumps(decision))
            base=dict(build='GRCh38',threads=2,java_heap_gb=4,keep_derived_bam=False,
                      callers=['mutect','lofreq','sage'],sage_extra_args=['-skip_bqr','-skip_msi_jitter'])
            config,receipt=prepare(base,decision,evidence,release,out,source_base_sha256='0'*64)
            self.assertEqual(config['sample_workers'],12);self.assertEqual(config['cpus'],16)
            self.assertEqual(config['memory_gib'],96)
            self.assertEqual(config['source_base_config_sha256'],'0'*64)
            self.assertEqual(receipt['sage_variant'],'skip_bqr_and_skip_msi_jitter')
            self.assertFalse(receipt['lofreq_parallelized'])
            self.assertEqual(receipt['protocol_projection_sha256'],protocol_fingerprint(config))
            self.assertEqual(config['caller_protocol_receipt_sha256'],sha(out/'caller_protocol_unchanged.json'))
            with self.assertRaisesRegex(ValueError,'requires_user_decision'):
                prepare(base,decision|{'requires_user_decision':True},evidence,release,out)
            (evidence/'commands.log').write_text('changed')
            with self.assertRaises(ValueError):prepare(base,decision,evidence,release,out)


if __name__=='__main__':unittest.main()

"""Compare actual argument arrays against the preserved v3.4 runtime source."""
import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from aging import runtime
from aging.common import atomic

spec=importlib.util.spec_from_file_location('aging._protocol_baseline',
     Path(__file__).parent/'fixtures/protocol_baseline_runtime.py')
baseline=importlib.util.module_from_spec(spec);spec.loader.exec_module(baseline)


class ProtocolTests(unittest.TestCase):
    def capture(self,module,caller,root,external=False):
        directory=root/'calls';directory.mkdir(exist_ok=True)
        qc=root/'qc';qc.mkdir(exist_ok=True)
        for name in ('contamination.table','segments.table'):atomic(qc/name,'fixture')
        commands=[]
        def run(argv,log,stdout=None):
            commands.append([str(x) for x in argv])
        config=dict(profile='exploratory',threads=2,java_heap_gb=4,callers=[caller],
                    keep_derived_bam=True,sage_build='38',sage_ensembl_dir='/ensembl',
                    evidence={'max_depth':100000},sage_extra_args=['-skip_bqr','-skip_msi_jitter'],
                    resources={k:{'path':'/'+k} for k in ('fasta','common_bed','gnomad',
                               'contamination_sites','sage_jar','sage_hotspots','sage_high_confidence','sage_driver_panel')})
        with patch.object(module,'validate_bam',return_value='TUMOR'), \
             patch.object(module,'bam_for_tools',return_value='/input.bam'), \
             patch.object(module,'bed_intervals',return_value=[('chr1',0,10)]), \
             patch.object(module,'run',side_effect=run),patch.object(module,'normalize') as norm:
            kwargs={'contamination_directory':qc} if external else {}
            module.call_sample({},config,directory,'/pon',**kwargs)
            normalized=[list(map(str,call.args[:2])) for call in norm.call_args_list]
        args=(directory/'sage.arguments.txt').read_text() if caller=='sage' else None
        return commands,normalized,args

    def test_sage_lofreq_and_mutect_arguments_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for caller in ('sage','lofreq','mutect'):
                old=self.capture(baseline,caller,root)
                new=self.capture(runtime,caller,root)
                self.assertEqual(old,new,caller)

    def test_mutect_consumes_qc_without_reestimating(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            old,old_norm,_=self.capture(baseline,'mutect',root)
            new,new_norm,_=self.capture(runtime,'mutect',root,external=True)
            self.assertFalse(any('GetPileupSummaries' in cmd or 'CalculateContamination' in cmd for cmd in new))
            expected=[cmd for cmd in old if 'GetPileupSummaries' not in cmd and 'CalculateContamination' not in cmd]
            expected=[[arg.replace(str(root/'calls/contamination.table'),str(root/'qc/contamination.table'))
                           .replace(str(root/'calls/segments.table'),str(root/'qc/segments.table')) for arg in cmd]
                      for cmd in expected]
            self.assertEqual(expected,new);self.assertEqual(old_norm,new_norm)

    def test_qc_contamination_commands_match_original(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            old,_,_=self.capture(baseline,'mutect',root)
            expected=[cmd for cmd in old if 'GetPileupSummaries' in cmd or 'CalculateContamination' in cmd]
            commands=[]
            config=dict(java_heap_gb=4,resources={'contamination_sites':{'path':'/contamination_sites'}})
            with patch.object(runtime,'validate_bam',return_value='TUMOR'), \
                 patch.object(runtime,'bam_for_tools',return_value='/input.bam'), \
                 patch.object(runtime,'run',side_effect=lambda argv,log:commands.append(list(map(str,argv)))):
                runtime.contamination({},config,root/'calls')
            self.assertEqual(expected,commands)


if __name__=='__main__':unittest.main()

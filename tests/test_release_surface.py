import unittest
from pathlib import Path


class ReleaseSurfaceTests(unittest.TestCase):
    def test_builder_contains_files_required_by_release_surface_audit(self):
        dockerfile=(Path(__file__).parents[1]/'Dockerfile').read_text()
        builder=dockerfile.split('\nFROM genomic-aging-lfs:3.0-callers-permfix2\n',1)[0]
        for source in ('Dockerfile','finalize_bqr_feasibility.py','prepare_final_run.sh'):
            self.assertIn(f'COPY {source} /opt/aging/{source}',builder)

    def test_final_image_exposes_only_staged_entrypoint(self):
        dockerfile=(Path(__file__).parents[1]/'Dockerfile').read_text()
        final=dockerfile.rsplit('FROM genomic-aging-lfs:3.0-callers-permfix2',1)[1]
        self.assertIn('aging.staged_production',final)
        self.assertNotIn('aging.production',final)
        copied=final.split('RUN chmod',1)[0]
        for module in ('reuse_callers.py','single_pass.py','optimized.py','window_recount.py',
                       'callers.py','evidence.py','cli.py','explore.py'):
            self.assertNotIn('/opt/aging/aging/'+module+' ',copied)
        self.assertIn('test ! -e /opt/aging/aging/reuse_callers.py',final)
        self.assertIn('test ! -e /opt/aging/aging/single_pass.py',final)

    def test_final_release_has_no_bqr_relaunch_helper(self):
        root=Path(__file__).parents[1]
        self.assertFalse((root/'run_bqr_benchmark.sh').exists())
        self.assertTrue((root/'finalize_bqr_feasibility.py').exists())
        self.assertTrue((root/'prepare_final_run.sh').exists())


if __name__=='__main__':unittest.main()

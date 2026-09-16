"""Dedicated FULL1329 production entrypoint; no recount command is exposed."""
import argparse
import json
from pathlib import Path
from .common import sha
from .callers import execute


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('workbook','linkage','config','out'): p.add_argument('--'+name,required=True)
    a=p.parse_args(); c=json.loads(Path(a.config).read_text())
    if c.get('execution_mode')!='CALLERS_ONLY' or c.get('keep_derived_bam') is not False:
        raise ValueError('Production requires CALLERS_ONLY and disposal of derived LoFreq BAMs')
    if not c.get('compact_qc') or c.get('threads')!=2 or c.get('java_heap_gb')!=4:
        raise ValueError('Frozen caller resource profile / compact QC mismatch')
    if not 1<=c.get('sample_workers',0)<=4:
        raise ValueError('Choose 1 to 4 workers after capacity collection')
    if sha(c['compact_qc']['mask_path'])!=c['compact_qc']['mask_sha256']:
        raise ValueError('Coverage mask changed')
    a.command='callers-only'; a.samples=None; a.sample_ids=None
    return execute(a)


if __name__=='__main__': raise SystemExit(main())

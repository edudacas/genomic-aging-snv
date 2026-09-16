"""Incremental snapshots; all state and temporary files stay outside the run."""
import argparse
import fcntl
import hashlib
import io
import json
import os
import tarfile
import tempfile
import time
from pathlib import Path

from .common import safe_id, sha, write_json
from .stages import STAGES, stage_fingerprint, verify_task
from .stage_contract import validate_outputs


def export_stage(run, destination, stage, fingerprint):
    run, destination = Path(run).resolve(), Path(destination).resolve()
    if stage not in STAGES:
        raise ValueError('unknown stage')
    if destination == run or run in destination.parents:
        raise ValueError('exports must be outside the run directory')
    destination.mkdir(parents=True, exist_ok=True)
    with (destination / '.export.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        exported = {}
        for receipt_path in destination.glob('export-*.json'):
            receipt = json.loads(receipt_path.read_text())
            if receipt['config_fingerprint'] != fingerprint or receipt['run'] != str(run):
                raise ValueError('export destination belongs to a different frozen run')
            if sha(destination / receipt['archive']) != receipt['archive_sha256']:
                raise ValueError('prior export archive missing or corrupted')
            if receipt['stage'] == stage:
                exported.update(receipt['samples'])
        selected, files = {}, {}
        index=STAGES.index(stage)
        previous=sha(run/(STAGES[index-1].upper()+'_COMPLETE.json')) if index else None
        expected_fingerprint=stage_fingerprint(fingerprint,stage,previous)
        for marker in sorted((run / stage).glob('*/COMPLETE.json')):
            sid = safe_id(marker.parent.name)
            digest = sha(marker)
            if sid in exported:
                if exported[sid] != digest:
                    raise ValueError('previously exported sample changed')
                continue
            task = json.loads(marker.read_text())
            verify_task(marker.parent, expected_fingerprint)
            # Completion is the publication boundary; never include an in-progress directory.
            directories = [marker.parent]
            if stage != 'qc':
                directories.append(run / 'qc' / sid)
            for directory in directories:
                qc_marker = directory / 'COMPLETE.json'
                metadata = json.loads(qc_marker.read_text())
                expected=expected_fingerprint if directory==marker.parent else stage_fingerprint(fingerprint,'qc')
                verify_task(directory, expected)
                validate_outputs(directory,stage if directory==marker.parent else 'qc',metadata)
                for name, expected in metadata['outputs'].items():
                    path = directory / name
                    files[str(path.relative_to(run))] = expected
                files[str(qc_marker.relative_to(run))] = sha(qc_marker)
            selected[sid] = digest
        if not selected:
            return dict(status='NO_NEW_SAMPLES', stage=stage)
        manifest = dict(stage=stage, run=str(run), config_fingerprint=fingerprint,
                        samples=selected, files=files, status='EXPORT_COMPLETE')
        label = 'export-' + stage + '-' + str(time.time_ns())
        archive = destination / (label + '.tar.gz')
        fd, temporary = tempfile.mkstemp(dir=destination, prefix='.partial-')
        os.close(fd)
        try:
            with tarfile.open(temporary, 'w:gz') as tar:
                for name, expected in sorted(files.items()):
                    path = run / name
                    # Read the exact bytes to be archived and validate against the task receipt.
                    data = path.read_bytes()
                    if hashlib.sha256(data).hexdigest() != expected:
                        raise ValueError('source changed during export: ' + name)
                    item = tarfile.TarInfo(name);item.size=len(data);item.mode=0o600
                    tar.addfile(item, io.BytesIO(data))
                payloads = {'manifest.json': json.dumps(manifest, indent=2, sort_keys=True) + '\n'}
                checks = dict(files)
                checks['manifest.json'] = hashlib.sha256(payloads['manifest.json'].encode()).hexdigest()
                payloads['SHA256SUMS'] = ''.join(f'{digest}  {name}\n' for name,digest in sorted(checks.items()))
                for name, text in payloads.items():
                    data=text.encode();item=tarfile.TarInfo(name);item.size=len(data);item.mode=0o600
                    tar.addfile(item, io.BytesIO(data))
            os.replace(temporary, archive)
            # Only a committed receipt advances the incremental cursor. A crash before
            # this point leaves an orphan archive; retry includes the samples again.
            manifest.update(archive=archive.name, archive_sha256=sha(archive))
            write_json(destination / (label + '.json'), manifest)
            return manifest
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=STAGES)
    parser.add_argument('--run', required=True)
    parser.add_argument('--destination', required=True)
    args=parser.parse_args()
    inputs=json.loads((Path(args.run)/'RUN_INPUTS.json').read_text())
    receipt=export_stage(args.run,args.destination,args.stage,inputs['fingerprint'])
    print(json.dumps(dict(status=receipt['status'],stage=args.stage,
                         samples=len(receipt.get('samples',{})),archive=receipt.get('archive'))))


if __name__=='__main__':
    main()

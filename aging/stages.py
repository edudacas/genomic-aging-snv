"""Stage-major execution primitives. No discovery union or recount entry point.

The caller owns the exclusive run lock and the frozen input fingerprint.
Completed sample directories are immutable and safe for concurrent export.
"""
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .common import safe_id, sha, write_json

STAGES = ('qc', 'sage', 'mutect', 'lofreq')


def resource_budget(config):
    workers = config['sample_workers']
    cpus = config['cpus']
    memory = config['memory_gib']
    if type(workers) is not int or not 1 <= workers <= 12:
        raise ValueError('workers must be an integer from 1 to 12')
    if type(cpus) not in (int, float) or not workers * 1.25 <= cpus <= 16:
        raise ValueError('require 1.25 * workers <= cpus <= 16')
    if type(memory) not in (int, float) or not workers * 4 < memory <= 96:
        raise ValueError('memory must exceed aggregate Java heaps and be <= 96 GiB')
    if config['threads'] != 2 or config['java_heap_gb'] != 4:
        raise ValueError('preserve 2 threads and 4 GiB Java heap per sample')
    return workers


def verify_task(directory, fingerprint):
    directory = Path(directory)
    if directory.is_symlink():
        raise ValueError('completed task directory cannot be a symlink')
    receipt = json.loads((directory / 'COMPLETE.json').read_text())
    if receipt.get('status') != 'COMPLETE' or receipt.get('fingerprint') != fingerprint:
        raise ValueError('stale or incomplete sample task')
    if not receipt.get('outputs'):
        raise ValueError('empty completed sample task')
    for name, digest in receipt['outputs'].items():
        path = directory / name
        if Path(name).is_absolute() or '..' in Path(name).parts or path.is_symlink():
            raise ValueError('unsafe completed output')
        if directory.resolve() not in path.resolve().parents:
            raise ValueError('completed output escapes task directory')
        if not path.is_file() or sha(path) != digest:
            raise ValueError('completed output corrupted: ' + name)
    return receipt


def sample_task(directory, fingerprint, callback):
    directory = Path(directory)
    if (directory / 'COMPLETE.json').exists():
        return verify_task(directory, fingerprint), True
    if directory.exists():
        directory.rename(directory.with_name(directory.name + '.failed.' + str(time.time_ns())))
    directory.mkdir(parents=True)
    start = time.monotonic()
    try:
        callback(directory)
        outputs = {str(p.relative_to(directory)): sha(p)
                   for p in directory.rglob('*') if p.is_file() and not p.is_symlink()}
        if not outputs:
            raise ValueError('empty sample task')
        receipt = dict(status='COMPLETE', fingerprint=fingerprint, outputs=outputs,
                       elapsed_seconds=time.monotonic() - start)
        write_json(directory / 'COMPLETE.json', receipt)
        return receipt, False
    except Exception as exc:
        write_json(directory / 'FAILURE.json', dict(error=type(exc).__name__,
                   message=str(exc), elapsed_seconds=time.monotonic() - start))
        raise


def stage_fingerprint(fingerprint, stage, previous_digest=None):
    return hashlib.sha256(json.dumps([fingerprint, stage, previous_digest]).encode()).hexdigest()


def verify_stage(out, stage, fingerprint, sample_ids):
    """Verify the full dependency chain, including all completed sample outputs."""
    out = Path(out)
    index = STAGES.index(stage)
    previous = None
    if index:
        predecessor = STAGES[index - 1]
        verify_stage(out, predecessor, fingerprint, sample_ids)
        previous = sha(out / (predecessor.upper() + '_COMPLETE.json'))
    sfp = stage_fingerprint(fingerprint, stage, previous)
    receipt = json.loads((out / (stage.upper() + '_COMPLETE.json')).read_text())
    if (receipt.get('status') != stage.upper() + '_COMPLETE'
            or receipt.get('fingerprint') != sfp
            or sorted(receipt.get('samples', {})) != sorted(sample_ids)):
        raise ValueError('stage receipt does not match frozen run')
    for sid in sample_ids:
        path = out / stage / safe_id(sid)
        if receipt['samples'][sid]['receipt_sha256'] != sha(path / 'COMPLETE.json'):
            raise ValueError('sample receipt changed')
        verify_task(path, sfp)
    return receipt


def execute_stage(out, stage, rows, fingerprint, workers, callback):
    """Run just one stage; missing predecessors fail closed, never auto-execute."""
    if stage not in STAGES:
        raise ValueError('unknown caller stage')
    if type(workers) is not int or not 1 <= workers <= 12:
        raise ValueError('invalid workers')
    sample_ids = [safe_id(r['sample_id']) for r in rows]
    if not sample_ids or len(set(sample_ids)) != len(sample_ids):
        raise ValueError('empty or duplicate sample registry')
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    index = STAGES.index(stage)
    previous = None
    if index:
        predecessor = STAGES[index - 1]
        verify_stage(out, predecessor, fingerprint, sample_ids)
        previous = sha(out / (predecessor.upper() + '_COMPLETE.json'))
    marker = out / (stage.upper() + '_COMPLETE.json')
    if marker.exists():
        receipt = verify_stage(out, stage, fingerprint, sample_ids)
    else:
        sfp = stage_fingerprint(fingerprint, stage, previous)
        start = time.monotonic()
        completed, failures = {}, {}
        initial=dict(stage=stage,completed=0,failed=0,total=len(rows),status='RUNNING')
        write_json(out/'progress.json',initial)
        print(json.dumps(initial),flush=True)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            pending = {pool.submit(sample_task, out / stage / row['sample_id'], sfp,
                       lambda d, row=row: callback(row, d)): row['sample_id'] for row in rows}
            for future in as_completed(pending):
                sid = pending[future]
                try:
                    task, reused = future.result()
                    completed[sid] = dict(receipt_sha256=sha(out / stage / sid / 'COMPLETE.json'),
                                         elapsed_seconds=task['elapsed_seconds'], reused=reused)
                except Exception as exc:
                    failures[sid] = dict(error=type(exc).__name__, message=str(exc))
                progress=dict(stage=stage,completed=len(completed),failed=len(failures),total=len(rows),
                              sample_id=sid,sample_status='FAILED' if sid in failures else 'COMPLETE')
                write_json(out/'progress.json',progress)
                print(json.dumps(progress),flush=True)
        receipt = dict(status=stage.upper() + ('_FAILED' if failures else '_COMPLETE'),
                       config_fingerprint=fingerprint, fingerprint=sfp,
                       predecessor_receipt_sha256=previous, samples=completed, failures=failures,
                       attempt_wall_elapsed_seconds=time.monotonic() - start,
                       timing_scope='CURRENT_ATTEMPT; sample elapsed times are not additive wall time')
        if failures:
            write_json(out / (stage.upper() + '_FAILED.json'), receipt)
            write_json(out/'status.json',dict(status=stage.upper()+'_FAILED',failed_samples=sorted(failures)))
            raise RuntimeError(stage + ' failed for ' + str(len(failures)) + ' samples')
        write_json(marker, receipt)
    write_json(out/'status.json',dict(status=stage.upper()+'_COMPLETE',selected_samples=len(rows),
                                    union_performed=False,recount_performed=False))
    if stage == 'lofreq':
        # This also recovers a crash between LOFREQ_COMPLETE and CALLERS_COMPLETE.
        for name in STAGES:
            verify_stage(out, name, fingerprint, sample_ids)
        write_json(out / 'CALLERS_COMPLETE.json', dict(status='CALLERS_COMPLETE',
                   config_fingerprint=fingerprint, selected_samples=len(rows),
                   union_performed=False, recount_performed=False,
                   stages={s: sha(out / (s.upper() + '_COMPLETE.json')) for s in STAGES}))
        write_json(out/'status.json',dict(status='CALLERS_COMPLETE',selected_samples=len(rows),
                                         union_performed=False,recount_performed=False))
    return receipt

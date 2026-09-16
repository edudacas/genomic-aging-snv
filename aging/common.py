"""Private, atomic I/O and auditable subprocess execution."""
import csv
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

os.umask(0o077)

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def atomic(path, text):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix='.' + path.name)
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(text); f.flush(); os.fsync(f.fileno())
        os.replace(tmp, path)
        path.chmod(0o600)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)

def write_json(path, data):
    atomic(path, json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n')

def write_tsv(path, rows, fields=None):
    import io
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    if not fields: raise ValueError('explicit fields required for empty table')
    buf = io.StringIO()
    w = csv.DictWriter(buf, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    w.writeheader(); w.writerows(rows); atomic(path, buf.getvalue())

def run(argv, log, stdout=None):
    log = Path(log); log.parent.mkdir(parents=True, exist_ok=True)
    with log.open('a') as err:
        err.write(json.dumps([str(x) for x in argv]) + '\n'); err.flush()
        if stdout:
            with open(stdout, 'w') as out:
                subprocess.run(list(map(str, argv)), check=True, stdout=out, stderr=err)
        else:
            subprocess.run(list(map(str, argv)), check=True, stdout=err, stderr=err)

def unique_index(rows, key):
    out = {}
    for r in rows:
        if r.get(key): out.setdefault(r[key], []).append(r)
    return out

def number(value):
    import math
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError): return None

def load_config(path):
    return json.loads(Path(path).read_text())

def safe_id(s):
    import re
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,100}', s):
        raise ValueError('unsafe sample identifier')
    return s

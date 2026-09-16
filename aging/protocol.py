"""Small non-circular fingerprint for the frozen caller protocol."""
import hashlib
import json
from pathlib import Path
from .common import sha


def projection(config):
    return {
      'callers':config['callers'],
      'sage_extra_args':config['sage_extra_args'],
      'threads':config['threads'],
      'java_heap_gb':config['java_heap_gb'],
      'keep_derived_bam':config['keep_derived_bam'],
      'runtime_sha256':config['runtime_sha256'],
      'bqr_evidence_sha256':config['bqr_decision']['evidence_sha256'],
      'analysis_plan_sha256':config['analysis_plan_sha256'],
      'source_base_config_sha256':config.get('source_base_config_sha256'),
    }


def fingerprint(config):
    return hashlib.sha256(json.dumps(projection(config),sort_keys=True,separators=(',',':')).encode()).hexdigest()


def verify_implementation(config):
    runtime=Path(__file__).parent/'runtime.py'
    if sha(runtime)!=config.get('runtime_sha256'):
        raise ValueError('caller runtime implementation changed')

"""GitHub tests/test_fragment_pilot_release.py: discover the isolated review pilot."""
import sys,unittest
from pathlib import Path

def load_tests(loader,tests,pattern):
    root=Path(__file__).resolve().parents[1]/'experiments'/'fragment_recount_pilot'
    sys.path.insert(0,str(root))
    return loader.discover(str(root),pattern='test_*.py',top_level_dir=str(root))

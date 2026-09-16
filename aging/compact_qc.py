"""Loss-accounted depth summaries; source depth is discarded only by runtime.qc."""
import bisect
import collections
import csv
import gzip
import hashlib
import json
import os
import tempfile
from pathlib import Path
from .common import sha, write_json, write_tsv


COMPARTMENTS = ('AUTOSOMES', 'X_NONPAR', 'X_PAR', 'Y_NONPAR_NONXTR', 'Y_XTR', 'Y_PAR', 'OTHER')
DEFAULT_THRESHOLDS = [10, 20, 30, 50, 80, 100, 200]


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def split_territory(bed, mask_path):
    mask = json.loads(Path(mask_path).read_text())
    if mask.get('build') != 'GRCh38' or mask.get('coordinates') != '0-based-half-open':
        raise ValueError('Coverage mask build/coordinate mismatch')
    masks = collections.defaultdict(list)
    for m in mask['masks']:
        if m['compartment'] not in COMPARTMENTS or not 0 <= m['start'] < m['end']:
            raise ValueError('Invalid compartment mask')
        if (m['chrom'], m['compartment']) not in [('chrX','X_PAR'), ('chrY','Y_PAR'), ('chrY','Y_XTR')]:
            raise ValueError('Unexpected chromosome/compartment mask')
        masks[m['chrom']].append(m)
    for chrom, values in masks.items():
        values.sort(key=lambda m: m['start'])
        if any(a['end'] > b['start'] for a,b in zip(values, values[1:])):
            raise ValueError('Overlapping masks')
    parents = collections.defaultdict(list)
    with open(bed) as stream:
        for line in stream:
            if not line.strip() or line.startswith(('#','track','browser')):
                continue
            c,s,e = line.split()[:3]; s,e = int(s),int(e)
            if not 0 <= s < e:
                raise ValueError('Invalid BED interval')
            parents[c].append((s,e))
    segments=[]
    for chrom, intervals in parents.items():
        intervals.sort()
        if any(a[1] > b[0] for a,b in zip(intervals, intervals[1:])):
            raise ValueError('Common BED must be nonoverlapping')
        if chrom in ('X','Y'):
            raise ValueError('Mask requires chrX/chrY contig names')
        default = ('AUTOSOMES' if chrom in {'chr'+str(n) for n in range(1,23)} else
                   'X_NONPAR' if chrom=='chrX' else 'Y_NONPAR_NONXTR' if chrom=='chrY' else 'OTHER')
        for start,end in intervals:
            borders = sorted({start,end} | {p for m in masks[chrom] for p in (m['start'],m['end']) if start < p < end})
            for s,e in zip(borders,borders[1:]):
                label=next((m['compartment'] for m in masks[chrom] if m['start'] <= s and e <= m['end']), default)
                segments.append(dict(chrom=chrom,start=s,end=e,parent_start=start,parent_end=end,
                                     compartment=label,territory_bases=e-s))
    if not segments:
        raise ValueError('Empty common BED')
    return segments, mask


def summarize(depth, bed, mask_path, out, filters, thresholds=None):
    """Stream a file or text iterator of samtools depth; omitted positions are zero."""
    thresholds = DEFAULT_THRESHOLDS if thresholds is None else thresholds
    if (not thresholds or any(type(x) is not int or x < 1 for x in thresholds)
            or thresholds != sorted(set(thresholds))):
        raise ValueError('Thresholds must be sorted unique positive integers')
    segments, mask = split_territory(bed, mask_path)
    starts=collections.defaultdict(list); indices=collections.defaultdict(list)
    for i,r in enumerate(segments):
        starts[r['chrom']].append(r['start']); indices[r['chrom']].append(i)
    sums=[0]*len(segments); covered=[0]*len(segments)
    buckets=[[0]*(len(thresholds)+1) for _ in segments]
    hist=collections.Counter(); last={}; source_hash=hashlib.sha256(); count=0
    own_stream=isinstance(depth,(str,Path)); stream=open(depth) if own_stream else depth
    try:
        for line in stream:
            source_hash.update(line.encode())
            fields=line.rstrip().split('\t')
            if len(fields)!=3:
                raise ValueError('Expected one BAM / three depth columns')
            c,p,d=fields; p=int(p)-1; d=int(d)
            if d<0 or p<=last.get(c,-1):
                raise ValueError('Negative depth or duplicate/unsorted position')
            last[c]=p
            j=bisect.bisect_right(starts.get(c,()),p)-1
            if j<0:
                raise ValueError('Depth outside BED')
            i=indices[c][j]
            if p>=segments[i]['end']:
                raise ValueError('Depth outside BED')
            sums[i]+=d; covered[i]+=int(d>0); hist[d]+=1; count+=1
            buckets[i][bisect.bisect_right(thresholds,d)]+=1
    finally:
        if own_stream: stream.close()
    total=sum(r['territory_bases'] for r in segments)
    if count>total: raise ValueError('Depth exceeds territory')
    hist[0]+=total-count
    measures=['territory_bases','depth_sum','covered_bases']+['bases_ge_'+str(t) for t in thresholds]
    totals={c:dict.fromkeys(measures,0) for c in COMPARTMENTS}
    for i,r in enumerate(segments):
        r.update(depth_sum=sums[i],covered_bases=covered[i])
        tail=0
        for k in range(len(thresholds)-1,-1,-1):
            tail+=buckets[i][k+1]; r['bases_ge_'+str(thresholds[k])]=tail
        for k in measures: totals[r['compartment']][k]+=r[k]
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    fields=['chrom','start','end','parent_start','parent_end','compartment']+measures
    dest=out/'coverage.intervals.tsv.gz'
    fd,tmp=tempfile.mkstemp(dir=out,prefix='.coverage.'); os.close(fd)
    try:
        with gzip.open(tmp,'wt',newline='') as f:
            w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(segments)
        # Reconstruct totals independently from the persisted interval summary.
        restored={c:dict.fromkeys(measures,0) for c in COMPARTMENTS}
        with gzip.open(tmp,'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                for k in measures: restored[r['compartment']][k]+=int(r[k])
        if restored!=totals: raise ValueError('Interval summary roundtrip failed')
        os.replace(tmp,dest)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)
    if sum(v['depth_sum'] for v in totals.values()) != sum(d*n for d,n in hist.items()):
        raise ValueError('Histogram/depth sum disagreement')
    for t in thresholds:
        if sum(v['bases_ge_'+str(t)] for v in totals.values())!=sum(n for d,n in hist.items() if d>=t):
            raise ValueError('Histogram/threshold disagreement')
    write_tsv(out/'depth_histogram.tsv',[dict(depth=d,bases=n) for d,n in sorted(hist.items())])
    write_tsv(out/'coverage.compartments.tsv',[dict(compartment=c,**totals[c]) for c in COMPARTMENTS])
    means={c:(v['depth_sum']/v['territory_bases'] if v['territory_bases'] else None) for c,v in totals.items()}
    auto=means['AUTOSOMES']
    ratio=lambda c: means[c]/auto if means[c] is not None and auto else None
    receipt=dict(status='COMPACT_COVERAGE_PASS',schema_version=1,coordinates='0-based-half-open',
        source_depth_coordinates='1-based',thresholds=thresholds,territory_bases=total,
        depth_sum=sum(d*n for d,n in hist.items()),covered_bases=total-hist[0],
        source_depth_sha256=source_hash.hexdigest(),bed_sha256=sha(bed),mask_sha256=sha(mask_path),
        mask_version=mask['version'],mask_xtr_status=mask.get('xtr_status'),
        filters=filters,filters_sha256=digest_json(filters),
        x_auto_ratio=ratio('X_NONPAR'),y_auto_ratio=ratio('Y_NONPAR_NONXTR'),
        sex_inferred='NOT_ASSESSED_GATE_REQUIRED',sex_classification_performed=False,
        interval_count=len(segments),outputs={p:sha(out/p) for p in
            ('coverage.intervals.tsv.gz','coverage.compartments.tsv','depth_histogram.tsv')})
    write_json(out/'coverage.summary.json',receipt)
    return hist,receipt


def discard_verified_depth(depth, directory):
    """Called after successful QC, before the task COMPLETE manifest is created."""
    directory=Path(directory); receipt=json.loads((directory/'coverage.summary.json').read_text())
    if receipt['status']!='COMPACT_COVERAGE_PASS' or sha(depth)!=receipt['source_depth_sha256']:
        raise ValueError('Depth source changed; retained')
    for name,digest in receipt['outputs'].items():
        if sha(directory/name)!=digest: raise ValueError('Compact output changed; depth retained')
    if not (directory/'qc.json').is_file(): raise ValueError('QC not complete; depth retained')
    Path(depth).unlink()

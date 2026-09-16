"""Fail-closed resource locks and scientific gate; no default calibration is a PASS."""
import json
import re
from bisect import bisect_right
from pathlib import Path
from .common import *

def territory_index(path):
    index={}
    for chrom,start,end in bed_intervals(path):
        starts,ends=index.setdefault(chrom,([],[]))
        starts.append(start);ends.append(end)
    return index

def within_territory(index,chrom,pos,ref):
    starts,ends=index.get(chrom,([],[]))
    i=bisect_right(starts,int(pos)-1)-1
    return i>=0 and int(pos)-1+len(ref)<=ends[i]

def bed_intervals(path, lengths=None):
    last = {}; seen = set(); current = None
    total = 0; rows = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(('#','track','browser')): continue
        c,s,e,*_ = line.split(); s,e = int(s),int(e)
        if c != current:
            if c in seen: raise ValueError('BED contigs must be contiguous')
            seen.add(c); current = c
        if s < 0 or e <= s or s < last.get(c,0): raise ValueError('BED invalid/overlapping/unsorted')
        if lengths and (c not in lengths or e > lengths[c]): raise ValueError('BED/reference mismatch')
        last[c] = e; total += e-s; rows.append((c,s,e))
    if not total: raise ValueError('empty BED')
    return rows

def freeze_resources(config, output):
    resources = config['resources']; locked = {}
    required = ['fasta','fai','dict','common_bed','driver_bed','driver_whitelist',
                'gnomad','gnomad_index','abraom','abraom_index','contamination_sites',
                'contamination_sites_index','sage_jar','sage_hotspots','sage_high_confidence',
                'sage_driver_panel','sage_ensembl_manifest']
    for name in required:
        entry = resources[name]; path = Path(entry['path']).resolve(strict=True)
        if not path.is_file() or path.stat().st_size == 0: raise ValueError('missing resource '+name)
        digest = sha(path)
        if not re.fullmatch('[0-9a-f]{64}',entry['sha256']) or digest != entry['sha256']:
            raise ValueError('resource checksum mismatch '+name)
        if entry['build'] != config['build']: raise ValueError('resource build mismatch '+name)
        locked[name] = dict(path=str(path), sha256=digest)
    # SAGE cache tree must be a supplied, checksummed list, not a mutable directory.
    cache = Path(config['sage_ensembl_dir']).resolve(strict=True)
    cache_rows=read_tsv(resources['sage_ensembl_manifest']['path'])
    if not cache_rows:raise ValueError('empty SAGE cache manifest')
    for row in cache_rows:
        p = (cache/row['relative_path']).resolve(strict=True)
        if not p.is_relative_to(cache) or sha(p) != row['sha256']: raise ValueError('SAGE cache mismatch')
    fasta=Path(resources['fasta']['path'])
    if Path(str(fasta)+'.fai').resolve()!=Path(resources['fai']['path']).resolve():raise ValueError('FAI must be adjacent to FASTA')
    if fasta.with_suffix('.dict').resolve()!=Path(resources['dict']['path']).resolve():raise ValueError('dictionary must be adjacent to FASTA')
    for key in ['gnomad','abraom','contamination_sites']:
        index=Path(resources[key+'_index']['path']).resolve()
        if index not in [Path(resources[key]['path']+ext).resolve() for ext in ['.tbi','.csi']]:
            raise ValueError('VCF index must be adjacent '+key)
    lengths = {f[0]:int(f[1]) for f in (l.split('\t') for l in Path(resources['fai']['path']).read_text().splitlines())}
    common = bed_intervals(resources['common_bed']['path'],lengths)
    driver = bed_intervals(resources['driver_bed']['path'],lengths)
    bychrom = {}
    for c,s,e in common: bychrom.setdefault(c,[]).append((s,e))
    for c,s,e in driver:
        if not any(a <= s and e <= b for a,b in bychrom.get(c,[])):
            raise ValueError('driver BED must be wholly inside common capture territory')
    for assay,entry in config['assays'].items():
        p=Path(entry['capture_bed']).resolve(strict=True)
        if sha(p)!=entry['capture_sha256']: raise ValueError('capture checksum '+assay)
        cap={}
        for c,s,e in bed_intervals(p,lengths):cap.setdefault(c,[]).append((s,e))
        for c,s,e in common:
            if not any(a<=s and e<=b for a,b in cap.get(c,[])):raise ValueError('common BED outside assay '+assay)
        locked['capture:'+assay]={'path':str(p),'sha256':sha(p)}
    write_json(output,dict(build=config['build'],resources=locked,config_sha256=None))
    return locked

def validate_gate(gate, samples_path, resource_lock, config_path):
    g = json.loads(Path(gate).read_text())
    if g.get('status') != 'FROZEN' or not g.get('reviewer') or not g.get('signed_at_utc'):
        raise ValueError('calibration gate not frozen and signed')
    for k,p in [('samples_sha256',samples_path),('resource_lock_sha256',resource_lock),('config_sha256',config_path)]:
        if g.get(k)!=sha(p): raise ValueError('gate input changed: '+k)
    for k in ['battery','power','statistical_plan','pilot_qc']:
        entry=g['evidence'][k]
        if sha(entry['path'])!=entry['sha256']: raise ValueError('gate evidence changed '+k)
    battery=json.loads(Path(g['evidence']['battery']['path']).read_text())
    if battery['status']!='CALIBRATION_REVIEW_READY':raise ValueError('calibration incomplete')
    if set(g['assays'])-set(battery['assays']):raise ValueError('assay absent from calibration battery')
    for entry in battery['evidence']:
        if sha(entry['path'])!=entry['sha256']:raise ValueError('nested calibration evidence changed')
    power=json.loads(Path(g['evidence']['power']['path']).read_text())
    if power.get('status')!='PILOT_POWER_REVIEW_REQUIRED' or not power.get('grid'):
        raise ValueError('pilot-based power evidence missing')
    if power.get('scope')!='FINAL_DESIGN_WITH_EXPECTED_CALLABILITY':
        raise ValueError('gate requires power on the final design, not pilot sample size')
    qc=json.loads(Path(g['evidence']['pilot_qc']['path']).read_text())
    if qc.get('status')!='PILOT_QC_REVIEWED' or not qc.get('reviewer') or not qc.get('umi_decision'):
        raise ValueError('pilot QC/UMI review missing')
    if g.get('orthogonal_status') not in ['VALIDATED','MODEL_BASED_NO_ORTHOGONAL']:
        raise ValueError('orthogonal validation declaration missing')
    if g.get('primary_population') != 'TREATMENT_NAIVE':
        raise ValueError('v1 implements treatment-naive primary only; change plan before fallback')
    for assay,v in g['assays'].items():
        f=float(v['exome_vaf_floor']); branch='A' if f<=0.03 else 'B'
        if v['branch']!=branch or not 0<f<0.5: raise ValueError('invalid endpoint switching rule')
        if branch=='B' and not 0<float(v['driver_vaf_floor'])<0.5:raise ValueError('driver floor missing')
        if int(v['min_alt'])<1 or not 0<float(v['sensitivity'])<1:raise ValueError('invalid callability rule')
        if g['orthogonal_status']=='MODEL_BASED_NO_ORTHOGONAL' and not v.get('conservative_floor_justification'):
            raise ValueError('conservative model-based floor requires justification')
    branches={v['branch'] for v in g['assays'].values()}
    if len(branches)!=1:raise ValueError('mixed assay branches cannot define a single primary outcome')
    if not 1<=g['tier']['min_callers']<=3:raise ValueError('invalid caller concordance')
    for k in ['min_mapq','min_baseq','min_alt_each_strand','min_alt_each_orientation','max_population_af','germline_alpha','min_pass_callers']:
        if k not in g['tier']: raise ValueError('tier incomplete '+k)
    return g

def selected_samples(path, stage, pilot_ids=None):
    rows=read_tsv(path)
    for r in rows:safe_id(r['sample_id'])
    if len({r['sample_id'] for r in rows})!=len(rows):raise ValueError('duplicate sample IDs')
    if stage=='qc': return [r for r in rows if r['technical_status']!='HOLD']
    if stage=='pilot':
        ids=set(Path(pilot_ids).read_text().split())
        if not ids or len(ids)>40:raise ValueError('pilot requires 1..40 explicit sample IDs')
        chosen=[r for r in rows if r['sample_id'] in ids]
        if len(chosen)!=len(ids):raise ValueError('pilot ID absent')
    else:
        chosen=[r for r in rows if r['analysis_status']=='METADATA_READY']
        if any(r.get('origin','COHORT')!='COHORT' for r in chosen):raise ValueError('derived calibration sample in biological cohort')
    if not chosen:raise ValueError('no eligible samples')
    if any(r['analysis_status']!='METADATA_READY' for r in chosen):raise ValueError('pilot contains HOLDs')
    if len({r['person_id'] for r in chosen})!=len(chosen):raise ValueError('duplicate people')
    return chosen

"""Input selection, hashing and resource freezing for the staged-only runner."""
import hashlib
import json
from pathlib import Path
from .common import read_tsv, safe_id, sha, write_json
from .contracts import bed_intervals


def fingerprint(rows,c,extra=None):
    code={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}
    files={}
    for row in rows:
        for key in ('bam_path','bai_path','germline_vcf'):
            if row.get(key):
                path=Path(row[key]);stat=path.stat()
                files[str(path.resolve())]=[stat.st_size,stat.st_mtime_ns,
                    row.get('bam_sha256','') if key=='bam_path' else sha(path)]
    return hashlib.sha256(json.dumps(dict(rows=rows,config=c,code=code,files=files,extra=extra),
                                        sort_keys=True).encode()).hexdigest()


def selected(rows,ids=None):
    seen=set();requested=set(ids) if ids is not None else None
    for row in rows:
        sid=safe_id(row['sample_id'])
        if sid in seen:raise ValueError('duplicate sample ID')
        seen.add(sid)
    if requested is not None and requested-seen:raise ValueError('unknown requested sample ID')
    result=[]
    for original in rows:
        row=dict(original)
        reasons=set(filter(None,row.get('technical_reasons','').split(';')))-{'STRATUM_MISSING'}
        if requested is not None and row['sample_id'] not in requested:state='NOT_REQUESTED'
        elif reasons or not row.get('bam_path') or not row.get('bai_path'):state='TECHNICAL_INPUT_PENDING'
        elif row.get('technical_status')=='HOLD' and not row.get('technical_reasons'):state='TECHNICAL_HOLD_UNSPECIFIED'
        else:state='SELECTED_NOT_SERVER_VALIDATED'
        row['explore_status']=state;result.append(row)
    return result


def freeze(c,output=None):
    resources=c['resources'];callers=c.get('callers',['mutect','lofreq','sage'])
    if callers!=['mutect','lofreq','sage']:
        raise ValueError('frozen staged protocol requires mutect, lofreq, sage')
    required={'fasta','fai','dict','common_bed','sage_jar','sage_hotspots',
              'sage_high_confidence','sage_driver_panel','sage_ensembl_manifest'}
    if required-set(resources):raise ValueError('missing technical resources')
    locked={}
    for name,entry in resources.items():
        path=Path(entry['path']).resolve(strict=True)
        if not path.is_file() or not path.stat().st_size:raise ValueError('empty resource '+name)
        digest=sha(path)
        if entry.get('sha256') and entry['sha256']!=digest:raise ValueError('checksum mismatch '+name)
        if entry.get('build',c['build'])!=c['build']:raise ValueError('reference build mismatch '+name)
        locked[name]={'path':str(path),'sha256':digest}
    fasta=Path(resources['fasta']['path'])
    if Path(str(fasta)+'.fai').resolve()!=Path(resources['fai']['path']).resolve():raise ValueError('FAI path mismatch')
    if fasta.with_suffix('.dict').resolve()!=Path(resources['dict']['path']).resolve():raise ValueError('dictionary path mismatch')
    lengths={x[0]:int(x[1]) for x in (line.split('\t') for line in Path(resources['fai']['path']).read_text().splitlines())}
    bed_intervals(resources['common_bed']['path'],lengths)
    for name in ('gnomad','abraom','contamination_sites'):
        if name in resources:
            index=resources.get(name+'_index',{}).get('path','')
            if not index or Path(index).resolve() not in [Path(resources[name]['path']+suffix).resolve() for suffix in ('.tbi','.csi')]:
                raise ValueError('indexed optional resource required: '+name)
    cache=Path(c['sage_ensembl_dir']).resolve(strict=True)
    entries=read_tsv(resources['sage_ensembl_manifest']['path'])
    if not entries:raise ValueError('empty SAGE cache manifest')
    for row in entries:
        path=(cache/row['relative_path']).resolve(strict=True)
        if not path.is_relative_to(cache) or sha(path)!=row['sha256']:raise ValueError('SAGE cache mismatch')
    for assay,path in c.get('pon_by_assay',{}).items():
        locked['pon:'+assay]={'path':path,'sha256':sha(path)}
        indexes=[Path(path+suffix) for suffix in ('.tbi','.csi') if Path(path+suffix).is_file()]
        if not indexes:raise ValueError('PoN index missing')
        for index in indexes:locked[str(index)]={'sha256':sha(index)}
    if output:write_json(output,{'build':c['build'],'resources':locked,'profile':'exploratory'})
    return locked

"""Prepare exactly three pilot candidate-union inputs; never launch BAM recount."""
import argparse, bisect, csv, gzip, hashlib, json, os, re
from collections import Counter,defaultdict
from pathlib import Path

CALLERS=('sage','mutect','lofreq')
STAGES=('qc',)+CALLERS

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def read(p):return json.loads(Path(p).read_text())

def sfp(fp,stage,previous):
    return hashlib.sha256(json.dumps([fp,stage,previous]).encode()).hexdigest()

def bed_index(path):
    raw=defaultdict(list)
    with Path(path).open() as f:
        for line in f:
            if not line.strip() or line.startswith(('#','track','browser')):continue
            c,s,e=line.split()[:3];s=int(s);e=int(e)
            if s<0 or e<=s:raise ValueError('Invalid BED coordinates')
            raw[c].append((s,e))
    merged={}
    for c,spans in raw.items():
        out=[]
        for s,e in sorted(spans):
            if out and s<=out[-1][1]:out[-1]=(out[-1][0],max(e,out[-1][1]))
            else:out.append((s,e))
        merged[c]=([s for s,e in out],[e for s,e in out])
    if not merged:raise ValueError('Empty common territory')
    return merged

def inside(bed,chrom,pos):
    if chrom not in bed:return False
    starts,ends=bed[chrom];idx=bisect.bisect_right(starts,pos-1)-1
    return idx>=0 and pos-1<ends[idx]

def parse_snv(line,bed):
    fields=line.rstrip('\r\n').split('\t')
    if len(fields)<8:raise ValueError('Malformed VCF')
    c,pos,_,ref,alts=fields[:5];pos=int(pos)
    if pos<1:raise ValueError('Invalid VCF POS')
    if not inside(bed,c,pos):return [],'outside_common_bed'
    rows=[]
    for alt in alts.split(','):
        if len(ref)==1 and len(alt)==1 and ref in 'ACGT' and alt in 'ACGT' and ref!=alt:
            rows.append((c,pos,ref,alt))
    return rows, 'non_snv_alleles' if len(rows)!=len(alts.split(',')) else None

def validate(stages,pilots,expected):
    with (stages/'samples.private.tsv').open(encoding='utf-8-sig',newline='') as f:
        registry=list(csv.DictReader(f,delimiter='\t'))
    ids={r['sample_id'] for r in registry}
    if any(not re.fullmatch(r'[A-Za-z0-9_-]+',sid) for sid in ids):raise ValueError('Unsafe sample ID')
    if len(registry)!=expected or len(ids)!=expected:raise ValueError('Invalid study registry')
    if len(pilots)!=3 or len(set(pilots))!=3 or not set(pilots)<=ids:raise ValueError('Require three distinct registered pilot IDs')
    chain=read(stages/'CALLERS_COMPLETE.json');fp=chain['config_fingerprint']
    if chain.get('status')!='CALLERS_COMPLETE' or chain.get('selected_samples')!=expected:raise ValueError('CALLERS_COMPLETE required')
    if chain.get('union_performed') is not False or chain.get('recount_performed') is not False:raise ValueError('Unexpected run boundary')
    previous=None;inputs={};stage_hashes={}
    for stage in STAGES:
        p=stages/(stage.upper()+'_COMPLETE.json');r=read(p);digest=sha(p)
        fingerprint=sfp(fp,stage,previous)
        if chain['stages'][stage]!=digest or r.get('fingerprint')!=fingerprint or r.get('config_fingerprint')!=fp:
            raise ValueError('Stage fingerprint mismatch')
        if r.get('status')!=stage.upper()+'_COMPLETE' or set(r.get('samples',{}))!=ids or r.get('failures'):
            raise ValueError('Incomplete stage')
        for sid in pilots:
            directory=stages/stage/sid;marker=directory/'COMPLETE.json';task=read(marker)
            if sha(marker)!=r['samples'][sid]['receipt_sha256'] or task.get('fingerprint')!=fingerprint or task.get('status')!='COMPLETE':
                raise ValueError('Pilot task receipt mismatch')
            if stage in CALLERS:
                filename=stage+'.norm.vcf.gz';path=directory/filename
                expected_hash=task.get('outputs',{}).get(filename)
                if not expected_hash or path.is_symlink() or sha(path)!=expected_hash:raise ValueError('Caller VCF hash mismatch')
                inputs[(sid,stage)]=(path,expected_hash)
        stage_hashes[stage]=digest;previous=digest
    return registry,fp,stage_hashes,inputs

def main():
    os.umask(0o077)
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stages',type=Path,required=True);ap.add_argument('--bed',type=Path,required=True)
    ap.add_argument('--pilot-sample',action='append',required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--expected',type=int,default=1329);a=ap.parse_args()
    stages=a.stages.resolve();out=a.out.resolve()
    if out==stages or stages in out.parents:raise ValueError('Output must be outside immutable production stages')
    if out.exists():raise ValueError('Use a new output directory')
    reg,fp,hashes,inputs=validate(stages,a.pilot_sample,a.expected)
    bed=bed_index(a.bed);out.mkdir(parents=True)
    union={};skips=Counter();native=0
    with gzip.open(out/'caller_records.PRIVATE.tsv.gz','wt',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','caller','variant_id','filter','native_vcf_record'])
        for (sid,caller),(path,digest) in sorted(inputs.items()):
            with gzip.open(path,'rt') as vcf:
                for line in vcf:
                    if line.startswith('#'):continue
                    native+=1;alleles,reason=parse_snv(line,bed)
                    if reason:skips[caller+'_'+reason]+=1
                    for allele in alleles:
                        key='{}:{}:{}:{}'.format(*allele)
                        item=union.setdefault(allele,{'samples':set(),'callers':set()})
                        item['samples'].add(sid);item['callers'].add(caller)
                        w.writerow([sid,caller,key,line.split('\t')[6],line.rstrip('\r\n')])
    with (out/'candidate_union.pilot.vcf').open('w') as f, (out/'site_scopes.PRIVATE.tsv').open('w',newline='') as s:
        f.write('##fileformat=VCFv4.2\n##source=GenomicAgingPilotPreparation\n##reference=GRCh38_UNVERIFIED_UNTIL_FASTA_CHECK\n')
        f.write('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Discovery caller names, not a positivity rule">\n')
        f.write('##INFO=<ID=STATE,Number=1,Type=String,Description="Evidence evaluation state">\n')
        for c in sorted({x[0] for x in union}):f.write('##contig=<ID='+c+'>\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        w=csv.writer(s,delimiter='\t');w.writerow(['variant_id','discovery_sample_ids','planned_pilot_sample_ids','evaluation_state'])
        def order(a):
            c=a[0][3:] if a[0].startswith('chr') else a[0];return (int(c) if c.isdigit() else 100,c,a[1],a[2],a[3])
        for allele,item in sorted(union.items(),key=lambda x:order(x[0])):
            c,pos,ref,alt=allele;key='{}:{}:{}:{}'.format(*allele)
            f.write(f'{c}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tCALLERS={",".join(sorted(item["callers"]))};STATE=RECOUNT_PENDING\n')
            w.writerow([key,';'.join(sorted(item['samples'])),';'.join(sorted(a.pilot_sample)),'RECOUNT_PENDING'])
    with (out/'pilot_samples.PRIVATE.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(reg[0]),delimiter='\t');w.writeheader();w.writerows(r for r in reg if r['sample_id'] in a.pilot_sample)
    files=['candidate_union.pilot.vcf','site_scopes.PRIVATE.tsv','pilot_samples.PRIVATE.tsv','caller_records.PRIVATE.tsv.gz']
    receipt=dict(status='PILOT_CANDIDATE_INPUTS_PREPARED_NOT_RECOUNTED',pilot_samples=a.pilot_sample,
        candidate_snv_alleles=len(union),planned_sample_allele_pairs=3*len(union),native_records_scanned=native,
        skip_counters=dict(skips),config_fingerprint=fp,stage_receipt_sha256=hashes,
        bed_sha256=sha(a.bed),script_sha256=sha(__file__),registry_sha256=sha(stages/'samples.private.tsv'),
        source_vcf_sha256={sid+'/'+caller:digest for (sid,caller),(path,digest) in inputs.items()},
        scope_gate_frozen=False,cohort_launch_authorized=False,recount_performed=False,callers_rerun=False,
        reference_ref_validation='PENDING',population_pon_germline='NOT_EVALUATED',
        note='Unclassified SNV union, all caller FILTER/VAF states; this is not the full cohort union',
        outputs={name:sha(out/name) for name in files})
    (out/'PREPARATION_RECEIPT.PRIVATE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))

if __name__=='__main__':main()

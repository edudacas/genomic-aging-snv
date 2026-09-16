"""Truthset is authoritative; registry is physical linkage, never case ascertainment."""
import collections
import re
from pathlib import Path
from .common import *

CLINICAL = ['cohort_person_id', 'age_at_collection', 'tissue', 'sex', 'treatment_status',
            'cancer_history', 'PC1', 'PC2', 'PC3', 'clinical_evidence', 'reviewer']

def load_truthset(path):
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True, data_only=False)
    if wb.sheetnames != ['truthset_wes_v2']: raise ValueError('unexpected workbook sheets')
    values=list(wb.active.values);wb.close()
    it = iter(values); headers = next(it)
    if len(headers) != len(set(headers)): raise ValueError('duplicate columns')
    rows = []
    for vals in it:
        if any(isinstance(v, str) and v.startswith('=') for v in vals):
            raise ValueError('formulas prohibited in frozen truthset')
        rows.append(dict(zip(headers, ['' if v is None else str(v).strip() for v in vals])))
    wb.close()
    for k in ['cohort_person_id', 'franklin_analysis_id']:
        vs = [r[k] for r in rows]
        if any(not v for v in vs) or len(vs) != len(set(vs)):
            raise ValueError('truthset requires unique nonempty ' + k)
    return rows

def prepare(truthset, registry, out, clinical=None, germline_index=None, assay_aliases=None):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    truth = load_truthset(truthset)
    regpath = Path(registry).resolve(strict=True)
    reg = unique_index(read_tsv(regpath), 'analysis_id')
    aliases={}
    if assay_aliases:
        for x in read_tsv(assay_aliases):
            if x['alias'] in aliases:raise ValueError('duplicate assay alias')
            aliases[x['alias']]=x['assay_id']
    clin = unique_index(read_tsv(clinical), 'cohort_person_id') if clinical else {}
    vcfs = read_tsv(germline_index) if germline_index else []
    vi = unique_index(vcfs, 'analysis_id_franklin')
    # FAP fallback is exact 12 digits, and must resolve to a single canonical record.
    vf = unique_index(vcfs, 'FAP')
    rows = []
    for t in truth:
        holds = []; tech = []
        sid = 'GA_' + hashlib.sha256(t['cohort_person_id'].encode()).hexdigest()[:16]
        matches = reg.get(t['franklin_analysis_id'], [])
        r = matches[0] if len(matches) == 1 else {}
        if len(matches) != 1: tech.append('REGISTRY_ABSENT' if not matches else 'REGISTRY_AMBIGUOUS')
        for a, b in [('franklin_run_id','batch_id')]:
            if t.get(a) and r.get(b) and t[a] != r[b]: tech.append('CONFLICT_' + b.upper())
        assay=aliases.get(t.get('assay_id'),t.get('assay_id')) or r.get('assay_id','')
        if assay and r.get('assay_id') and assay!=r['assay_id']:tech.append('CONFLICT_ASSAY_ID')
        if re.fullmatch(r'\d{12}',r.get('fap','')) and t.get('fap') and t['fap']!=r['fap']:
            tech.append('CONFLICT_FAP')
        exact_fap = t.get('fap','')
        if not re.fullmatch(r'\d{12}', exact_fap): holds.append('NO_EXACT_FAP')
        cr = clin.get(t['cohort_person_id'], [])
        if len(cr) > 1: raise ValueError('duplicate clinical person')
        c = cr[0] if cr else {}
        if number(c.get('age_at_collection')) is None or not 0 <= float(c['age_at_collection']) <= 120:
            holds.append('AGE_AT_COLLECTION_UNCONFIRMED')
        if c.get('tissue') != 'PERIPHERAL_BLOOD': holds.append('TISSUE_UNCONFIRMED')
        if c.get('sex') not in ['F', 'M']: holds.append('SEX_UNCONFIRMED')
        if c.get('treatment_status') not in ['NAIVE', 'EXPOSED']: holds.append('TREATMENT_UNKNOWN')
        if not c.get('clinical_evidence') or not c.get('reviewer'): holds.append('CLINICAL_PROVENANCE_MISSING')
        if c.get('cancer_history') not in ['YES','NO']: holds.append('CANCER_HISTORY_UNKNOWN')
        if any(number(c.get(k)) is None for k in ['PC1','PC2','PC3']): holds.append('ANCESTRY_MISSING')
        expected = 'TP53_NEGATIVE_BY_CURATION' if t['study_group'] == 'TP53_negative' else 'CONFIRMED_GERMLINE'
        if t['genetic_eligibility_status'] != expected: holds.append('GENETIC_PENDING')
        if t['batch_matched_design_member'] != 'YES': holds.append('OUTSIDE_FROZEN_DESIGN')
        if t['study_group'] != 'TP53_negative' and not t.get('variant_id'): holds.append('CARRIER_VARIANT_MISSING')
        vm = vi.get(t['franklin_analysis_id'], [])
        if not vm and exact_fap and re.fullmatch(r'\d{12}',exact_fap): vm = vf.get(exact_fap, [])
        vm = [v for v in vm if v.get('is_canonical_caller','').lower() in ['1','true','yes']]
        v = vm[0] if len(vm) == 1 else {}
        if len(vm) != 1: holds.append('GERMLINE_VCF_ABSENT_OR_AMBIGUOUS')
        batch = t['franklin_run_id']
        if not assay or not batch: tech.append('STRATUM_MISSING')
        row = dict(sample_id=sid, origin='COHORT', person_id=t['cohort_person_id'], analysis_id=t['franklin_analysis_id'],
            fap=exact_fap, group=t['study_group'], assay_id=assay, batch_id=batch,
            source_assay_id=t.get('assay_id',''),registry_fap_raw=r.get('fap',''),
            stratum=assay+'::'+batch, bam_path=r.get('bam_path',''), bai_path=r.get('bai_path',''),
            germline_vcf=v.get('vcf_path_host',''), germline_sample='',
            carrier_variant=t.get('variant_id',''), age_source_value=t['age_years'],
            source_sex=t.get('sex',''),source_sex_provenance=t.get('sex_source',''),
            genetic_status=t['genetic_eligibility_status'],
            bam_sha256=r.get('content_sha256_if_available',''),
            technical_status='HOLD' if tech else 'LINKED_NOT_SERVER_VALIDATED',
            technical_reasons=';'.join(tech), analysis_status='HOLD' if holds or tech else 'METADATA_READY',
            hold_reasons=';'.join(holds + tech))
        row.update({k:c.get(k,'') for k in CLINICAL if k != 'cohort_person_id'})
        rows.append(row)
    write_tsv(out/'samples.private.tsv', rows)
    template=out/'clinical_enrichment.template.private.tsv'
    if not template.exists():
        write_tsv(template,[dict(cohort_person_id=t['cohort_person_id']) for t in truth], CLINICAL)
    write_tsv(out/'download_recovery.private.tsv',
              [r for r in rows if r['technical_status']=='HOLD'], list(rows[0]))
    # No implicit prefix or subject ID resolution; crosswalk remains private.
    summary = dict(status='PREPARED_NOT_SCANNED', truthset_rows=len(rows),
        groups=dict(collections.Counter(r['group'] for r in rows)),
        linked=sum(r['technical_status']!='HOLD' for r in rows),
        metadata_ready=sum(r['analysis_status']=='METADATA_READY' for r in rows),
        hold_reasons=dict(collections.Counter(x for r in rows for x in r['hold_reasons'].split(';') if x)),
        truthset_sha256=sha(truthset), registry_sha256=sha(regpath), registry_resolved_path=str(regpath),
        clinical_sha256=sha(clinical) if clinical else None,
        germline_index_sha256=sha(germline_index) if germline_index else None,
        assay_aliases_sha256=sha(assay_aliases) if assay_aliases else None,
        server_accessed=False, real_bams_scanned=0)
    write_json(out/'cohort.audit.json',summary)
    return summary

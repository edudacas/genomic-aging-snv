"""Read the user-reviewed XLSX as the authoritative study selection and BAM input."""
import argparse
import hashlib
import json
import posixpath
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from .common import read_tsv, sha, write_json, write_tsv

NS = {'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}

def review_rows(path):
    # Standard library reader: no additional Excel package required in Docker.
    with zipfile.ZipFile(path) as archive:
        shared=[]
        if 'xl/sharedStrings.xml' in archive.namelist():
            root=ET.fromstring(archive.read('xl/sharedStrings.xml'))
            shared=[''.join(t.text or '' for t in si.findall('.//s:t',NS)) for si in root]
        sheets=ET.fromstring(archive.read('xl/workbook.xml')).findall('s:sheets/s:sheet',NS)
        sheet=next((s for s in sheets if s.attrib['name']=='Review'),None)
        if sheet is None:raise ValueError('Review sheet missing')
        rid=sheet.attrib['{'+NS['r']+'}id']
        rels=ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        target=next(r.attrib['Target'] for r in rels if r.attrib['Id']==rid)
        target=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
        def cells(row):
            result={}
            for cell in row.findall('s:c',NS):
                if cell.find('s:f',NS) is not None:raise ValueError('Formula in Review input; use explicit values')
                col=''.join(c for c in cell.attrib['r'] if c.isalpha())
                kind=cell.get('t'); value=cell.find('s:v',NS)
                if kind=='inlineStr':text=''.join(t.text or '' for t in cell.findall('.//s:t',NS))
                elif kind=='s':text=shared[int(value.text)] if value is not None else ''
                else:text=value.text or '' if value is not None else ''
                result[col]=text
            return result
        root=ET.fromstring(archive.read(target))
        allrows=root.findall('s:sheetData/s:row',NS)
        if not allrows:raise ValueError('Empty Review sheet')
        header=cells(allrows[0])
        if len(set(header.values()))!=len(header):raise ValueError('Duplicate workbook column headers')
        required={'individual_id','analysis_id','fap','BAM_path','assay_id','study_group','study_design_member'}
        if required-set(header.values()):raise ValueError('Required workbook columns missing')
        return [{name:data.get(col,'') for col,name in header.items()}
                for row in allrows[1:] if any((data:=cells(row)).values())]

def load_workbook_samples(workbook, linkage, expected=1329, check_files=False):
    rows=review_rows(workbook); links=read_tsv(linkage)
    by_analysis={r['analysis_id']:r for r in links}
    if len(by_analysis)!=len(links):raise ValueError('Duplicate linkage analysis IDs')
    if len(rows)!=expected:raise ValueError(f'Expected {expected} study rows; found {len(rows)}')
    for field in ('individual_id','analysis_id','fap','BAM_path'):
        values=[r[field].strip() for r in rows]
        if not all(values) or len(set(values))!=len(values):raise ValueError('Missing/duplicate '+field)
    if {r['analysis_id'] for r in rows}!=set(by_analysis):raise ValueError('Workbook membership differs from linked study')
    samples=[]
    for row in rows:
        link=by_analysis[row['analysis_id']]
        if any(row[f].strip()!=link[f] for f in ('individual_id','fap')):raise ValueError('Workbook identity differs from exact linkage')
        if row['study_design_member'].strip().lower()!='yes':raise ValueError('Study membership must be explicit for every row')
        path=row['BAM_path'].strip()
        if not Path(path).is_absolute() or path!=link['bam_path']:raise ValueError('Workbook BAM_path differs from linked manifest')
        if row['assay_id'].strip() not in ('',link['workbook_assay_id'],link['assay_id']):raise ValueError('Workbook assay changed; reconcile metadata')
        if check_files:
            for filename in (path,link['bai_path']):
                if not Path(filename).is_file() or not Path(filename).stat().st_size:raise ValueError('BAM/BAI absent or empty: '+filename)
            if Path(path).stat().st_size!=int(link['bam_bytes']):raise ValueError('BAM size differs from source manifest')
        samples.append(dict(sample_id=link['sample_id'],
                            person_id=link['person_id'],analysis_id=row['analysis_id'],fap=row['fap'],
                            bam_path=path,bai_path=link['bai_path'],assay_id=link['assay_id'],
                            batch_id=row.get('batch_id',''),study_group=row['study_group'],
                            genetic_eligibility_status=row.get('genetic_eligibility_status',''),
                            age_years=row.get('age_years',''),sex=row.get('sex',''),
                            technical_status='INPUT_LINKED',technical_reasons='',germline_vcf='',
                            carrier_variant='',target_id=link['target_id']))
    receipt=dict(status='FILES_PRESENT' if check_files else 'LINKED_NOT_SERVER_CHECKED',
                 rows=len(samples),workbook_sha256=sha(workbook),linkage_sha256=sha(linkage),
                 source_sheet='Review',path_column='BAM_path',
                 clinical_review_used_for_exclusion=False,external_controls_added=False)
    return samples,receipt

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--workbook',required=True);p.add_argument('--linkage',required=True)
    p.add_argument('--out',required=True);p.add_argument('--check-files',action='store_true')
    a=p.parse_args();rows,receipt=load_workbook_samples(a.workbook,a.linkage,check_files=a.check_files)
    out=Path(a.out)
    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/'samples.private.tsv',rows);write_json(out/'workbook_input.json',receipt)
    print(json.dumps(receipt))

if __name__=='__main__':main()

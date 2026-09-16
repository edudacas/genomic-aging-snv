import tempfile
import unittest
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape
from aging.common import write_tsv
from aging.workbook_input import load_workbook_samples

class WorkbookInputTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.xlsx=self.root/'study.xlsx';self.linkage=self.root/'linkage.tsv'
        self.bam=self.root/'one.bam';self.bam.write_bytes(b'fixture')
        self.bai=self.root/'one.bam.bai';self.bai.write_bytes(b'index')
        self.row=dict(individual_id='person1',analysis_id='analysis1',fap='123456789012',
                      BAM_path=str(self.bam),assay_id='',study_group='LFS_R337H',study_design_member='yes',
                      genetic_eligibility_status='PENDING_GERMLINE_VAF')
        write_tsv(self.linkage,[dict(individual_id='person1',analysis_id='analysis1',fap='123456789012',
                  bam_path=str(self.bam),bai_path=str(self.bai),bam_bytes=7,assay_id='resolved_assay',
                  workbook_assay_id='',target_id='target1',sample_id='GA_1',person_id='person1')])

    def save(self, rows=None):
        rows=rows or [self.row];headers=list(self.row)
        def column(n):
            result=''
            while n:n,k=divmod(n-1,26);result=chr(65+k)+result
            return result
        xmlrows=[]
        for i,values in enumerate([headers]+[[r.get(h,'') for h in headers] for r in rows],1):
            cells=''.join(f'<c r="{column(j)}{i}" t="inlineStr"><is><t>{escape(str(v))}</t></is></c>' for j,v in enumerate(values,1))
            xmlrows.append(f'<row r="{i}">{cells}</row>')
        with zipfile.ZipFile(self.xlsx,'w') as z:
            z.writestr('xl/workbook.xml','<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Review" sheetId="1" r:id="rId1"/></sheets></workbook>')
            z.writestr('xl/_rels/workbook.xml.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
            z.writestr('xl/worksheets/sheet1.xml','<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'+''.join(xmlrows)+'</sheetData></worksheet>')

    def test_reads_workbook_and_preserves_pending_member(self):
        self.save();rows,receipt=load_workbook_samples(self.xlsx,self.linkage,expected=1,check_files=True)
        self.assertEqual(rows[0]['bam_path'],str(self.bam))
        self.assertEqual(rows[0]['genetic_eligibility_status'],'PENDING_GERMLINE_VAF')
        self.assertEqual(rows[0]['assay_id'],'resolved_assay')
        self.assertEqual(rows[0]['sample_id'],'GA_1')
        self.assertEqual(receipt['rows'],1)
        self.assertNotIn('Name',rows[0])

    def test_rejects_altered_path(self):
        self.row['BAM_path']=str(self.root/'other.bam');self.save()
        with self.assertRaisesRegex(ValueError,'BAM_path'):load_workbook_samples(self.xlsx,self.linkage,expected=1)

    def test_rejects_missing_index(self):
        self.save();self.bai.unlink()
        with self.assertRaisesRegex(ValueError,'absent'):load_workbook_samples(self.xlsx,self.linkage,expected=1,check_files=True)

    def test_rejects_duplicate_members(self):
        self.save([self.row,self.row])
        with self.assertRaisesRegex(ValueError,'duplicate'):load_workbook_samples(self.xlsx,self.linkage,expected=2)

    def test_rejects_changed_identity(self):
        self.row['fap']='999999999999';self.save()
        with self.assertRaisesRegex(ValueError,'identity'):load_workbook_samples(self.xlsx,self.linkage,expected=1)

    def test_rejects_unexpected_denominator(self):
        self.save()
        with self.assertRaisesRegex(ValueError,'1329'):load_workbook_samples(self.xlsx,self.linkage)

"""Minimum published artifacts; shared by producer and read-only exporter."""
QC_OUTPUTS=('qc.json','coverage.intervals.tsv.gz','coverage.compartments.tsv',
            'coverage.summary.json','depth_histogram.tsv','pileups.table',
            'contamination.table','segments.table')


def required_outputs(stage):
    if stage=='qc':return QC_OUTPUTS
    if stage not in ('sage','mutect','lofreq'):raise ValueError('unknown stage')
    outputs=(stage+'.norm.vcf.gz',stage+'.norm.vcf.gz.tbi','commands.log')
    if stage=='sage':outputs+=('sage.vcf.gz','sage.vcf.gz.tbi','sage.arguments.txt','sage.timing.json')
    if stage=='mutect':outputs+=('mutect.raw.vcf.gz','mutect.filtered.vcf.gz','orientation.tar.gz')
    if stage=='lofreq':outputs+=('lofreq.vcf',)
    return outputs


def validate_outputs(directory,stage,receipt=None):
    for name in required_outputs(stage):
        path=directory/name
        if not path.is_file() or path.is_symlink() or path.stat().st_size==0:
            raise ValueError('required '+stage+' output missing or empty: '+name)
        if receipt is not None and name not in receipt['outputs']:
            raise ValueError('required output missing from receipt: '+name)

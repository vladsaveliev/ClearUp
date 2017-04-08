import ngs_utils.logger as log
from fingerprinting.model import Project


def get_snp_record(snps_dict, snp_a, snp_b, snp_index):
    seq_a, seq_b = snp_a.usercall or snp_a.genotype, snp_b.usercall or snp_b.genotype
    seq_a, seq_b = seq_a.replace('N', ''), seq_b.replace('N', '')
    snp_record = {'index': snp_index,
                  'chrom': snp_a.location.chrom,
                  'pos': snp_a.location.pos,
                  'rsid': snp_a.location.rsid,
                  'gene': snp_a.location.gene,
                  'depthA': snp_a.depth,
                  'depthB': snp_b.depth,
                  'snpA': seq_a,
                  'snpB': seq_b,
                  'usercallA': 'usercall' if snp_a.usercall else '',
                  'usercallB': 'usercall' if snp_b.usercall else '',
                  'penalty': 0,
                  'class': ''}
    if not seq_a or not seq_b:
        snps_dict['snp_missing'] += 1
        snp_record['class'] += ' nocall'
        return snp_record
    if seq_a == seq_b:
        snps_dict['matches'] += 1
        snp_record['class'] += ' match'
        snp_record['penalty'] = 0
    elif seq_a[0] == seq_b[0] or seq_a[1] == seq_b[1]:
        snps_dict['het_matches'] += 1
        snp_record['class'] += ' het_match'
        snp_record['penalty'] = 1
    else:
        snps_dict['mismatches'] += 1
        snp_record['class'] += ' mistmatch'
        snp_record['penalty'] = 2
    return snp_record


def get_sample_by_name(sample_name, project_name):
    project = Project.query.get(project_name)
    if not project:
        log.err('Project ' + project_name + ' not found in database')
        return None
    sample = project.samples.get(sample_name)
    if not sample:
        log.err('Sample ' + sample_name + ' not found in ' + project_name)
        return None
    return sample

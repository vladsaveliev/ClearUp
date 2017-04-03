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

    is_match = seq_a == seq_b
    is_homozygous = False
    if len(seq_a) == 1:
        is_homozygous = True
    elif seq_a[0] == seq_a[1] and seq_b[0] == seq_b[1]:
        is_homozygous = True
    if is_match and is_homozygous:
        snps_dict['hom_matches'] += 1
    elif is_match and not is_homozygous:
        snps_dict['het_matches'] += 1
    elif not is_match and is_homozygous:
        snps_dict['hom_mismatches'] += 1
        snp_record['class'] += ' hom_mismatch'
        snp_record['penalty'] = 3
    elif not is_match and not is_homozygous:
        snps_dict['het_mismatches'] += 1
        snp_record['class'] += ' het_mismatch'
        snp_record['penalty'] = 3
    return snp_record


def get_sample_by_name(sample_name, project_name):
    project = Project.query.filter_by(name=project_name).first()
    if not project:
        log.err('Project ' + project_name + ' not found in database')
        return None
    sample = project.samples.filter_by(name=sample_name).first()
    if not project:
        log.err('Sample ' + sample_name + ' not found in ' + project_name)
        return None
    return sample

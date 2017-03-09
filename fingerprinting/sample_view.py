import json
import os
import re
import subprocess
import time

from collections import defaultdict
from os.path import abspath, join, dirname, splitext, basename
from flask import Flask, render_template, send_from_directory, abort, redirect, send_file
from flask import Response, request

from ngs_utils import logger as log
from ngs_utils.file_utils import safe_mkdir, file_transaction, can_reuse

from fingerprinting import config
from fingerprinting.model import Project, db, Sample, PairedSample
from fingerprinting.lookups import get_snp_record, get_sample_by_name


def render_closest_comparison_page(run_id, sample_id, selected_idx=None):
    sample = Sample.query.filter_by(id=sample_id).first()
    if not sample:
        log.err('Sample ' + sample_id + ' not found in ' + run_id)
        abort(404)
    sample_name = sample.name
    matching_sample_id = sample.paired_samples.filter_by(run_id=run_id).first().sample_id
    if not matching_sample_id:
        log.err('No matched sample for ' + sample_name)
        abort(404)
    matching_sample = Sample.query.filter_by(id=matching_sample_id).first()
    snps_dict = defaultdict(int)
    # snps_dict['total_score'] = 0
    # snps_dict['confidence'] = 'Low'
    snp_tables = []
    snp_records = []
    for snp_a, snp_b in zip(sample.fingerprints, matching_sample.fingerprints):
        snp_records.append(get_snp_record(snps_dict, snp_a, snp_b))
        if snp_a.index % 41 == 0:
            snp_tables.append(snp_records)
            snp_records = []
    snp_tables.append(snp_records)

    bam_fpath_a = '/%s/bamfiles/%s' % (sample.project.name, sample_name + '.bam')
    bam_fpath_b = '/%s/bamfiles/%s' % (matching_sample.project.name, matching_sample.name + '.bam')
    snps_bed = '/snps/snps_bed'
    sample_a = {
        'id': sample.id,
        'name': sample.name,
        'project': sample.project.name,
        'bam': bam_fpath_a,
    }
    sample_b = {
        'id': matching_sample.id,
        'name': matching_sample.name,
        'project': matching_sample.project.name,
        'bam': bam_fpath_b,
    }
    t = render_template(
        'sample.html',
        run_id=run_id,
        genome=sample.project.genome,
        sampleA=sample_a,
        sampleB=sample_b,
        snps_data=snps_dict,
        snp_tables=snp_tables,
        snps_bed=snps_bed,
        selected_idx=selected_idx or "null",
    )
    return t


def send_file_for_igv(fpath):
    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get('Range', None)
    if not range_header:
        return send_file(fpath)

    m = re.search('(\d+)-(\d*)', range_header)
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        log.err(error_msg)
        return error_msg

    size = os.path.getsize(fpath)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset

    with open(fpath, 'rb') as f:
        f.seek(offset)
        data = f.read(length)

    rv = Response(data, 206, mimetype="application/octet-stream", direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(offset, offset + length - 1, size))

    log.info("GET range request: %s-%s %s" % (m.group(1), m.group(2), fpath))
    return rv


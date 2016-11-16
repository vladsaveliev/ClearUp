import json
import os
import re
import subprocess
import time

from collections import defaultdict
from os.path import abspath, join, dirname, splitext, basename
from flask import Flask, render_template, send_from_directory, abort, redirect
from flask import Response, request

from Utils import logger as log
from Utils.file_utils import safe_mkdir, file_transaction, can_reuse

from fingerprinting import config
from fingerprinting.model import Project, db, Sample, PairedSample
from fingerprinting.lookups import get_snp_record, get_sample_by_name


def render_sample_page(run_id, sample_id):
    sample = Sample.query.filter_by(id=sample_id).first()
    if not sample:
        log.err('Sample not found in ' + run_id)
        abort(404)
    sample_name = sample.name
    matching_sample_id = sample.paired_samples.filter_by(run_id=run_id).first().sample_id
    if not matching_sample_id:
        log.err('No matched sample for ' + sample_name)
        abort(404)
    matching_sample = Sample.query.filter_by(id=matching_sample_id).first()
    snps_dict = defaultdict(int)
    snps_dict['total_score'] = 0
    snps_dict['confidence'] = 'Low'
    snp_tables = []
    snp_records = []
    for snp_a, snp_b in zip(sample.fingerprints, matching_sample.fingerprints):
        snp_records.append(get_snp_record(snps_dict, snp_a, snp_b))
        if snp_a.index % 49 == 0:
            snp_tables.append(snp_records)
            snp_records = []
    snp_tables.append(snp_records)

    bam_fpath_a = '/%s/bamfiles/%s' % (sample.project.name, sample_name + '.bam')
    bam_fpath_b = '/%s/bamfiles/%s' % (matching_sample.project.name, matching_sample.name + '.bam')
    sample_a = {
        'id': sample.id,
        'name': sample.name,
        'project': sample.project.name,
        'bam': bam_fpath_a
    }
    sample_b = {
        'id': matching_sample.id,
        'name': matching_sample.name,
        'project': matching_sample.project.name,
        'bam': bam_fpath_b
    }
    t = render_template(
        'sample.html',
        run_id=run_id,
        genome=sample.project.genome,
        sampleA=sample_a,
        sampleB=sample_b,
        snps_data=snps_dict,
        snp_tables=snp_tables
    )
    return t


def send_bam_files(run_id, bam_fname):
    project_work_dirpath = safe_mkdir(join(config.DATA_DIR, run_id))
    bam_fpath = join(project_work_dirpath, bam_fname)
    full_path = abspath(bam_fpath)

    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get('Range', None)
    if not range_header:
        return send_from_directory(project_work_dirpath, bam_fname)

    m = re.search('(\d+)-(\d*)', range_header)
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        log.err(error_msg)
        return error_msg

    size = os.path.getsize(full_path)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset

    with open(full_path, 'rb') as f:
        f.seek(offset)
        data = f.read(length)

    rv = Response(data, 206, mimetype="application/octet-stream", direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(offset, offset + length - 1, size))

    log.info("GET range request: %s-%s %s" % (m.group(1), m.group(2), full_path))
    return rv

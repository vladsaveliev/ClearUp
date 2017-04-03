import os
from os.path import abspath, join, dirname
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from ngs_utils.utils import is_us
from ngs_utils.parallel import ParallelCfg
import az


sys_cfg = az.init_sys_cfg()
# parallel_cfg = ParallelCfg(sys_cfg.get('scheduler'), sys_cfg.get('queue'),
#                            sys_cfg.get('resources'), sys_cfg.get('threads'))
parallel_cfg = ParallelCfg()


DATA_DIR = abspath(join(dirname(__file__), '..', 'data'))


app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:////' + join(DATA_DIR, 'projects.db')
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)


if is_us():
    HOST_IP = '172.18.72.171'
    PORT = 5003
else:
    HOST_IP = 'localhost'
    PORT = 5004


def get_version():
    try:
        from fingerprinting import version
    except ImportError:
        version = None
    else:
        version = version.__version__
    return version


DEPTH_CUTOFF = 5



from os.path import abspath, join, dirname
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from ngs_utils import logger
from ngs_utils.utils import is_us
from ngs_utils.parallel import ParallelCfg


try:
    import az
except ImportError:
    sys_cfg = None
else:
    sys_cfg = az.init_sys_cfg()


if is_us():
    HOST_IP = 'rask.usbod.astrazeneca.net'
    PORT = 5003
    parallel_cfg = ParallelCfg(sys_cfg.get('scheduler'), sys_cfg.get('queue'),
                               sys_cfg.get('resources'), sys_cfg.get('threads'))
    # parallel_cfg = ParallelCfg(threads=20)
else:
    HOST_IP = 'localhost'
    PORT = 5004
    parallel_cfg = ParallelCfg()


DATA_DIR = abspath(join(dirname(__file__), '..', 'data'))


app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:////' + join(DATA_DIR, 'projects.db')
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)


def get_version():
    try:
        from clearup import version
    except ImportError:
        version = None
    else:
        version = version.__version__
    return version


DEPTH_CUTOFF = 5
AF_CUTOFF = 0.30

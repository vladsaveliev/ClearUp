from os.path import abspath, join, dirname
from ngs_utils.utils import is_us


if is_us():
    HOST_IP = '172.18.72.171'
    PORT = 5001
else:
    HOST_IP = 'localhost'
    PORT = 5003

DATA_DIR = abspath(join(dirname(__file__), '..', 'data'))


def get_version():
    try:
        from fingerprinting import version
    except ImportError:
        version = None
    else:
        version = version.__version__
    return version

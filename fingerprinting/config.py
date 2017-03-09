from os.path import abspath, join, dirname
import ngs_utils.logger as log


IS_DEBUG = log.is_debug = True


if log.is_local():
    HOST_IP = 'localhost'
    PORT = 5003
else:
    HOST_IP = '172.18.72.171'
    PORT = 5001


DATA_DIR = abspath(join(dirname(__file__), '..', 'data'))


def get_version():
    try:
        from fingerprinting import version
    except ImportError:
        version = None
    else:
        version = version.__version__
    return version

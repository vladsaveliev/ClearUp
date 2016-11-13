from os.path import abspath, join, dirname
import Utils.logger as log


IS_DEBUG = log.is_debug = True


if log.is_local:
    HOST_IP = 'localhost'
else:
    HOST_IP = '172.18.72.171'


DATA_DIR = abspath(join(dirname(__file__), '..', 'data'))



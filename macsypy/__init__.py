from time import strftime, localtime

__version__ = 'master-{}'.format(strftime("%Y%m%d", localtime()))

__MACSY_CONF__ = '$MACSYCONF'
__MACSY_DATA__ = '$MACSYDATA'
from time import strftime, localtime

__version__ = 'master-{}'.format(strftime("%Y-%m-%d", localtime()))

__MACSY_CONF__ = '$PREFIXCONF'
__MACSY_DATA__ = '$PREFIXDATA'
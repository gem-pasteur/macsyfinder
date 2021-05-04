#!/usr/bin/env bash

CMD="$1"
shift
ARGS=${@}

case ${CMD} in
	macsyfinder )
		exec /usr/local/bin/macsyfinder ${ARGS} ;;
	macsydata )
		exec /usr/local/bin/macsydata ${ARGS} ;;
  macsyprofile)
    exec /usr/local/bin/macsyprofile ${ARGS} ;;
	* )
		echo -e "command \"${CMD}\" is not supported.\nAvailable commands: macsyfinder | macsydata | macsyprofile"
    exit 127
    ;;
esac

#!/usr/bin/env bash

CMD="$1"
shift
ARGS=${@}

case ${CMD} in
	macsyfinder )
		exec /usr/local/bin/macsyfinder ${ARGS} ;;
	msf_data )
		exec /usr/local/bin/msf_data ${ARGS} ;;
  msf_profile)
    exec /usr/local/bin/msf_profile ${ARGS} ;;
	* )
		echo -e "command \"${CMD}\" is not supported.\nAvailable commands: macsyfinder | msf_data | msf_profile"
    exit 127
    ;;
esac

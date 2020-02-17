#!/bin/bash

#https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
 
if  [ -z ${1+x} ];
then
	echo Specify directory
	exit
else
	echo Building $1
	for DEF in "$1"/*.def; do
		if [ -e ${DEF/def/simg} ]; then
			echo " - Skipping $DEF"
		else
			echo " - Building $DEF"
			set -x
			sudo singularity build "${DEF/def/simg}" "$DEF"
			set +x
		fi
	done
fi

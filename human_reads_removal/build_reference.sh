#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

if [ ! -e "$DIR/nohost/hg19_main_mask_ribo_animal_allplant_allfungus.fa" ]; then
	echo "Download the reference from 'https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit'"
	echo "and place it in nohost/hg19_main_mask_ribo_animal_allplant_allfungus.fa (gunzipped)"
	exit
elif [ ! -d "$DIR/nohost/ref" ] && [ ! -e "$DIR/nohost/ref/genome/1/summary.txt" ]; then
	echo "Building index"
	set -x
	bbmap.sh ref="$DIR/nohost/hg19_main_mask_ribo_animal_allplant_allfungus.fa" path="$DIR/nohost" -Xmx24000m
else
	echo "All done"
fi


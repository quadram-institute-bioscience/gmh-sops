#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
set -euxo pipefail
seqfu='/Users/telatin/git/seqfu/bin/seqfu_mac'
mkdir -p "$DIR"/refs
mkdir -p "$DIR"/reads
mkdir -p "$DIR"/assembly

gunzip "$DIR"/src/*/*.* || true

for SAMPLE in "$DIR"/sample*txt;
do
	ID=$(basename $SAMPLE | cut -f1 -d.)
	perl make_sample.pl "$SAMPLE" "$DIR"/src 2> refs/${ID}.fasta > reads/${ID}.fq
	$seqfu dei -o reads/${ID}_ reads/${ID}.fq
	spades.py --12 reads/${ID}.fq -o "$DIR"/assembly/${ID}/ -t 4 --meta
	rm reads/${ID}.fq
	mv  "$DIR"/assembly/${ID}/contigs.fasta "$DIR"/assembly/${ID}_assembly.fasta

done
rm -rf "$DIR"/assembly/
gzip "$DIR"/src/*/*.* || true
gzip "$DIR"/{reads,refs}/*.*

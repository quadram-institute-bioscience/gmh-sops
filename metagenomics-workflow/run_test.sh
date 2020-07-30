#!/bin/bash 
set -euxo pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

if [[ $USER == 'ubuntu' ]];
then
 OUT="$HOME/volume/_test/"
 PROF='vm'
else
 OUT="$DIR"/_test/"$USER"
 PROF='test'
fi
mkdir -p "$OUT"
mkdir -p "$OUT/tmp"

nextflow run test.nf  -with-timeline "$OUT"/timeline.html \
  -with-report "$OUT"/report.html  \
  -with-dag "$OUT"/dag.png -profile $PROF \
  -work-dir "$OUT"/ \
  --tempdir "$OUT"/tmp \
  --outdir "$OUT"/output -resume

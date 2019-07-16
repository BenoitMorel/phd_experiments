#!/bin/bash

# source: https://github.com/lczech/genesis-gappa-paper/blob/master/plot_to_pdf.sh

for filename in `ls $1/*.svg` ; do
    echo "Converting $filename ..."
    inkscape -D -z --file=${filename} --export-pdf=${filename//svg/pdf}
  done



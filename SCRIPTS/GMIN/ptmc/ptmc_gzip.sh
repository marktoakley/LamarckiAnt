#!/bin/bash
# zip the output files from a PTMC / BSPT simulation

echo "zipping GMIN_out.* files"
gzip GMIN_out.*

echo "zipping dumpstruct.* files"
gzip dumpstruct.*

echo "zipping overlap.* files"
gzip overlap.*

for d in `seq 100`; do
  if [ ! -d $d ]; then
    break
  fi
  echo "zipping $d/Visits.his* files"
  gzip $d/Visits.his*
done

#!/bin/bash
sed -n 's/^Final Quench.*energy=[ ]*\([^ ]*\).*$/\1/p' output | awk "{print sqrt((\$1-($1))**2) <= $2 }"

#!/bin/bash

#set -x
set -u

function fatal () {
    echo $1
    exit 1
}

which datamash || fatal "FATAL: datamash must be installed"

TMPIN=`mktemp`
TMPOUT=`mktemp`

function cleanup {
    rm -f $TMPIN $TMPOUT
}

trap cleanup EXIT

# build executables we want to test
dune build src/gauss_test.exe
dune build src/welford_test.exe

# some wanted distributions MU and SIGMA params
cat <<EOF > $TMPIN
0.0	1.0
-0.1	0.5
-0.5	0.25
1.2	2.0
-3.0	3.125
10.0	5.0
EOF

while read MU SIGMA; do
    printf "%.3f\t%.3f\tREF\n" $MU $SIGMA
    _build/default/src/gauss_test.exe $MU $SIGMA 100000 > $TMPOUT
    _build/default/src/welford_test.exe $TMPOUT
    cat $TMPOUT | datamash mean 1 sstdev 1 | \
        awk '{printf("%.3f\t%.3f\n",$1,$2)}'
done < $TMPIN

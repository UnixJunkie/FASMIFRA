#!/bin/bash

set -x
set -u

TMP=`mktemp`

function cleanup {
    rm -f $TMP
}

trap cleanup EXIT

dune build src/gauss_test.exe

_build/default/src/gauss_test.exe 0.0 1.0 10000 > $TMP
cat $TMP | datamash mean 1 sstdev 1

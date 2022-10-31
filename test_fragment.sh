#!/bin/bash

set -u

# get some molecules
unxz -k data/CHEMBL_100k.smi.xz
head -20000 data/CHEMBL_100k.smi > chembl_20000.smi

# tag cut bonds
time ./bin/fragment_SA.py -i chembl_20000.smi -o chembl_20000_frags.smi

# create fragments dictionary and encode molecules
for n in `echo 64 128 256 512 1024 2048 4096 8192 16384`; do
    echo "bits: "$n
    time ./bin/fragment_SA_dict.py -i chembl_20000_frags.smi \
         -o chembl_20000_dict$n.smi -e chembl_20000_b$n.txt -n $n 2> err_$n.log
done

# representable molecules Vs num_bits
wc -l chembl_20000_b*.txt

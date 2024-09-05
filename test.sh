#!/usr/bin/env bash

set -u
set -x # DEBUG

# compile the exe
# if people have executed ./install.sh; then the software is already installed
# and should be in their PATH
#make 

# get a molecular training set
xzcat data/CHEMBL_100k.smi.xz | head -1000 > chembl_1k.smi

# fragment them
./bin/fasmifra_fragment.py -i chembl_1k.smi -o chembl_1k_frags.smi

# generate molecules from those fragments
# -f: overwrite fragments cache, if any
./fasmifra -f -n 1000 -i chembl_1k_frags.smi -o gen_1k.smi

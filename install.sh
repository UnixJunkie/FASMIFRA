#!/usr/bin/env bash

set -u

sudo apt install -y opam
opam init -y # user-install of ocaml-4.13.1 as of May 16th 2023
pip3 install rdkit
eval `opam config env`
opam install --fake conf-rdkit
opam install -y fasmifra

# test install was successful
which fasmifra_fragment.py
which fasmifra

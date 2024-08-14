#!/usr/bin/env python3
#
# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, Tokyo University,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# Read a SMILES fragments file, then output on stdout
# the corresponding fragment dictionary.
# Output lines are of the form:
# '^<smiles:str>\t<cano_smiles:str>\t<frag_id:int>\n'
# i.e. each SMILES fragment has its canonical SMILES computed
# and each different canonical SMILES is assigned a new
# fragment identifier.

import rdkit, sys, typing
from rdkit import Chem

input_fn = sys.argv[1]

def dict_has_key(d, k) -> bool:
    return (d.get(k) != None)

with open(input_fn) as input:
    smi2cansmi = {}
    cansmi2id = {}
    count = 0
    for line in input.readlines():
        smi = line.strip()
        # only one SMILES per input line allowed
        assert(len(smi.split()) == 1)
        mol = Chem.MolFromSmiles(smi)
        cano_smi = Chem.MolToSmiles(mol)
        if dict_has_key(smi2cansmi, smi):
            # check rdkit canonical SMILES are stable
            # (they are supposed to)
            assert(smi2cansmi[smi] == cano_smi)
        else:
            # insert new binding
            smi2cansmi[smi] = cano_smi
        if not dict_has_key(cansmi2id, cano_smi):
            cansmi2id[cano_smi] = len(cansmi2id)
        # how many frags were read in
        count += 1
    # output created dictionary to stdout
    # mean and stddev columns will be used by Thompson sampling
    # they are all NaNs initially
    print('#smi\tcano_smi\tid\tmean\tstddev') # format header line
    for smi, cano_smi in smi2cansmi.items():
        frag_id = cansmi2id[cano_smi]
        print('%s\t%s\t%d\tnan\tnan' % (smi, cano_smi, frag_id))
    # user feedback
    print('%s: %d SMILES; %d unique canoSMILES' % \
          (input_fn, count, len(cansmi2id)),
          file=sys.stderr)

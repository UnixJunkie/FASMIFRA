#!/usr/bin/env python3

import rdkit, sys, typing
from rdkit import Chem

input_fn = sys.argv[1]

def dict_has_key(d, k) -> bool:
    return (d.get(k) != None)

with open(input_fn) as input:
    smi2cansmi = {}
    cansmi2id = {}
    for line in input.readlines():
        smi = line.strip()
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
    # output created dictionary to stdout
    for smi, cano_smi in smi2cansmi.items():
        frag_id = cansmi2id[cano_smi]
        print('%s\t%s\t%d' % (smi, cano_smi, frag_id))

#!/usr/bin/env python3
#
# Copyright (C) 2022, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# Create fragments dictionary for molecules fragmented by fragment_SA.py

import argparse
import random
import rdkit
import sys
import time

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import RWMol
from rdkit.Chem.AtomPairs import Pairs

def RobustSmilesSupplier(filename):
    with open(filename) as f:
        for line in f:
            smile, name = line.strip().split("\t") # enforce TAB-separated
            yield (smile, name)

def snd(a_b):
    a, b = a_b
    return b

def key_values(dico):
    return list(zip(dico.keys(), dico.values()))

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "Create fragments dictionary")
    # can also output as bitstrings the input molecules
    parser.add_argument("-i", metavar = "input.smi", dest = "smi_in_fn",
                        help = "input molecules")
    parser.add_argument("-o", metavar = "dico.txt", dest = "dico_out_fn",
                        help = "fragments output dictionary")
    parser.add_argument("-e", metavar = "encoded.smi", dest = "smi_out_fn",
                        help = "output encoded molecules")
    parser.add_argument("--seed", dest = "seed", default = -1,
                        type = int, help = "RNG seed")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    smi_in_fn = args.smi_in_fn
    dico_out_fn = args.dico_out_fn
    smi_out_fn = args.smi_out_fn
    rng_seed = args.seed
    randomize = True
    if rng_seed != -1:
        # only if the user asked for it, we make experiments repeatable
        random.seed(rng_seed)
    output = open(args.dico_out_fn, 'w')
    count = 0
    errors = 0
    # fragments indexing ------------------------------------------------------
    max_frags = 1000 # FBR: put on the CLI; default should be -1 (all frags seen)
    frags_count = {} # how many time each canonical fragment was seen
    for smi, name in RobustSmilesSupplier(smi_in_fn):
        # really cut the bonds
        mol = Chem.MolFromSmiles(smi)
        count += 1        
        if not mol:
            errors += 1
            print("cannot parse: %s" % smi, file=sys.stderr)
            continue
        rw_mol = None
        rw_mol = Chem.RWMol(mol)
        for b in mol.GetBonds():
            left = b.GetBeginAtom()
            right = b.GetEndAtom()
            if left.GetAtomMapNum() > 0 and right.GetAtomMapNum() > 0:
                # two artificially introduced atoms flagging a cut bond
                left_i = left.GetIdx()
                right_i = right.GetIdx()
                rw_mol.RemoveBond(left_i, right_i)
        cut_mol = rw_mol.GetMol()
        cut_mol_smi = Chem.MolToSmiles(cut_mol)
        # print("cut_mol: %s" % cut_mol_smi)
        # cano SMILES for each frag
        for frag_smi in cut_mol_smi.split('.'):
            frag_mol = Chem.MolFromSmiles(frag_smi)
            frag_cano_smi = Chem.MolToSmiles(frag_mol)
            try:
                frags_count[frag_cano_smi] += 1
            except KeyError:
                frags_count[frag_cano_smi] = 1
    # dico and count
    kvs = key_values(frags_count)
    # decr. count sort
    kvs.sort(key=snd, reverse=True)
    for i, (cano_smi, count) in enumerate(kvs):
        print("%s\t%d\t%d" % (cano_smi, i, count), file=output)
    print("seen_frags: %d" % len(kvs), file=sys.stderr)
    after = time.time()
    dt = after - before
    print("read %d molecules at %.2f Hz; %d errors" %
          (count, count / dt, errors), file=sys.stderr)
    output.close()

# FBR: remove the dummy atoms? could correct the kekul errors
#      - will require updating FASMIFRA


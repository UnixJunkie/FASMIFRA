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

def canonicalize(smi):
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol, canonical=True)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "Create fragments dictionary")
    # can also output as bitstrings the input molecules
    parser.add_argument("-i", metavar = "input.smi", dest = "smi_in_fn",
                        help = "input molecules")
    parser.add_argument("-o", metavar = "dico.txt", dest = "dico_out_fn",
                        help = "fragments output dictionary")
    parser.add_argument("-e", metavar = "encoded.txt", dest = "bits_out_fn",
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
    bits_out_fn = args.bits_out_fn
    rng_seed = args.seed
    randomize = True
    if rng_seed != -1:
        # only if the user asked for it, we make experiments repeatable
        random.seed(rng_seed)
    count = 0
    errors = 0
    # fragments indexing ------------------------------------------------------
    max_frags = 1000 # FBR: put on the CLI; default should be -1 (all frags seen)
    frags_count = {} # how many time each canonical fragment was seen
    # 1) build the fragments dictionary ---------------------------------------
    cut_mols = []
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
        fragments = cut_mol_smi.split('.')
        cano_frags = list(map(canonicalize, fragments))
        cut_mols.append((cano_frags, name))
        for frag_cano_smi in cano_frags:
            try:
                frags_count[frag_cano_smi] += 1
            except KeyError:
                frags_count[frag_cano_smi] = 1
    # 2) store fragments dictionary on disk -----------------------------------
    # dico and count
    kvs = key_values(frags_count)
    # decr. count sort
    kvs.sort(key=snd, reverse=True)
    frag_to_index = {}
    with open(dico_out_fn, 'w') as dico_out:
        print("#frag_cano_smi\tindex\tcount", file=dico_out) # header
        for i, (cano_smi, count) in enumerate(kvs):
            frag_to_index[cano_smi] = i
            print("%s\t%d\t%d" % (cano_smi, i, count), file=dico_out)
    dico_size = len(kvs)
    print("seen_frags: %d" % dico_size, file=sys.stderr)
    # 3) encode input molecules -----------------------------------------------
    with open(bits_out_fn, 'w') as bits_out:
        for cano_frags, name in cut_mols:
            # a molecule can only be encoded if it is made of distinct fragments
            # only information loss is how fragments are connected; however, there
            # might be only very few possibilities
            frags_set = set()
            bitstring = list("0" * dico_size)
            num_frags = len(cano_frags)
            for frag_cano_smi in cano_frags:
                frag_i = frag_to_index[frag_cano_smi]
                frags_set.add(frag_i)
                bitstring[frag_i] = "1"
            if num_frags == len(frags_set):
                print("%s\t%s" % ("".join(bitstring), name), file=bits_out)
            else:
                print("duplicated frags in %s" % name, file=sys.stderr)
    # post-processing ---------------------------------------------------------
    after = time.time()
    dt = after - before
    print("read %d molecules at %.2f Hz; %d errors" %
          (count, count / dt, errors), file=sys.stderr)

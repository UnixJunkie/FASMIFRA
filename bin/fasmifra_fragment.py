#!/usr/bin/env python3
#
# Copyright (C) 2021, Francois Berenger
# Tsuda laboratory, Tokyo University,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# Ultra-fast generator of only valid molecules using (Deep)SMILES fragments

import argparse
import fasmifra_common as common
import random
import rdkit
import sys
import time

from fasmifra_common import RobustSmilesMolSupplier
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import RWMol

def get_name(mol):
    return mol.GetProp("name")

def set_name(mol, name):
    mol.SetProp("name", name)

def index_for_atom_type(atom_types_dict, atom_type):
    try:
        return atom_types_dict[atom_type]
    except KeyError:
        # want indexes to start at 1; so the isotope number is
        # always explicit in the fragments SMILES output
        v = len(atom_types_dict) + 1
        atom_types_dict[atom_type] = v
        return v

def get_cut_bonds(frag_weight, mol):
    cuttable_bonds = common.find_cuttable_bonds(mol)
    cut_bonds_indexes = [b.GetIdx() for b in cuttable_bonds]
    total_weight = Descriptors.MolWt(mol)
    nb_frags = round(total_weight / frag_weight)
    max_cuts = min(len(cut_bonds_indexes), nb_frags - 1)
    # print("mol %s; cut %d bonds" % (mol.GetProp("name"), max_cuts),
    #       file=sys.stderr)
    random.shuffle(cut_bonds_indexes)
    return cut_bonds_indexes[0:max_cuts]

def random_reorder_atoms(randomize, mol):
    if randomize:
        rand_order = list(range(mol.GetNumAtoms()))
        random.shuffle(rand_order)
        rand_mol = Chem.RenumberAtoms(mol, newOrder=rand_order)
        set_name(rand_mol, get_name(mol))
        return rand_mol
    else:
        return mol

def tag_cut_bonds(frag_weight, randomize, atom_types_dict, input_mol):
    name = get_name(input_mol)
    mol = random_reorder_atoms(randomize, input_mol)
    to_cut = get_cut_bonds(frag_weight, mol)
    if len(to_cut) == 0:
        # molecule too small: not fragmented
        # still, we output it so that input and output SMILES files can be
        # visualized side-by-side
        # also, this makes the molecular generator reading this fragments file
        # able to generate any of the molecules that were fragmented (good point)
        smi = Chem.MolToSmiles(mol)
        return (smi, name)
    else:
        rw_mol = Chem.RWMol(mol)
        # atom type bonds end points in the input molecule
        # (atom typing does not work on a RWMol)
        # tag and replace all cut bonds in the RWMol
        for i in to_cut:
            b_i = mol.GetBondWithIdx(i)
            a_j = b_i.GetBeginAtom()
            j = a_j.GetIdx()
            a_j_t = common.type_atom(mol.GetAtomWithIdx(j))
            a_k = b_i.GetEndAtom()
            k = a_k.GetIdx()
            a_k_t = common.type_atom(mol.GetAtomWithIdx(k))
            rw_mol.RemoveBond(j, k)
            # we will replace this bond by -['*][''*]-
            # ' and '' will be isotope numbers encoding the bond
            # start and end atom types
            start_o = rw_mol.AddAtom(Chem.Atom(0)) # wildcard atom has atomic number 0
            start_a = rw_mol.GetAtomWithIdx(start_o)
            # encode attached atom's type with an isotope number
            start_a.SetIsotope(index_for_atom_type(atom_types_dict, a_j_t))
            end_o = rw_mol.AddAtom(Chem.Atom(0)) # wildcard atom
            stop_o = rw_mol.GetAtomWithIdx(end_o)
            # encode attached atom's type
            stop_o.SetIsotope(index_for_atom_type(atom_types_dict, a_k_t))
            rw_mol.AddBond(j, start_o, Chem.BondType.SINGLE)
            rw_mol.AddBond(start_o, end_o, Chem.BondType.SINGLE)
            rw_mol.AddBond(end_o, k, Chem.BondType.SINGLE)
        new_mol = rw_mol.GetMol()
        # forbid the generated SMILES to start w/ an unspecified atom
        assert(new_mol.GetAtomWithIdx(0).GetAtomicNum() != 0);
        smi = Chem.MolToSmiles(new_mol, canonical=(not randomize))
        return (smi, name)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "fragment molecules (tag cut bonds)")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "fragments output file")
    parser.add_argument("--seed", dest = "seed", default = -1,
                        type = int, help = "RNG seed")
    parser.add_argument("-n", dest = "nb_passes", default = 1,
                        type = int, help = "number of fragmentation passes")
    # 150 Da: D. Rognan's suggested max fragment weight
    parser.add_argument("-w", dest = "frag_weight", default = 150.0,
                        type = float, help = "fragment weight (default=150Da)")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    nb_passes = args.nb_passes
    rng_seed = args.seed
    frag_weight = args.frag_weight
    randomize = True
    if rng_seed != -1:
        # only if the user asked for it, we make experiments repeatable
        random.seed(rng_seed)
    output = open(args.output_fn, 'w')
    count = 0
    # fragmenting ---------------------------------------------------------
    mol_supplier = RobustSmilesMolSupplier(input_fn)
    seen_types_dict = {}
    for name, mol in mol_supplier:
        for i in range(nb_passes):
            tagged_bonds_smi, parent_name = tag_cut_bonds(frag_weight, randomize, seen_types_dict, mol)
            if i == 0:
                print("%s\t%s" %
                      (tagged_bonds_smi, parent_name), file=output)
            else:
                print("%s\t%sp%d" %
                      (tagged_bonds_smi, parent_name, i), file=output)
        count += 1
    print("seen_types: %d" % len(seen_types_dict), file=sys.stderr)
    after = time.time()
    dt = after - before
    print("read %d molecules at %.2f mol/s" %
          (count, count / dt), file=sys.stderr)
    output.close()
